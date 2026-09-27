#!/usr/bin/env python3
r"""
Step 055: Leg-permutation null + realistic non-factorizable scintillation null
==============================================================================

Addresses the open falsification gap that all step_008 simulations implement
scalar-factorizable screens (X_ij = a_i a_j*, for which the closure phase
psi = arg(X_01 X_12 X_02*) = 0 is a mathematical identity), while the measured
observable is extracted from a single realization of the secondary spectrum of
an intensity (chi^2) field whose cross-term phases are generically NOT
factorizable under standard scintillation physics.

Two controls are implemented:

(A) Leg-permutation null (the control deferred in Section 4.4).
    For each viable epoch the three cross-term leg phases
    (phi_01, phi_12, phi_02) are recomputed through the identical
    step_002 -> step_003 chain (secondary spectrum, Hough arcs, arclet
    detection, cross-screen triplet selection, sub-pixel cross-term
    measurement, 3x3 complex-patch phase extraction).  The closure phase
    psi = phi_01 + phi_12 - phi_02 is then reconstructed after permuting one
    leg's phases across the triplets of the same epoch.  This destroys any
    genuine per-loop phase association while exactly preserving the marginal
    distribution of every leg, the estimator chain, and the epoch weighting.
    If the measured concentration is a true closure (loop) property, the
    permuted circular concentration collapses toward the uniform value; if it
    survives, the signal is a per-leg marginal property (screen structure).
    A stronger cross-epoch variant draws each leg from a different epoch,
    removing all realization-specific association.

(B) End-to-end non-factorizable null (squared-Gaussian scintillation field).
    Synthetic dynamic spectra I(nu,t) = |E(nu,t)|^2 are synthesized on the
    identical (n_time, n_freq) grid with matched sampling, where E is a
    complex Gaussian field whose conjugate spectrum rho(tau, f_D) concentrates
    on the measured J0437 scintillation arcs: two parabolas
    tau = eta_s f_D^2 with empirically measured curvatures, f_D extents, and
    band thickness, plus additive white noise matched to the observed flux
    statistics.  For a chi^2 field the intensity cross-terms are convolutions
    of the field spectrum -- not factorizable products -- so this model carries
    the standard-physics channel the audit identified.  The synthetic spectra
    are pushed through step_002 and step_003 unchanged (arc detection, arclet
    selection, cross-screen triplet requirement, cross-term measurement, phase
    extraction, epoch aggregation).  The resulting psi distribution is the
    pipeline-level null for a realistic persistently-structured screen.

Outputs: results/step_055_leg_permutation_null.json
"""
import sys, json, contextlib, io
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

from scripts.utils.logger import TEPLogger, print_status, set_step_logger
from scripts.utils.parallel_workers import worker_count

logger = TEPLogger("step_055_leg_permutation_null")

PROC_DIR = PROJECT_ROOT / "data" / "processed" / "j0437"
RESULTS = PROJECT_ROOT / "results"

N_PERM = 400          # leg-permutation replicates per leg choice
N_SYNTH = 1500        # synthetic epochs for the end-to-end null
RNG_SEED = 20260924
MIN_TRIPLETS_EPOCH = 5  # matches step_003 aggregation cut


# ---------------------------------------------------------------------------
# Shared epoch machinery (identical calls to step_002 / step_003 functions)
# ---------------------------------------------------------------------------

def _epoch_secondary(dyn, dt_s, freq_MHz):
    """Run the unchanged step_002 chain on one dynamic spectrum."""
    from steps.step_002_secondary_spectra import (
        compute_secondary_spectrum,
        detect_arcs,
        detect_arclets,
        detect_peaks_direct,
    )

    ss = compute_secondary_spectrum(dyn, dt_s, freq_MHz)
    arcs = detect_arcs(ss["secondary"], ss["tau_us"], ss["fD_mHz"])
    arclets = detect_arclets(
        ss["secondary"], ss["tau_us"], ss["fD_mHz"], arcs, min_snr=1.5
    )
    if len(arclets) < 3:
        dp = detect_peaks_direct(
            ss["secondary"], ss["tau_us"], ss["fD_mHz"], min_snr=1.5
        )
        if len(dp) > len(arclets):
            arclets = dp
    eta1 = arcs[0]["eta"] if len(arcs) >= 1 else 0.0
    eta2 = arcs[1]["eta"] if len(arcs) >= 2 else 0.0
    return ss, arclets, eta1, eta2


def _leg_phase(S_complex, fD_mHz, tau_us, ct):
    """Replicates step_003.compute_closure.get_phase_at_peak exactly."""
    f_idx = np.argmin(np.abs(fD_mHz - ct.fD_meas))
    t_idx = np.argmin(np.abs(tau_us - ct.tau_meas))
    f_start = max(0, f_idx - 1)
    f_end = min(len(fD_mHz), f_idx + 2)
    t_start = max(0, t_idx - 1)
    t_end = min(len(tau_us), t_idx + 2)
    region = S_complex[f_start:f_end, t_start:t_end]
    return float(np.angle(np.mean(region)))


def _epoch_legs(dynspec, dt_s, freq_MHz, mjd):
    """Return per-triplet leg phases and triplet SNRs for one epoch, or None."""
    from steps.step_003_closure_delays_final import (
        find_best_triplets,
        calculate_velocity_vector,
    )

    ss, arclets, eta1, eta2 = _epoch_secondary(dynspec, dt_s, freq_MHz)
    if len(arclets) < 3:
        return None
    v_eff = calculate_velocity_vector(mjd, pulsar_name="J0437-4715", verbose=False)
    res = find_best_triplets(
        ss["secondary"],
        ss["secondary_complex"],
        ss["tau_us"],
        ss["fD_mHz"],
        arclets,
        v_eff,
        eta1=eta1,
        eta2=eta2,
        mjd=mjd,
        nu_ref_mhz=float(np.mean(freq_MHz)),
        pulsar_name="J0437-4715",
    )
    if not res:
        return None
    legs, snrs, deltas = [], [], []
    for r in res:
        phi01 = _leg_phase(ss["secondary_complex"], ss["fD_mHz"], ss["tau_us"], r.tau_01)
        phi12 = _leg_phase(ss["secondary_complex"], ss["fD_mHz"], ss["tau_us"], r.tau_12)
        phi02 = _leg_phase(ss["secondary_complex"], ss["fD_mHz"], ss["tau_us"], r.tau_02)
        legs.append((phi01, phi12, phi02))
        snrs.append(r.snr)
        deltas.append(r.geometric_delta_us)
    return {
        "legs": np.asarray(legs),          # (n_trip, 3)
        "snr": np.asarray(snrs),           # triplet SNR
        "geo_delta": np.asarray(deltas),   # for epoch weight (unchanged by permutation)
        "mjd": mjd,
    }


def _circ_mean(angles, weights):
    z = np.sum(weights * np.exp(1j * angles))
    return float(np.angle(z)), float(np.clip(np.abs(z) / np.sum(weights), 0, 1))


def _aggregate(epoch_list):
    """Reproduce step_003's epoch-level aggregation from leg-phase records.

    epoch_list: list of dicts with 'legs' (n_trip,3), 'snr', 'geo_delta'.
    Returns psi_mean, rbar, n_epochs (epochs with >= MIN_TRIPLETS_EPOCH triplets).
    """
    psi_epoch, w_epoch = [], []
    for e in epoch_list:
        legs, snr, gd = e["legs"], e["snr"], e["geo_delta"]
        if len(legs) < MIN_TRIPLETS_EPOCH:
            continue
        psi = (legs[:, 0] + legs[:, 1] - legs[:, 2] + np.pi) % (2 * np.pi) - np.pi
        m, _ = _circ_mean(psi, np.square(np.maximum(snr, 1e-6)))
        sem = np.std(gd, ddof=1) / np.sqrt(len(gd)) if len(gd) > 1 else 1e-3
        psi_epoch.append(m)
        w_epoch.append(1.0 / (sem**2 + 0.001**2))
    if not psi_epoch:
        return None
    psi_epoch = np.asarray(psi_epoch)
    w_epoch = np.asarray(w_epoch)
    mean, rbar = _circ_mean(psi_epoch, w_epoch)
    return {"psi_mean": mean, "rbar": rbar, "n_epochs": len(psi_epoch),
            "epoch_psi": psi_epoch, "epoch_w": w_epoch}


def _permute_legs(epoch_list, rng, leg_mode="01"):
    """One permutation replicate: permute one leg's phases across triplets.

    leg_mode '01', '12', '02' permutes that leg within each epoch;
    'any' picks a random leg per epoch; 'xepoch' draws the leg-01 phase
    from a random triplet of a different epoch.
    """
    legs_by_epoch = [e["legs"] for e in epoch_list]
    pooled01 = np.concatenate([L[:, 0] for L in legs_by_epoch])
    out = []
    for e in epoch_list:
        legs = e["legs"].copy()
        n = len(legs)
        if leg_mode == "xepoch":
            idx = rng.integers(0, len(pooled01), size=n)
            legs[:, 0] = pooled01[idx]
        else:
            col = {"01": 0, "12": 1, "02": 2}.get(leg_mode)
            if col is None:
                col = int(rng.integers(0, 3))
            legs[:, col] = legs[rng.permutation(n), col]
        out.append({"legs": legs, "snr": e["snr"], "geo_delta": e["geo_delta"]})
    return out


def _worker_real(npz_name):
    import logging
    logging.disable(logging.CRITICAL)
    p = PROC_DIR / npz_name
    if not p.exists():
        return None
    d = np.load(p)
    with contextlib.redirect_stdout(io.StringIO()):
        return _epoch_legs(
            d["dynspec"], float(d["dt_s"]), d["freq_MHz"], float(d["mjd_start"])
        )


# ---------------------------------------------------------------------------
# Synthetic squared-Gaussian scintillation field
# ---------------------------------------------------------------------------

def _synth_dynspec(rng, n_t, n_f, dt_s, freq_MHz, eta_us_mhz2, sig_fD_mHz,
                   w_tau_us, flux_mean, flux_std, n_arcs=2, eta_ratio=0.7,
                   asymmetric=False):
    """Synthesize I(nu,t) = |E|^2 with conjugate spectrum on measured arcs.

    rho(tau, fD) is concentrated on parabolas tau = eta_s fD^2 (one per
    screen), with Gaussian fD extent sig_fD and arc thickness w_tau.  The
    field E = IFFT2(sqrt(rho) * g), g iid CN(0,1), is complex Gaussian, so
    I = |E|^2 is a chi^2 speckle field whose secondary-spectrum cross-terms
    are convolutions -- the non-factorizable standard-physics channel.
    """
    dnu_Hz = abs(np.median(np.diff(freq_MHz))) * 1e6
    u_t = np.fft.fftfreq(n_t, d=dt_s)          # Hz -> fD axis
    u_f = np.fft.fftfreq(n_f, d=dnu_Hz)        # s  -> tau axis
    fD_Hz = u_t[:, None]
    tau_s_ax = u_f[None, :]

    rho = np.zeros((n_t, n_f))
    etas = [eta_us_mhz2, eta_us_mhz2 * eta_ratio][:n_arcs]
    for s, eta in enumerate(etas):
        eta_si = eta * 10.0                    # us/mHz^2 -> s^3
        sig = sig_fD_mHz * 1e-3                # mHz -> Hz
        w = w_tau_us * 1e-6                    # us -> s
        arc = eta_si * fD_Hz ** 2
        amp = 1.0
        if asymmetric:
            amp = np.where(fD_Hz >= 0, 1.0, 0.55)  # uneven arms
        rho += amp * np.exp(-0.5 * (fD_Hz / sig) ** 2) * np.exp(
            -0.5 * ((tau_s_ax - arc) / w) ** 2
        )
    # Bright unscattered/core image at the origin: core x scattered-image
    # pairs trace the parabola tau = eta fD^2 in the secondary spectrum,
    # which is how the observed arc is generated physically.
    w0 = w_tau_us * 1e-6
    sig0 = sig_fD_mHz * 1e-3
    rho += 6.0 * np.exp(-0.5 * (fD_Hz / (0.3 * sig0)) ** 2) * np.exp(
        -0.5 * (tau_s_ax / w0) ** 2
    )
    rho = np.fft.ifftshift(rho)
    g = (rng.standard_normal((n_t, n_f)) + 1j * rng.standard_normal((n_t, n_f))) / np.sqrt(2)
    E = np.fft.ifft2(np.sqrt(rho) * g)
    E /= np.abs(E).mean()
    I = np.abs(E) ** 2
    # radiometer noise matched to observed scatter about the speckle mean
    noise = np.sqrt(max(flux_std**2 - np.var(I), 0.01 * flux_mean**2))
    I = I * flux_mean + rng.standard_normal((n_t, n_f)) * noise
    return np.clip(I, 0.0, None).astype(np.float32)


def _worker_synth(args):
    import logging
    logging.disable(logging.CRITICAL)
    (seed, n_t, n_f, dt_s, freq_MHz, eta, eratio, sigf, wtau,
     fmean, fstd, mjd, asym) = args
    rng = np.random.default_rng(seed)
    dyn = _synth_dynspec(rng, n_t, n_f, dt_s, freq_MHz, eta, sigf, wtau,
                       fmean, fstd, eta_ratio=eratio, asymmetric=asym)
    with contextlib.redirect_stdout(io.StringIO()):
        return _epoch_legs(dyn, dt_s, freq_MHz, mjd)


# ---------------------------------------------------------------------------
# Calibration of synthetic parameters from real data
# ---------------------------------------------------------------------------

def _calibrate(epoch_names, rng, n_cal=40):
    """Measure dynspec stats, arc curvatures and arclet fD extent."""
    from steps.step_002_secondary_spectra import (
        compute_secondary_spectrum, detect_arcs, detect_arclets, detect_peaks_direct,
    )
    flux_m, flux_s, dts, fref = [], [], [], None
    eta1s, eta2s, fDmax = [], [], []
    for nm in rng.choice(epoch_names, min(n_cal, len(epoch_names)), replace=False):
        d = np.load(PROC_DIR / nm)
        dyn, freq, dt = d["dynspec"], d["freq_MHz"], float(d["dt_s"])
        flux_m.append(float(dyn.mean())); flux_s.append(float(dyn.std())); dts.append(dt)
        fref = freq
        with contextlib.redirect_stdout(io.StringIO()):
            ss = compute_secondary_spectrum(dyn, dt, freq)
            arcs = detect_arcs(ss["secondary"], ss["tau_us"], ss["fD_mHz"])
            arclets = detect_arclets(ss["secondary"], ss["tau_us"], ss["fD_mHz"], arcs, min_snr=1.5)
            if len(arclets) < 3:
                dp = detect_peaks_direct(ss["secondary"], ss["tau_us"], ss["fD_mHz"], min_snr=1.5)
                if len(dp) > len(arclets):
                    arclets = dp
        if arcs:
            eta1s.append(arcs[0]["eta"])
            if len(arcs) > 1:
                eta2s.append(arcs[1]["eta"])
        if len(arclets):
            fDmax.append(float(np.abs(arclets[:, 1]).max()))
    return {
        "flux_mean": float(np.median(flux_m)),
        "flux_std": float(np.median(flux_s)),
        "dt_s": float(np.median(dts)),
        "freq_MHz": fref,
        "eta1": float(np.median(eta1s)) if eta1s else 0.0074,
        "eta2": float(np.median(eta2s)) if eta2s else 0.0067,
        "fD_arclet_max_mHz": float(np.median(fDmax)) if fDmax else 1.3,
    }


def main():
    set_step_logger(logger)
    rng = np.random.default_rng(RNG_SEED)

    stored = json.load(open(RESULTS / "step_003_closure_final_per_epoch_j0437.json"))
    epoch_names = [e["epoch"].replace("_secondary", ".npz") for e in stored]
    epoch_names = [n for n in epoch_names if (PROC_DIR / n).exists()]
    print_status(f"Recomputing leg phases for {len(epoch_names)} stored viable epochs", "INFO")

    workers = worker_count(role="cpu_bound", reserve=2)
    legs_cache = RESULTS / "step_055_leg_phases_cache.npz"
    epochs = []
    if legs_cache.exists():
        z = np.load(legs_cache, allow_pickle=False)
        legs_all = z["legs"]; snr_all = z["snr"]; gd_all = z["geo_delta"]
        mjd_all = z["mjd"]; idx_all = z["idx"]
        for e_i, (i, j) in enumerate(idx_all):
            epochs.append({"legs": legs_all[i:j], "snr": snr_all[i:j],
                           "geo_delta": gd_all[i:j], "mjd": float(mjd_all[e_i])})
        print_status(f"Loaded leg-phase cache: {len(epochs)} epochs", "INFO")
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for out in ex.map(_worker_real, epoch_names, chunksize=16):
                if out is not None:
                    epochs.append(out)
        legs_cat = np.concatenate([e["legs"] for e in epochs])
        snr_cat = np.concatenate([e["snr"] for e in epochs])
        gd_cat = np.concatenate([e["geo_delta"] for e in epochs])
        idx = np.zeros((len(epochs), 2), dtype=int)
        c = 0
        for e_i, e in enumerate(epochs):
            idx[e_i] = (c, c + len(e["legs"]))
            c += len(e["legs"])
        np.savez_compressed(
            legs_cache, legs=legs_cat, snr=snr_cat, geo_delta=gd_cat,
            mjd=np.array([e["mjd"] for e in epochs]), idx=idx,
        )
        print_status(f"Recomputed {len(epochs)} viable epochs (cached)", "INFO")

    agg_real = _aggregate(epochs)
    print_status(
        f"Real aggregation: psi={agg_real['psi_mean']:.3f} rad, "
        f"Rbar={agg_real['rbar']:.4f}, n={agg_real['n_epochs']}",
        "RESULT",
    )

    # --- (A) leg-permutation nulls -------------------------------------------
    perm = {}
    for mode in ["01", "12", "02", "any", "xepoch"]:
        rbars, means = [], []
        for rep in range(N_PERM):
            a = _aggregate(_permute_legs(epochs, rng, mode))
            if a:
                rbars.append(a["rbar"]); means.append(a["psi_mean"])
        rbars = np.asarray(rbars); means = np.asarray(means)
        perm[mode] = {
            "rbar_mean": float(rbars.mean()),
            "rbar_std": float(rbars.std()),
            "psi_mean_of_means": float(np.angle(np.mean(np.exp(1j * means)))),
            "p_ge_real": float((np.sum(rbars >= agg_real["rbar"]) + 1) / (len(rbars) + 1)),
        }
        print_status(
            f"perm[{mode}]: Rbar={rbars.mean():.4f}+-{rbars.std():.4f} "
            f"(real {agg_real['rbar']:.4f}, p={perm[mode]['p_ge_real']:.4f})",
            "RESULT",
        )

    # --- (B) end-to-end non-factorizable null ---------------------------------
    cal = _calibrate(epoch_names, rng)
    print_status(f"Synthetic calibration: {cal}", "DATA")

    mjds = np.array([e["mjd"] for e in epochs])
    ref = np.load(PROC_DIR / epoch_names[0])
    n_t, n_f = ref["dynspec"].shape
    sig_fD = max(cal["fD_arclet_max_mHz"] / 2.0, 0.5)  # half the observed span
    w_tau = 0.003                                     # ~1.5 delay bins
    eta_ratio = cal["eta2"] / cal["eta1"] if cal["eta1"] > 0 else 0.7
    tasks = []
    for i in range(N_SYNTH):
        mjd = float(rng.choice(mjds))
        tasks.append((RNG_SEED + 1000 + i, n_t, n_f, cal["dt_s"], cal["freq_MHz"],
                      cal["eta1"], eta_ratio, sig_fD, w_tau,
                      cal["flux_mean"], cal["flux_std"], mjd, False))
    # asymmetric-screen variant: uneven arc arms break parity
    for i in range(N_SYNTH):
        mjd = float(rng.choice(mjds))
        tasks.append((RNG_SEED + 5000 + i, n_t, n_f, cal["dt_s"], cal["freq_MHz"],
                      cal["eta1"], eta_ratio, sig_fD, w_tau,
                      cal["flux_mean"], cal["flux_std"], mjd, True))

    synth = []
    synth_flags = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for out, t in zip(ex.map(_worker_synth, tasks, chunksize=8), tasks):
            if out is not None:
                synth.append(out)
                synth_flags.append(bool(t[-1]))
    synth_flags = np.asarray(synth_flags)
    agg_syn = _aggregate(synth)
    if agg_syn:
        print_status(
            f"Synthetic null: n={agg_syn['n_epochs']} viable, "
            f"psi_mean={agg_syn['psi_mean']:.3f} rad, Rbar={agg_syn['rbar']:.4f}",
            "RESULT",
        )

    # Matched-n comparison: resample the real epoch-psi distribution down to
    # the synthetic viable count so Rbar is compared at equal sample size.
    matched = None
    if agg_syn and agg_syn["n_epochs"] > 0:
        n_s = agg_syn["n_epochs"]
        rbar_sub = []
        ep_real, w_real = agg_real["epoch_psi"], agg_real["epoch_w"]
        for rep in range(400):
            idx = rng.integers(0, len(ep_real), size=n_s)
            _, rb = _circ_mean(ep_real[idx], w_real[idx])
            rbar_sub.append(rb)
        rbar_sub = np.asarray(rbar_sub)
        matched = {
            "n_synth_viable": int(n_s),
            "rbar_real_subsample_mean": float(rbar_sub.mean()),
            "rbar_real_subsample_std": float(rbar_sub.std()),
            "p_synth_ge_real": float(
                (np.sum(rbar_sub <= agg_syn["rbar"]) + 1) / (len(rbar_sub) + 1)
            ),
        }
        # epoch-mean direction: fraction of synthetic epoch means near +0.984
        ep_syn = agg_syn["epoch_psi"]
        near = np.abs((ep_syn - agg_real["psi_mean"] + np.pi) % (2 * np.pi) - np.pi) < 0.5
        matched["frac_synth_epoch_within_0p5rad_of_real"] = float(np.mean(near))
        # uniformity of synthetic epoch-mean directions (Rayleigh on means)
        z = np.abs(np.sum(np.exp(1j * ep_syn))) ** 2 / len(ep_syn)
        matched["synth_epoch_mean_rayleigh_z"] = float(z)
        matched["synth_epoch_mean_rayleigh_p"] = float(np.exp(-z) * (
            1 + (2 * z - z**2) / (4 * len(ep_syn))
        ))

    # split symmetric / asymmetric variants
    split = {}
    for flag, name in [(False, "symmetric"), (True, "asymmetric")]:
        sub = [s for s, fl in zip(synth, synth_flags) if fl == flag]
        a = _aggregate(sub)
        # unweighted direction persistence of the epoch means
        dir_rbar = None
        if a:
            ep = a["epoch_psi"]
            dir_rbar = float(np.abs(np.mean(np.exp(1j * ep))))
        split[name] = {
            "n_viable": a["n_epochs"] if a else 0,
            "psi_mean": a["psi_mean"] if a else None,
            "rbar": a["rbar"] if a else None,
            "epoch_dir_rbar_unweighted": dir_rbar,
        }

    # leg-permutation applied to the synthetic set (control validation)
    rbars_syn = []
    for rep in range(50):
        a = _aggregate(_permute_legs(synth, rng, "any"))
        if a:
            rbars_syn.append(a["rbar"])
    rbars_syn = np.asarray(rbars_syn)

    out = {
        "control_A_leg_permutation": {
            "description": (
                "Per-triplet leg phases (phi_01, phi_12, phi_02) recomputed "
                "through the unchanged step_002->step_003 chain on every "
                "stored viable epoch; psi = phi_01 + phi_12 - phi_02 "
                "reconstructed after permuting one leg across triplets "
                "within each epoch ('xepoch': leg drawn from a different "
                "epoch). Destroys loop-level association while preserving "
                "every leg marginal, the estimator chain, and epoch weights."
            ),
            "n_perm_replicates": N_PERM,
            "real": {
                "psi_mean": agg_real["psi_mean"],
                "rbar": agg_real["rbar"],
                "n_epochs": agg_real["n_epochs"],
            },
            "permuted": perm,
        },
        "control_B_nonfactorizable_screen": {
            "description": (
                "I(nu,t)=|E|^2 with E a complex Gaussian field whose "
                "conjugate spectrum concentrates on the measured J0437 arc "
                "curvatures (eta1=%.4f, eta2=%.4f us/mHz^2, arclet |fD|max "
                "%.2f mHz modelled as a Gaussian of sigma=%.2f mHz, arc "
                "thickness %.4f us) plus white noise matched to the "
                "observed flux statistics. Pushed through step_002/step_003 "
                "unchanged. Intensity cross-terms of a chi^2 field are "
                "non-factorizable by construction -- the standard-physics "
                "channel the factorizable step_008 suite cannot reach."
                % (cal["eta1"], cal["eta2"], cal["fD_arclet_max_mHz"],
                   sig_fD, w_tau)
            ),
            "n_synth_tasks": len(tasks),
            "n_viable": agg_syn["n_epochs"] if agg_syn else 0,
            "psi_mean": agg_syn["psi_mean"] if agg_syn else None,
            "rbar": agg_syn["rbar"] if agg_syn else None,
            "epoch_psi_std": float(np.std(agg_syn["epoch_psi"])) if agg_syn else None,
            "epoch_dir_rbar_unweighted": (
                float(np.abs(np.mean(np.exp(1j * agg_syn["epoch_psi"]))))
                if agg_syn else None
            ),
            "matched_n_comparison": matched,
            "variant_split": split,
            "permuted_rbar_mean": float(rbars_syn.mean()) if len(rbars_syn) else None,
            "permuted_rbar_std": float(rbars_syn.std()) if len(rbars_syn) else None,
        },
        "interpretation": (
            "If the permuted Rbar collapses well below the real value, the "
            "concentration is a loop-level association, not a leg-marginal "
            "artifact. If the realistic screen null's Rbar/mean is far below "
            "the observed psi=0.984 rad, Rbar=0.308 concentration, standard "
            "non-factorizable screen structure does not reproduce it."
        ),
    }
    path = RESULTS / "step_055_leg_permutation_null.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print_status(f"saved: {path}", "INFO")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
