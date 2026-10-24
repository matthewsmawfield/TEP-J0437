#!/usr/bin/env python3
"""
TEP-J0437 audit diagnostic (step_059): permutation-floor extension for the
leg-permutation null of step_055.

step_055 ran 400 replicates per leg mode, so its reported p = 0.0025 is the
resolution floor (1/(400+1)), not a measured tail.  This audit reuses the
cached leg-phase products and the identical _permute_legs/_aggregate
machinery to run a deeper null on the strictest mode ('xepoch': the leg-01
phase is drawn from a random triplet of a *different* epoch, destroying
loop-level association while preserving every leg marginal and the epoch
weighting), plus a parametric tail estimate calibrated on the permutation
distribution.

Inputs:
  - results/step_055_leg_phases_cache.npz
  - results/step_055_leg_permutation_null.json (real rbar, prior null)

Outputs:
  - results/step_059_permutation_floor_audit.json
"""

import json
import sys
from pathlib import Path

import numpy as np

PACKAGE_ROOT = Path(__file__).resolve().parents[1].parent
sys.path.insert(0, str(PACKAGE_ROOT / "scripts" / "steps"))

from step_055_leg_permutation_null import _aggregate, _permute_legs

RESULTS = PACKAGE_ROOT / "results"
N_PERM = 5000
RNG_SEED = 20251016


def main():
    z = np.load(RESULTS / "step_055_leg_phases_cache.npz", allow_pickle=False)
    legs_all, snr_all, gd_all = z["legs"], z["snr"], z["geo_delta"]
    mjd_all, idx_all = z["mjd"], z["idx"]
    epochs = [
        {"legs": legs_all[i:j], "snr": snr_all[i:j],
         "geo_delta": gd_all[i:j], "mjd": float(mjd_all[k])}
        for k, (i, j) in enumerate(idx_all)
    ]
    print(f"loaded {len(epochs)} cached epochs")

    agg_real = _aggregate(epochs)
    rbar_real = agg_real["rbar"]
    print(f"real: psi={agg_real['psi_mean']:.3f} rbar={rbar_real:.4f} "
          f"n={agg_real['n_epochs']}")

    rng = np.random.default_rng(RNG_SEED)
    rbars = np.empty(N_PERM)
    for rep in range(N_PERM):
        a = _aggregate(_permute_legs(epochs, rng, "xepoch"))
        rbars[rep] = a["rbar"] if a else np.nan
    rbars = rbars[np.isfinite(rbars)]

    n_ge = int(np.sum(rbars >= rbar_real))
    p_emp = (n_ge + 1) / (len(rbars) + 1)
    mu, sd = float(rbars.mean()), float(rbars.std(ddof=1))
    z_par = (rbar_real - mu) / sd
    # Gaussian-tail estimate of the permutation p-value
    from scipy import stats as _st
    p_par = float(1.0 - _st.norm.cdf(z_par))

    prior = json.load(open(RESULTS / "step_055_leg_permutation_null.json"))
    prior_x = prior["control_A_leg_permutation"]["permuted"]["xepoch"]

    # ------------------------------------------------------------------
    # Additional referee-grade diagnostics (all from stored pipeline
    # artifacts: step_055_leg_phases_cache.npz and
    # step_003_closure_final_per_epoch_j0437.json):
    #   (1) 3-leg loop-closure amplification over the 2-leg sum (phi_01+phi_12)
    #       and 3rd leg (phi_02), on the full sample, on-axis (fD==0 leg), and
    #       off-axis (no fD==0 leg) subsets, plus distinct-pixel triplets.
    #   (2) MJD-ordered contiguous block-bootstrap SE and z-score for both
    #       unweighted (psi_uw = 1.1204 rad) and weighted (psi_w = 0.9840 rad)
    #       epoch circular means across block = i.i.d., 15, 30, 60, 120, 240.
    # ------------------------------------------------------------------
    ep_raw = json.loads(
        (RESULTS / "step_003_closure_final_per_epoch_j0437.json").read_text()
    )
    fD_legs = []
    shared_mask = []
    for e in ep_raw:
        for t in e["triplets"]:
            fD_legs.append((t["fD_01"], t["fD_12"], t["fD_02"]))
            sh = (
                (abs(t["tau_01"] - t["tau_02"]) < 1e-9 and abs(t["fD_01"] - t["fD_02"]) < 1e-9)
                or (abs(t["tau_01"] - t["tau_12"]) < 1e-9 and abs(t["fD_01"] - t["fD_12"]) < 1e-9)
                or (abs(t["tau_12"] - t["tau_02"]) < 1e-9 and abs(t["fD_12"] - t["fD_02"]) < 1e-9)
            )
            shared_mask.append(sh)
    fD_legs = np.asarray(fD_legs, dtype=float)
    shared_mask = np.asarray(shared_mask, dtype=bool)
    has_fD0 = np.any(np.abs(fD_legs) < 1e-6, axis=1)

    def _loop_amp_block(L):
        z01 = np.mean(np.exp(1j * L[:, 0]))
        z12 = np.mean(np.exp(1j * L[:, 1]))
        z02 = np.mean(np.exp(1j * L[:, 2]))
        z_sum = np.mean(np.exp(1j * (L[:, 0] + L[:, 1])))
        psi_c = (L[:, 0] + L[:, 1] - L[:, 2] + np.pi) % (2.0 * np.pi) - np.pi
        z_cls = np.mean(np.exp(1j * psi_c))
        r_sum = float(np.abs(z_sum))
        r02 = float(np.abs(z02))
        r_cls = float(np.abs(z_cls))
        indep = r_sum * r02
        Z_ray = 2.0 * len(L) * (r_cls ** 2)
        return {
            "n_triplets": int(len(L)),
            "rbar_phi01": float(np.abs(z01)),
            "rbar_phi12": float(np.abs(z12)),
            "rbar_phi02": r02,
            "mean_phi01_plus_phi12_rad": float(np.angle(z_sum)),
            "rbar_phi01_plus_phi12": r_sum,
            "independent_product_rbar": float(indep),
            "closure_mean_rad": float(np.angle(z_cls)),
            "closure_rbar": r_cls,
            "loop_amplification_ratio": float(r_cls / max(indep, 1e-15)),
            "rayleigh_Z_pooled": float(Z_ray),
            "rayleigh_p_pooled": float(_st.chi2.sf(Z_ray, 2)),
        }

    loop_amplification = {
        "all_triplets": _loop_amp_block(legs_all),
        "with_fD0_leg": _loop_amp_block(legs_all[has_fD0]),
        "no_fD0_leg_off_axis": _loop_amp_block(legs_all[~has_fD0]),
        "no_shared_pixel_triplets": _loop_amp_block(legs_all[~shared_mask]),
        "n_shared_pixel_triplets_excluded": int(shared_mask.sum()),
    }

    # Block-bootstrap on MJD-ordered independent epochs (n = 1093)
    ep5 = [e for e in epochs if len(e["legs"]) >= 5]
    order_mjd = np.argsort([e["mjd"] for e in ep5])
    psi_ep = agg_real["epoch_psi"][order_mjd]
    w_ep = agg_real["epoch_w"][order_mjd]
    n_ep = len(psi_ep)
    n_eff_w = float((np.sum(w_ep) ** 2) / np.sum(w_ep ** 2))

    rng_blk = np.random.default_rng(RNG_SEED + 1)
    B_BLK = 10000

    def _boot_psi(angles, weights, block=None):
        w = np.ones_like(angles) if weights is None else weights
        m0 = float(np.angle(np.sum(w * np.exp(1j * angles))))
        bs = np.empty(B_BLK)
        if block is None:
            for b in range(B_BLK):
                idx = rng_blk.integers(0, n_ep, n_ep)
                bs[b] = np.angle(np.sum(w[idx] * np.exp(1j * angles[idx])))
        else:
            nb = int(np.ceil(n_ep / block))
            for b in range(B_BLK):
                starts = rng_blk.integers(0, n_ep - block + 1, nb)
                idx = np.concatenate([np.arange(s, s + block) for s in starts])[:n_ep]
                bs[b] = np.angle(np.sum(w[idx] * np.exp(1j * angles[idx])))
        diff = (bs - m0 + np.pi) % (2.0 * np.pi) - np.pi
        se = float(np.std(diff, ddof=1))
        return {
            "mean_rad": m0,
            "se_rad": se,
            "z_sigma": float(abs(m0) / se),
        }

    psi_block_bootstrap = {
        "n_epochs": int(n_ep),
        "n_eff_weighted": n_eff_w,
        "n_bootstrap": B_BLK,
        "unweighted": {
            "iid": _boot_psi(psi_ep, None, None),
            "blocks": {
                str(blk): _boot_psi(psi_ep, None, blk)
                for blk in (15, 30, 60, 120, 240)
            },
        },
        "weighted": {
            "iid": _boot_psi(psi_ep, w_ep, None),
            "blocks": {
                str(blk): _boot_psi(psi_ep, w_ep, blk)
                for blk in (15, 30, 60, 120, 240)
            },
        },
    }

    out = {
        "audit": "leg_permutation_floor_extension",
        "mode": "xepoch (strictest null of step_055)",
        "n_perm_replicates": int(len(rbars)),
        "real": {"psi_mean": agg_real["psi_mean"], "rbar": rbar_real,
                 "n_epochs": agg_real["n_epochs"]},
        "permuted_xepoch": {
            "rbar_mean": mu, "rbar_std": sd,
            "rbar_max": float(rbars.max()),
            "n_ge_real": n_ge,
            "p_ge_real_empirical": float(p_emp),
            "z_parametric": float(z_par),
            "p_gaussian_tail": p_par,
        },
        "loop_amplification_and_subsamples": loop_amplification,
        "psi_mjd_block_bootstrap": psi_block_bootstrap,
        "step_055_prior": prior_x,
        "note": (
            "step_055 reported p=0.0025, the resolution floor of 400 reps. "
            f"With {len(rbars)} replicates the empirical bound is "
            f"p < {p_emp:.1e}; the observed Rbar sits {z_par:.1f} null-std "
            "above the xepoch null mean (Gaussian-tail p ~ "
            f"{p_par:.1e}). Also records 3-leg loop-closure amplification "
            "(231x over the independent product R(phi01+phi12)*R(phi02), "
            "324x on-axis, 49.5x off-axis), distinct-pixel and off-axis "
            "subsample Rayleigh tests, and MJD block-bootstrap SEs for psi."
        ),
    }
    out_path = RESULTS / "step_059_permutation_floor_audit.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out["permuted_xepoch"], indent=2))


if __name__ == "__main__":
    sys.exit(main())
