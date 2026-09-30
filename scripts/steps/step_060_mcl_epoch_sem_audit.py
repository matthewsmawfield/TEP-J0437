#!/usr/bin/env python3
"""
TEP-J0437 audit diagnostic (step_060): epoch-level SEM audit for the
delay-domain closure magnitude M_cl and its noise-subtracted excess.

step_003 reports H_sem_ns = 1/sqrt(sum w_i), an inverse-variance
measurement-noise SEM. The project's stated inference unit is the epoch.
This audit recomputes the standard error of the identical weighted
estimators by resampling epochs:

  - SEM of the raw IVW magnitude mean M_cl
  - SEM of the paired noise-subtracted excess M_cl - E[M_cl]. The folded-
    normal floor E[M_cl] is itself an IVW mean over the same epochs with
    the same weights (per-epoch MAD scale), so each resample recomputes
    BOTH terms and the SEM is taken of their difference. Within-epoch
    corr(m_i, f_i) ~ 0.9 makes this pairing essential: resampling only the
    numerator while holding the floor fixed double-charges epoch-dispersion
    variance the floor already absorbs (the defect corrected relative to
    the superseded step_058 audit, which reported excess_t against the
    raw-mean SEM).
  - Floor-convention robustness: the same audit repeated with a std-based
    floor (sigma = std(d) instead of 1.4826*MAD) bounds the model-choice
    sensitivity of the excess, which no resampling SEM captures.

Resampling modes:

  - i.i.d. epoch bootstrap (independent epoch means)
  - contiguous-block bootstrap in MJD order (moving-block resampling,
    ceil(n/blk) blocks truncated to n per resample), which additionally
    absorbs residual epoch-to-epoch correlation (screen-state /
    observing-campaign clustering)

The block scan doubles as a specification check: for J0437 the epoch-mean
ACF remains ~0.15-0.2 at lag 200, so the raw-mean SEM does not reach a
strict plateau before moving-block resampling breaks down (nb <~ 4
blocks); the paired-excess ACF is near-flat, so its blocked SEM
stabilizes and the excess significance is robust.

Inputs:
  - results/step_003_closure_final_per_epoch_{j0437,j1603}.json
  - results/step_003_closure_final_summary_{j0437,j1603}.json

Outputs:
  - results/step_060_mcl_epoch_sem_audit.json
"""

import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

PACKAGE_ROOT = Path(__file__).resolve().parents[1].parent
if not (PACKAGE_ROOT / "results").exists():
    PACKAGE_ROOT = Path(__file__).resolve().parents[2]

N_BOOT = 10000
RNG_SEED = 20251015

PULSARS = {
    "J0437-4715": {
        "per_epoch": "step_003_closure_final_per_epoch_j0437.json",
        "summary": "step_003_closure_final_summary_j0437.json",
        "blocks": (15, 30, 60, 120, 240),
    },
    "J1603-7202": {
        "per_epoch": "step_003_closure_final_per_epoch_j1603.json",
        "summary": "step_003_closure_final_summary_j1603.json",
        "blocks": (5, 10, 15, 30),
    },
}


def epoch_stats(e):
    """Per-epoch (m, w, f_mad, f_std): mean |d|, IVW weight, floors."""
    d = np.array([t["geometric_delta_us"] for t in e["triplets"]], dtype=float)
    m = float(np.mean(np.abs(d)))
    sd = float(np.std(d, ddof=1)) if len(d) > 1 else 1e-3
    sem = sd / np.sqrt(len(d))
    w = 1.0 / (sem**2 + 0.001**2)  # matches step_003 (1 ns floor, us units)
    med = float(np.median(d))
    mad = float(np.median(np.abs(d - med)))
    sig_mad = 1.4826 * mad
    if not np.isfinite(sig_mad) or sig_mad <= 0.0:
        sig_mad = sd
    sig_std = sd if np.isfinite(sd) and sd > 0.0 else 1e-3
    rt = np.sqrt(2.0 / np.pi)
    return m, w, sig_mad * rt, sig_std * rt


def main():
    rng = np.random.default_rng(RNG_SEED)
    B = N_BOOT
    results = {}

    for pulsar, cfg in PULSARS.items():
        ep = json.loads((PACKAGE_ROOT / "results" / cfg["per_epoch"]).read_text())
        summary = json.loads((PACKAGE_ROOT / "results" / cfg["summary"]).read_text())

        ep5 = [e for e in ep if len(e["triplets"]) >= 5]
        n = len(ep5)
        stats_ = np.array([epoch_stats(e) for e in ep5])
        ms, ws, fs_mad, fs_std = stats_.T

        order_sorted = np.asarray(np.argsort([e["mjd"] for e in ep5]))
        mjd_sorted = np.array([e["mjd"] for e in ep5])[order_sorted]

        obs = float(np.average(ms, weights=ws))
        floor_mad = float(np.average(fs_mad, weights=ws))
        floor_std = float(np.average(fs_std, weights=ws))
        excess_mad = obs - floor_mad
        excess_std = obs - floor_std

        def draw(idx):
            wm = np.average(ms[idx], weights=ws[idx])
            fm = np.average(fs_mad[idx], weights=ws[idx])
            fst = np.average(fs_std[idx], weights=ws[idx])
            return wm, fm, fst

        def scan(order, block=None):
            mean_b = np.empty(B)
            floor_b = np.empty(B)
            exc_mad_b = np.empty(B)
            exc_std_b = np.empty(B)
            if block is None:
                for b in range(B):
                    idx = order[rng.integers(0, n, n)]
                    mean_b[b], fm, fst = draw(idx)
                    floor_b[b] = fm
                    exc_mad_b[b] = mean_b[b] - fm
                    exc_std_b[b] = mean_b[b] - fst
            else:
                nb = int(np.ceil(n / block))
                for b in range(B):
                    starts = rng.integers(0, n - block + 1, nb)
                    sel = np.concatenate([order[s : s + block] for s in starts])[:n]
                    mean_b[b], fm, fst = draw(sel)
                    floor_b[b] = fm
                    exc_mad_b[b] = mean_b[b] - fm
                    exc_std_b[b] = mean_b[b] - fst
            return (
                float(mean_b.std(ddof=1) * 1e3),
                float(floor_b.std(ddof=1) * 1e3),
                float(exc_mad_b.std(ddof=1) * 1e3),
                float(exc_std_b.std(ddof=1) * 1e3),
            )

        sem_iid_mean, sem_iid_floor, sem_iid_exc, sem_iid_exc_std = scan(np.arange(n))

        block_mean, block_floor, block_exc, block_exc_std, block_span = {}, {}, {}, {}, {}
        for blk in cfg["blocks"]:
            sm, sf, sx, sx_std = scan(order_sorted, blk)
            block_mean[str(blk)] = sm
            block_floor[str(blk)] = sf
            block_exc[str(blk)] = sx
            block_exc_std[str(blk)] = sx_std
            spans = mjd_sorted[blk:] - mjd_sorted[:-blk]
            block_span[str(blk)] = float(np.median(spans))

        def acf(x, k):
            xc = x - np.mean(x)
            return float(np.mean(xc[:-k] * xc[k:]) / np.mean(xc * xc))

        lags = [1, 5, 15, 30, 60, 120, 200]
        lags = [k for k in lags if k < n // 3]
        m_ord = ms[order_sorted]
        d_ord = (ms - fs_mad)[order_sorted]
        acf_mean = {str(k): acf(m_ord, k) for k in lags}
        acf_excess = {str(k): acf(d_ord, k) for k in lags}

        # Calendar-half paired-excess audit (Jan-Jun vs Jul-Dec): the excess
        # must persist in both halves to exclude a seasonal subsample driver.
        doy = np.array(
            [
                float((datetime(1858, 11, 17) + timedelta(days=float(m))).timetuple().tm_yday)
                for m in [e["mjd"] for e in ep5]
            ]
        )
        halves = {}
        for name, sel_idx in (
            ("jan_jun", np.where(doy <= 181)[0]),
            ("jul_dec", np.where(doy > 181)[0]),
        ):
            nh = len(sel_idx)
            mh, wh, fh = ms[sel_idx], ws[sel_idx], fs_mad[sel_idx]
            obs_h = np.average(mh, weights=wh) - np.average(fh, weights=wh)
            bs = np.empty(B)
            for b in range(B):
                idx = rng.integers(0, nh, nh)
                bs[b] = np.average(mh[idx], weights=wh[idx]) - np.average(
                    fh[idx], weights=wh[idx]
                )
            sem_h = float(bs.std(ddof=1))
            halves[name] = {
                "n_epochs": int(nh),
                "magnitude_ns": float(np.average(mh, weights=wh) * 1e3),
                "floor_ns": float(np.average(fh, weights=wh) * 1e3),
                "excess_ns": float(obs_h * 1e3),
                "excess_sem_iid_ns": sem_h * 1e3,
                "excess_t_iid": float(obs_h / sem_h),
            }

        results[pulsar] = {
            "n_epochs": n,
            "n_bootstrap": B,
            "replicated_H_mean_ns": obs * 1e3,
            "pipeline_H_mean_ns": summary["H_magnitude_ns"],
            "replicated_noise_floor_ns": floor_mad * 1e3,
            "pipeline_noise_floor_ns": summary["H_noise_bias_ns"],
            "excess_mad_floor_ns": excess_mad * 1e3,
            "excess_std_floor_ns": excess_std * 1e3,
            "pipeline_sem_ns_ivw": summary["H_sem_ns"],
            "pipeline_excess_t": summary["H_excess_t_statistic"],
            "corr_epoch_mean_vs_floor": float(np.corrcoef(ms, fs_mad)[0, 1]),
            "epoch_sem_iid_ns": sem_iid_mean,
            "epoch_sem_iid_floor_ns": sem_iid_floor,
            "epoch_sem_iid_excess_ns": sem_iid_exc,
            "epoch_sem_iid_excess_std_floor_ns": sem_iid_exc_std,
            "epoch_sem_block_ns": block_mean,
            "epoch_sem_block_floor_ns": block_floor,
            "epoch_sem_block_excess_ns": block_exc,
            "epoch_sem_block_excess_std_floor_ns": block_exc_std,
            "block_span_days_median": block_span,
            "acf_epoch_mean_index_lag": acf_mean,
            "acf_epoch_excess_index_lag": acf_excess,
            "seasonal_halves": halves,
            "excess_t_iid": excess_mad * 1e3 / sem_iid_exc,
            "excess_t_block": {k: excess_mad * 1e3 / v for k, v in block_exc.items()},
            "excess_t_iid_std_floor": excess_std * 1e3 / sem_iid_exc_std,
            "excess_t_block_std_floor": {
                k: excess_std * 1e3 / v for k, v in block_exc_std.items()
            },
            "mean_t_iid_reference": obs * 1e3 / sem_iid_mean,
        }

    out = {
        "audit": "epoch_level_sem_for_mcl_excess",
        "n_bootstrap": N_BOOT,
        "rng_seed": RNG_SEED,
        "pulsars": results,
        "note": (
            "pipeline t uses 1/sqrt(sum w_i) (measurement noise only). Epoch "
            "resampling gives the raw-mean SEM and, for the excess, the SEM "
            "of the paired difference m_i - f_i (the floor is an IVW mean "
            "over the same epochs, corr ~0.9). Contiguous moving blocks in "
            "MJD order absorb campaign/screen-state correlation; the raw-mean "
            "blocked SEM is a lower envelope (epoch-mean ACF ~0.15-0.2 to lag "
            "200 for J0437; resampling collapses for blocks ~n/3 or larger), "
            "while the paired-excess SEM stabilizes. The std-floor variant "
            "bounds the folded-normal floor-convention sensitivity. "
            "Supersedes the step_058 audit, which divided the excess by the "
            "raw-mean SEM while holding the correlated floor fixed."
        ),
    }

    out_path = PACKAGE_ROOT / "results/step_060_mcl_epoch_sem_audit.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    sys.exit(main())
