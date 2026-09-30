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
        "step_055_prior": prior_x,
        "note": (
            "step_055 reported p=0.0025, the resolution floor of 400 reps. "
            f"With {len(rbars)} replicates the empirical bound is "
            f"p < {p_emp:.1e}; the observed Rbar sits {z_par:.1f} null-std "
            "above the xepoch null mean (Gaussian-tail p ~ "
            f"{p_par:.1e})."
        ),
    }
    out_path = RESULTS / "step_059_permutation_floor_audit.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out["permuted_xepoch"], indent=2))


if __name__ == "__main__":
    sys.exit(main())
