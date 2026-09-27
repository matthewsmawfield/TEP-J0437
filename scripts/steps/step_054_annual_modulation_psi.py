#!/usr/bin/env python3
"""
Step 054: Annual Modulation Test on Phase Closure and Signed Delays
====================================================================

Tests the unscreened kinematic channel identified in Section 2.2.2:
Earth's ~30 km/s orbital motion modulates the effective velocity u in
u·∇φ with a known period (1 yr) and known phase (ecliptic longitude).
Unlike the pulsar's orbital reflex (screened by the white-dwarf
companion), the observer-side velocity acts on the clock congruence
end of the comparison links and is not subject to companion screening.

Procedure:
  1. Reconstruct the per-epoch SNR^2-weighted circular-mean psi and the
     per-epoch signed geometric-delay mean exactly as in Step 003
     (>= 5 triplets per epoch).
  2. Least-squares fit a fixed-period sinusoid (T = 365.25 d) to each
     per-epoch series: y = c + A_c cos(2πt/T) + A_s sin(2πt/T).
  3. Report the amplitude, its standard error from the residual-scaled
     covariance, the F-test p-value against the constant model, and the
     peak phase as a fraction of the year.

Outputs:
  results/step_054_annual_modulation_psi.json
"""

import json
import numpy as np
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent.parent / "results"
T_YEAR_DAYS = 365.25


def load_epoch_series(per_epoch_path):
    mjd, psi, signed_ns = [], [], []
    with open(per_epoch_path) as f:
        epochs = json.load(f)
    for e in epochs:
        triplets = e.get("triplets", [])
        if len(triplets) < 5:
            continue
        p = np.array([t["phase_closure_rad"] for t in triplets], dtype=float)
        s = np.array([t["snr"] for t in triplets], dtype=float)
        w = np.square(np.maximum(s, 1e-6))
        z = np.sum(w * np.exp(1j * p))
        if not np.isfinite(z) or z == 0:
            continue
        g = np.array([t["geometric_delta_us"] for t in triplets], dtype=float)
        mjd.append(float(e["mjd"]))
        psi.append(float(np.angle(z)))
        signed_ns.append(float(np.mean(g)) * 1e3)
    return np.array(mjd), np.array(psi), np.array(signed_ns)


def sinusoid_fit(mjd, y):
    """Fixed-period least-squares sinusoid; returns amp, se_amp, F, p, peak phase."""
    n = len(y)
    th = 2.0 * np.pi * (mjd - mjd.min()) / T_YEAR_DAYS
    X = np.column_stack([np.ones(n), np.cos(th), np.sin(th)])
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - 3
    rss = float(np.sum(resid**2))
    rss0 = float(np.sum((y - y.mean()) ** 2))
    s2 = rss / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    amp = float(np.hypot(beta[1], beta[2]))
    se_amp = float(np.sqrt(cov[1, 1] + cov[2, 2]))
    F = ((rss0 - rss) / 2.0) / (rss / dof)
    try:
        from scipy import stats as st
        p = float(1 - st.f.cdf(F, 2, dof))
    except Exception:
        # Fallback: chi-square approximation of the F statistic
        p = float(np.exp(-0.5 * F * 2.0) * (1 + 0.5 * F))
    peak_frac = float((np.angle(beta[1] + 1j * beta[2]) % (2 * np.pi)) / (2 * np.pi))
    return {
        "n": int(n),
        "const": float(beta[0]),
        "amp": amp,
        "amp_se": se_amp,
        "amp_over_se": amp / se_amp if se_amp > 0 else None,
        "beta_cos": float(beta[1]),
        "beta_sin": float(beta[2]),
        "se_cos": float(se[1]),
        "se_sin": float(se[2]),
        "F": float(F),
        "p_value": p,
        "peak_fraction_of_year": peak_frac,
        "amp_95_upper": amp + 1.96 * se_amp,
    }


def main():
    out = {
        "step": "step_054_annual_modulation_psi",
        "period_days": T_YEAR_DAYS,
        "method": (
            "Per-epoch SNR^2-weighted circular-mean psi and per-epoch signed "
            "geometric-delay mean reconstructed as in Step 003 (>=5 triplets). "
            "Fixed-period (365.25 d) least-squares sinusoid on unwrapped psi "
            "about the global circular mean, and on the signed delay series. "
            "Amplitude SE from residual-scaled design-matrix covariance; "
            "F-test against the constant model."
        ),
        "results": {},
    }
    for tag in ["j0437", "j1603"]:
        path = RESULTS_DIR / f"step_003_closure_final_per_epoch_{tag}.json"
        if not path.exists():
            continue
        mjd, psi, signed_ns = load_epoch_series(path)
        zbar = np.mean(np.exp(1j * psi))
        psi_u = psi - np.angle(zbar)
        psi_u = (psi_u + np.pi) % (2 * np.pi) - np.pi
        out["results"][tag] = {
            "n_epochs": int(len(mjd)),
            "mjd_min": float(mjd.min()),
            "mjd_max": float(mjd.max()),
            "span_years": float((mjd.max() - mjd.min()) / 365.25),
            "psi_annual": sinusoid_fit(mjd, psi_u),
            "signed_delay_annual": sinusoid_fit(mjd, signed_ns),
            "psi_circular_mean_rad": float(np.angle(zbar)),
        }
    with open(RESULTS_DIR / "step_054_annual_modulation_psi.json", "w") as f:
        json.dump(out, f, indent=2)
    for tag, r in out["results"].items():
        pa = r["psi_annual"]
        print(
            f"{tag}: n={r['n_epochs']} span={r['span_years']:.1f} yr; "
            f"psi annual amp={pa['amp']:.4f}+-{pa['amp_se']:.4f} rad "
            f"(F={pa['F']:.2f}, p={pa['p_value']:.3f}); "
            f"signed delay amp={r['signed_delay_annual']['amp']:.3f} ns "
            f"(p={r['signed_delay_annual']['p_value']:.3f})"
        )


if __name__ == "__main__":
    main()
