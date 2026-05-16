"""
Binary pulsar reflex transverse velocity in an equatorial tangent basis.

Maps the pulsar's Keplerian reflex motion to (v_alpha_star, v_delta) in km/s,
using the same tangent convention as 4.74 * mu_mas_yr * D_pc / 1000 for proper
motion (mu_alpha* = mu_alpha cos delta, mu_delta).

Orbital rotation / Thiele–Innes P, Q basis follows the standard visual-binary
decomposition (e.g. Murray & Dermott 1999; exoplanet orbit literature).

References (J0437-4715):
    Reardon et al. 2024, ApJL 971 L19 (arXiv:2407.07132) — PPTA-DR3 timing:
    Pb, T0, x, e, omega, i, Omega used for consistency with that solution.
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np


def pq_unit_vectors(omega_deg: float, inc_deg: float, kom_deg: float) -> Tuple[np.ndarray, np.ndarray]:
    """Thiele–Innes orthonormal basis (P, Q) in equatorial Cartesian coordinates."""
    om = math.radians(omega_deg)
    i = math.radians(inc_deg)
    Om = math.radians(kom_deg)
    si, ci = math.sin(i), math.cos(i)
    so, co = math.sin(om), math.cos(om)
    SO, CO = math.sin(Om), math.cos(Om)
    Px = CO * co - SO * so * ci
    Py = SO * co + CO * so * ci
    Pz = so * si
    Qx = -CO * so - SO * co * ci
    Qy = -SO * so + CO * co * ci
    Qz = co * si
    P = np.array([Px, Py, Pz], dtype=float)
    Q = np.array([Qx, Qy, Qz], dtype=float)
    return P, Q


def _solve_kepler(M: float, e: float, max_iter: int = 40, tol: float = 1e-14) -> float:
    """Eccentric anomaly from mean anomaly (radians)."""
    if e <= 0.0:
        return M
    E = M if abs(M) < math.pi else math.copysign(math.pi, M)
    for _ in range(max_iter):
        f = E - e * math.sin(E) - M
        fp = 1.0 - e * math.cos(E)
        dE = f / fp
        E -= dE
        if abs(dE) < tol:
            break
    return E


def equatorial_xyz_to_pm_kms(v_xyz: np.ndarray, ra_rad: float, dec_rad: float) -> np.ndarray:
    """
    Project equatorial Cartesian velocity (km/s) onto (mu_alpha*, mu_delta) directions.

    Matches step_003 convention: first component is along increasing RA with cos(delta)
    folded into mu_alpha*; second is along increasing declination.
    """
    ca, sa = math.cos(ra_rad), math.sin(ra_rad)
    cd, sd = math.cos(dec_rad), math.sin(dec_rad)
    # Eastward unit vector (increasing RA) at the pulsar: (-sin alpha, cos alpha, 0)
    ex = np.array([-sa, ca, 0.0], dtype=float)
    # Northward unit vector on the sky
    en = np.array([-ca * sd, -sa * sd, cd], dtype=float)
    v_ra_star = float(np.dot(v_xyz, ex))
    v_dec = float(np.dot(v_xyz, en))
    return np.array([v_ra_star, v_dec], dtype=float)


def reflex_binary_transverse_velocity_kms(
    mjd: float,
    pb_days: float,
    t0_mjd: float,
    x_lt_s: float,
    e: float,
    omega_deg: float,
    inc_deg: float,
    kom_deg: float,
    c_light_km_s: float,
    ra_rad: float,
    dec_rad: float,
) -> np.ndarray:
    """
    Transverse reflex orbital velocity (km/s) in (v_mu_alpha*, v_mu_delta) basis.

    x_lt_s is the DD/T2 projected pulsar semi-major axis (a sin i)/c in light-seconds.
    """
    if pb_days <= 0.0 or x_lt_s <= 0.0:
        raise ValueError("Binary period and projected semi-major axis must be positive.")

    n_orb = 2.0 * math.pi / (pb_days * 86400.0)  # rad/s
    phase = ((mjd - t0_mjd) / pb_days) % 1.0
    M = 2.0 * math.pi * phase
    E = _solve_kepler(M, e)
    sin_i = math.sin(math.radians(inc_deg))
    if sin_i <= 1e-9:
        raise ValueError("Inclination too small: sin(i) ~ 0 is degenerate for projected orbit.")

    # Physical semi-major axis of the pulsar reflex orbit [km]
    a_km = (x_lt_s * c_light_km_s) / sin_i

    ce, se = math.cos(E), math.sin(E)
    den = 1.0 - e * ce
    if abs(den) < 1e-18:
        raise ValueError("Degenerate Kepler denominator 1 - e cos E.")

    sqrt_1me2 = math.sqrt(max(0.0, 1.0 - e * e))
    vx_orb = -a_km * n_orb * se / den
    vy_orb = a_km * n_orb * sqrt_1me2 * ce / den

    P, Q = pq_unit_vectors(omega_deg, inc_deg, kom_deg)
    v_xyz = vx_orb * P + vy_orb * Q
    return equatorial_xyz_to_pm_kms(v_xyz, ra_rad, dec_rad)
