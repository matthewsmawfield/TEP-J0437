#!/usr/bin/env python3
"""
STEP 057: Canonical-connection holonomy amplitude for the J0437 screen loop
==========================================================================

Purpose
-------
Close the psi -> H_resid measurement-model mapping (Appendix A.4).  The
weighted Phase Closure psi = +0.984 rad, interpreted under the measurement
model psi = omega_eff * H_resid with omega_eff = 2*pi*1380 MHz, implies an
equivalent transport delay H_implied ~ 1.1e-10 s per loop.  This step
evaluates the canonical disformal connection on the *measured* screen-triangle
geometry on the admissible B >= 0 branch and asks whether that amplitude is
reachable.

Connection (Paper 0, Appendix A3; Paper 0 step_15 implementation convention)
----------------------------------------------------------------------------
    sigma_i dx^i = f(x) du ,   f(x) = -b(u) R_H^2 (u.grad phi) / (A^2 N^2 c)  [m]
    H_resid(C) = (1/c) oint_C sigma . dx                                        [s]
    b(u) = B0 * u^2/(1+u^2) * exp(-u^4 / (2 sigma_B^4))   (corpus envelope;
           admissible branch B0 > 0; Paper 28 benchmark uses B0 = +1,
           sigma_B = 1.5; Paper 0 calibrated magnitude is |B0| = 3.2e-3)

Exactness bookkeeping
---------------------
If the prefactor f depends on position only through u, then f grad u = grad G(u)
is exact and the loop integral vanishes identically.  The disformal curvature
condition of Paper 0 Sec.7, d(delta sigma) = -(B/A^2) d(phi_dot/N^2) wedge d phi,
shows that a residual requires the drift/lapse prefactor to vary independently
of u around the loop.  On the ISM screen triangle the lapse is uniform (no deep
wells), so the operative channel is the spatial inhomogeneity of the drift
field phi_dot(x) -- the temporal landscape (Paper 0 Rule 7) -- parameterised
here by the fractional drift variation eta across the loop, with the
misalignment angle between grad phi_dot and grad u entering as sin(theta_m).

Loop geometry (measured, not assumed)
--------------------------------------
Each triplet leg (i,j) links two scattering-screen events separated by the
image-angle difference.  The measured leg delay |tau_leg| maps to a sky-angle
offset theta = sqrt(2 c |tau| / D_eff) and a screen-plane separation
Dx = theta * D_s, with D_eff = s(1-s) D_p = 37.5 pc and D_s = s D_p = 93.78 pc
(Reardon et al. 2021).  The leg-length distribution is taken from the actual
step_003 triplet products (median |tau_leg| ~ 0.018 us -> Dx ~ 9e9 m).

Field model on the screen (bracketed conventions)
-------------------------------------------------
(A) absolute-depth convention: u_screen ~ u_gal = v_c^2/c^2 ~ 6e-7 (Milky Way
    potential depth at the solar circle measured from the cosmological
    ambient); the coherent gradient is the galactic one,
    grad u ~ a_gal/c^2 = v_c^2/(c^2 R_0) ~ 2.4e-27 m^-1.
(B) ambient-excursion convention (the field origin selected by the Paper 0
    Earth-loop consistency argument): u_screen ~ delta u_ISM sourced by the
    screen overdensity itself, bounded at G rho L^2/c^2 ~ 5e-19 for n_e ~ 0.1
    cm^-3 over ~pc -- orders smaller, and reported as the strict reading.

Drift channels:
    (i)  cosmological drift partition: u.grad phi = Pi_bar = H0 on
         matter-hosting field values (W = 1 inside the galactic field);
    (ii) kinematic projection: u.grad phi = v_eff . grad u ~ 2.5e-22 s^-1,
         subdominant; enters the loop integral only through its own spatial
         variation and is carried as a separate term.

Outputs
-------
results/step_057_holonomy_amplitude_closure.json with:
    - implied H from the measured psi (weighted and unweighted variants)
    - H_canonical per (u convention, B0, eta) including the eta = 0 exactness
      control, B0 = 0 control, and traversal-reversal control
    - eta_required and B0_required to reach the implied H
    - terrestrial-isotropy consistency condition on B0_required
    - the delay-channel bound comparison (|H_signed| <= 0.459 ns)
    - J1603-7202 counterpart under identical conventions
"""

import json
import math
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.utils.json_numpy import NpEncoder
from scripts.utils.logger import TEPLogger, print_status, set_step_logger

RESULTS_DIR = PROJECT_ROOT / "results"
LOG_DIR = PROJECT_ROOT / "logs"

logger = TEPLogger("step_057_holonomy_amplitude_closure")

# ---------------------------------------------------------------------------
# Physical constants (CODATA / corpus conventions)
# ---------------------------------------------------------------------------
C_SI = 2.99792458e8          # m/s, exact
G_SI = 6.674e-11             # m^3 kg^-1 s^-2
H0 = 2.27e-18                # s^-1 (corpus drift scale)
R_H = C_SI / H0              # Hubble radius [m]
PC_M = 3.086e16              # m
AU_M = 1.496e11              # m

# Galactic field at the solar circle (absolute-depth convention A)
V_C = 2.33e5                 # m/s local circular speed
R_0 = 8.2e3 * PC_M           # solar galactocentric radius
U_GAL = V_C**2 / C_SI**2                        # ~ 6.0e-7
A_GAL = V_C**2 / R_0                             # ~ 2.1e-10 m/s^2
GRAD_U_GAL = A_GAL / C_SI**2                     # ~ 2.4e-27 m^-1

# Screen overdensity excursion bound (convention B)
N_E = 0.1e6                  # m^-3  (0.1 cm^-3, generous for a diffuse screen)
M_P = 1.673e-27              # kg
L_SCREEN = 1.0 * PC_M        # characteristic screen-cloud scale
DELTA_U_ISM = (4.0 * np.pi / 3.0) * G_SI * N_E * M_P * L_SCREEN**2 / C_SI**2

# J0437-4715 system geometry (Reardon et al. 2021; Deller et al. 2008)
D_P = 156.3 * PC_M           # pulsar distance
S_SCR = 0.6                  # screen fractional distance D_s/D_p
D_S = S_SCR * D_P            # screen distance from observer ~ 93.8 pc
D_EFF = S_SCR * (1.0 - S_SCR) * D_P   # 37.5 pc
V_EFF = 104.38e3             # m/s effective transverse velocity

# J1603-7202 (Walker et al. 2022; DM distance)
D_P2 = 250.0 * PC_M
S_SCR2 = 0.5
D_S2 = S_SCR2 * D_P2
D_EFF2 = S_SCR2 * (1.0 - S_SCR2) * D_P2

# Band centre for the narrowband-equivalent delay conversion
NU_EFF = 1380.0e6            # Hz
OMEGA_EFF = 2.0 * np.pi * NU_EFF

# Terrestrial bookkeeping for the isotropy consistency condition
M_E = 5.972e24
R_E = 6.371e6
M_SUN = 1.989e30
U_EARTH = (G_SI * M_E / (C_SI**2 * R_E)
           + G_SI * M_SUN / (C_SI**2 * AU_M))     # ~1.06e-8 (step_15 value)
SIG_BARE_EARTH = G_SI * M_E / (C_SI**2 * R_E**2)  # bare surface shear ~1.09e-16 m^-1
D_ISOTROPY_BOUND = 1e-18     # resonator/isotropy bound on D = B M^2 Sigma^2


def B_of_u(u, B0, sigma_B=1.5):
    """Corpus field-space envelope, admissible sign B0 > 0."""
    return B0 * u**2 / (1.0 + u**2) * np.exp(-u**4 / (2.0 * sigma_B**4))


def leg_dx_from_tau(tau_us, d_eff=D_EFF, d_s=D_S):
    """Map a measured leg delay |tau| [us] to the screen-plane separation [m]:
    theta = sqrt(2 c |tau| / D_eff),  Dx = theta * D_s."""
    tau = abs(tau_us) * 1e-6
    theta = np.sqrt(2.0 * C_SI * tau / d_eff)
    return theta * d_s


def loop_integral(vertices, u0, grad_u_vec, B0, udot0, eta, eta_dir,
                  udot_kin=0.0, npts=2000):
    """Numerically evaluate oint f(x) grad u . dx over a closed polygon.

    u(x)    = u0 + grad_u . x          (linear ramp across the loop)
    udot(x) = udot0 * (1 + eta * eta_dir . x / L_leg) + udot_kin
    f(x)    = b(u(x)) R_H^2 udot(x) / (A^2 N^2 c)   [m],  A = e^{-u}, N = 1

    With eta = 0 the prefactor depends on position only through u, so
    f grad u = grad G(u) is exact and the integral must vanish -- the
    exactness control.  eta != 0 supplies the drift-landscape inhomogeneity
    that breaks exactness (d(phi_dot) misaligned with d phi).
    """
    total = 0.0
    nv = len(vertices)
    for i in range(nv):
        xa = vertices[i]
        xb = vertices[(i + 1) % nv]
        seg = xb - xa
        L_leg = np.linalg.norm(seg)
        ts = np.linspace(0.0, 1.0, npts)
        pts = xa[None, :] + seg[None, :] * ts[:, None]
        pu = pts[:, 0] * grad_u_vec[0] + pts[:, 1] * grad_u_vec[1]
        pe = pts[:, 0] * eta_dir[0] + pts[:, 1] * eta_dir[1]
        u = u0 + pu
        b = B_of_u(u, B0)
        udot = udot0 * (1.0 + eta * pe / L_leg) + udot_kin
        f = b * R_H**2 * udot / (np.exp(-2.0 * u) * C_SI)  # A^2 = e^{-2u}, N=1
        # du along the leg = grad_u . dx (linear: constant per leg)
        du_dx = grad_u_vec @ (seg / L_leg)
        total += np.trapezoid(f * du_dx * L_leg, ts)
    return total


def canonical_H(L_leg, u0, grad_u, B0, udot0, eta, udot_kin=0.0):
    """H [s] for a triangle of characteristic leg L_leg.  eta is the
    fractional drift variation across the loop, applied along the direction
    maximally misaligned with grad u (sin theta_m = 1)."""
    # equilateral triangle in the screen plane, side L_leg
    v = [np.array([0.0, 0.0]),
         np.array([L_leg, 0.0]),
         np.array([0.5 * L_leg, 0.5 * np.sqrt(3.0) * L_leg])]
    gv = np.array([grad_u, 0.0])
    ed = np.array([0.0, 1.0])   # drift ramp transverse to grad u (sin=1)
    integ = loop_integral(v, u0, gv, B0, udot0, eta, ed, udot_kin)
    return integ / C_SI


def main():
    set_step_logger(logger)
    print_status("STEP 057: canonical holonomy amplitude on the J0437 screen loop",
                 "INFO")

    # ------------------------------------------------------------------
    # Measured inputs: leg scale and psi from the step_003 products
    # ------------------------------------------------------------------
    per_epoch = json.load(open(RESULTS_DIR
                               / "step_003_closure_final_per_epoch_j0437.json"))
    tau_legs = np.abs([t[e] for ep in per_epoch for t in ep["triplets"]
                       for e in ("tau_01", "tau_12", "tau_02")])
    fd_legs = np.abs([t[e] for ep in per_epoch for t in ep["triplets"]
                      for e in ("fD_01", "fD_12", "fD_02")])
    psi_all = np.array([t["phase_closure_rad"] for ep in per_epoch
                        for t in ep["triplets"]])
    psi_circ = float(np.angle(np.mean(np.exp(1j * psi_all))))

    tau_med = float(np.median(tau_legs))
    tau_p95 = float(np.percentile(tau_legs, 95))
    L_med = float(leg_dx_from_tau(tau_med))
    L_p95 = float(leg_dx_from_tau(tau_p95))
    # fD leg cross-check: fD = v_eff dtheta / lambda -> dtheta = fD lambda / v
    lam = C_SI / NU_EFF
    dx_fd = float(np.median(fd_legs) * 1e-3 * lam / V_EFF * D_S)

    # implied transport delay under the measurement-model postulate
    psi_w = 0.984     # inverse-variance weighted detection statistic (Sec.4.1)
    psi_uw = 1.120    # unweighted frame-invariant statistic
    H_implied_w = abs(psi_w) / OMEGA_EFF
    H_implied_uw = abs(psi_uw) / OMEGA_EFF
    H_implied_circ = abs(psi_circ) / OMEGA_EFF

    # delay-channel bound (hierarchical amplitude, issues register / step_049)
    H_delay_bound = 0.459e-9   # s,  |H| <= 0.46 ns (sigma-amplitude channel)

    res = {
        "step": "step_057_holonomy_amplitude_closure",
        "connection": "sigma_i dx^i = f(x) du; f = -b(u) R_H^2 udot/(A^2 N^2 c); "
                      "b(u) = B0 u^2/(1+u^2) exp(-u^4/2 sigma_B^4), sigma_B=1.5; "
                      "H = (1/c) oint sigma . dx",
        "measured_inputs": {
            "n_triplets": int(len(psi_all)),
            "tau_leg_median_us": tau_med,
            "tau_leg_p95_us": tau_p95,
            "L_leg_median_m": L_med,
            "L_leg_p95_m": L_p95,
            "L_leg_median_AU": L_med / AU_M,
            "L_leg_fD_crosscheck_m": dx_fd,
            "psi_circular_mean_rad": psi_circ,
            "psi_weighted_rad": psi_w,
            "psi_unweighted_rad": psi_uw,
        },
        "implied_H": {
            "omega_eff_rad_s": OMEGA_EFF,
            "H_implied_weighted_s": H_implied_w,
            "H_implied_unweighted_s": H_implied_uw,
            "H_implied_circular_s": H_implied_circ,
            "note": "narrowband-equivalent delay IF psi were a transport "
                    "phase; the identification is the postulate under test",
        },
        "field_model": {
            "u_gal_absolute": float(U_GAL),
            "grad_u_galactic_m-1": float(GRAD_U_GAL),
            "delta_u_screen_excursion_bound": float(DELTA_U_ISM),
            "udot_drift_H0_s-1": H0,
            "udot_kinematic_s-1": float(V_EFF * GRAD_U_GAL),
            "D_s_pc": float(D_S / PC_M),
            "D_eff_pc": float(D_EFF / PC_M),
        },
    }

    # ------------------------------------------------------------------
    # Controls: exactness (eta=0), B0=0, traversal reversal
    # ------------------------------------------------------------------
    L = L_med
    tri = [np.array([0.0, 0.0]), np.array([L, 0.0]),
           np.array([0.5 * L, 0.5 * np.sqrt(3.0) * L])]
    gv = np.array([GRAD_U_GAL, 0.0])
    ed = np.array([0.0, 1.0])

    I_eta0 = loop_integral(tri, U_GAL, gv, 1.0, H0, 0.0, ed)
    I_B0_0 = loop_integral(tri, U_GAL, gv, 0.0, H0, 1.0, ed)
    I_rev = loop_integral(tri[::-1], U_GAL, gv, 1.0, H0, 1.0, ed)
    I_fwd = loop_integral(tri, U_GAL, gv, 1.0, H0, 1.0, ed)
    res["controls"] = {
        "eta0_exactness_integral_m": float(I_eta0),
        "B0_zero_integral_m": float(I_B0_0),
        "reversal_ratio_Irev_over_Ifwd": float(I_rev / I_fwd) if I_fwd else None,
        "note": "eta=0 must return ~0 (prefactor then depends on position "
                "only through u, so f grad u is exact); reversal must flip "
                "the sign; B0=0 must return 0.",
    }

    # ------------------------------------------------------------------
    # Canonical amplitude on the screen geometry, admissible branch
    # ------------------------------------------------------------------
    scenarios = {}
    for name, u0 in [("absolute_depth", U_GAL),
                     ("ambient_excursion", DELTA_U_ISM)]:
        for B0 in (1.0, 3.2e-3):
            H_eta1 = canonical_H(L, u0, GRAD_U_GAL, B0, H0, 1.0)
            H_eta1_p95 = canonical_H(L_p95, u0, GRAD_U_GAL, B0, H0, 1.0)
            # kinematic channel alone: udot_kin spatially uniform is exact;
            # include its spatially-variable part via the same eta structure
            H_kin_eta1 = canonical_H(L, u0, GRAD_U_GAL, B0, 0.0, 1.0,
                                     udot_kin=V_EFF * GRAD_U_GAL)
            key = f"{name}__B0_{B0:g}"
            scenarios[key] = {
                "u0": float(u0),
                "B0": B0,
                "b_u": float(B_of_u(u0, B0)),
                "H_eta1_s": float(H_eta1),
                "H_eta1_p95leg_s": float(H_eta1_p95),
                "H_kinematic_eta1_s": float(H_kin_eta1),
                "ratio_implied_over_canonical":
                    float(H_implied_w / abs(H_eta1)) if H_eta1 else None,
                "B0_required_at_eta1":
                    float(B0 * H_implied_w / abs(H_eta1)) if H_eta1 else None,
                "eta_required_at_this_B0":
                    float(H_implied_w / abs(H_eta1)) if H_eta1 else None,
            }
    res["canonical_scenarios"] = scenarios

    # ------------------------------------------------------------------
    # Terrestrial-isotropy consistency of the required coupling
    # ------------------------------------------------------------------
    B0_req = scenarios["absolute_depth__B0_1"]["B0_required_at_eta1"]
    b_E = abs(B_of_u(U_EARTH, B0_req))
    BM2 = b_E * R_H**2
    s0_req = math.sqrt(D_ISOTROPY_BOUND / BM2) / SIG_BARE_EARTH
    b_E_adm = abs(B_of_u(U_EARTH, 1.0))
    s0_adm = math.sqrt(D_ISOTROPY_BOUND / (b_E_adm * R_H**2)) / SIG_BARE_EARTH
    res["isotropy_consistency"] = {
        "u_earth": float(U_EARTH),
        "sigma_bare_surface_m-1": float(SIG_BARE_EARTH),
        "B0_required_for_implied_H": float(B0_req),
        "s0_pinning_required_at_B0_req": float(s0_req),
        "s0_pinning_required_at_B0_1": float(s0_adm),
        "s0_pinning_required_at_B0_cal": float(
            math.sqrt(D_ISOTROPY_BOUND
                      / (abs(B_of_u(U_EARTH, 3.2e-3)) * R_H**2))
            / SIG_BARE_EARTH),
    }

    # ------------------------------------------------------------------
    # J1603 counterpart (identical conventions, absolute-depth)
    # ------------------------------------------------------------------
    per_epoch2_path = RESULTS_DIR / "step_003_closure_final_per_epoch_j1603.json"
    j1603 = {}
    if per_epoch2_path.exists():
        ep2 = json.load(open(per_epoch2_path))
        tau2 = np.abs([t[e] for ep in ep2 for t in ep["triplets"]
                       for e in ("tau_01", "tau_12", "tau_02")])
        psi2 = np.array([t["phase_closure_rad"] for ep in ep2
                         for t in ep["triplets"]])
        tau_med_j2 = float(np.median(tau2))
        L2 = float(leg_dx_from_tau(tau_med_j2, d_eff=D_EFF2, d_s=D_S2))
        H2 = canonical_H(L2, U_GAL, GRAD_U_GAL, 1.0, H0, 1.0)
        psi2_circ = float(np.angle(np.mean(np.exp(1j * psi2))))
        j1603 = {
            "tau_leg_median_us": tau_med_j2,
            "L_leg_median_m": L2,
            "psi_circular_mean_rad": psi2_circ,
            "H_implied_s": abs(psi2_circ) / OMEGA_EFF,
            "H_canonical_B0_1_eta1_s": float(H2),
            "ratio": float(abs(psi2_circ) / OMEGA_EFF / abs(H2)) if H2 else None,
        }
    res["j1603"] = j1603

    # ------------------------------------------------------------------
    # Delay-channel bound comparison
    # ------------------------------------------------------------------
    H_c1 = scenarios["absolute_depth__B0_1"]["H_eta1_s"]
    res["delay_channel"] = {
        "H_signed_bound_s": H_delay_bound,
        "canonical_ceiling_B0_1_eta1_s": float(H_c1),
        "B0_bound_from_delay_at_eta1": float(H_delay_bound / abs(H_c1)),
        "statement": (
            "the signed-delay amplitude |H| <= 0.46 ns bounds the literal "
            "circulation on this sightline; the canonical admissible-branch "
            "ceiling (~5e-13 s at B0=1, order-unity drift inhomogeneity) "
            "sits ~3 orders below the bound, so the delay channel does not "
            "exclude the admissible branch -- it constrains it"),
    }

    res["verdict"] = {
        "implied_H_s": float(H_implied_w),
        "canonical_ceiling_absolute_depth_B0_1_s": float(H_c1),
        "gap_orders_of_magnitude":
            float(math.log10(H_implied_w / abs(H_c1))) if H_c1 else None,
        "landscape_smoothness_note": (
            "the corpus's own spatial landscape contrast, bounded by the "
            "observed redshift scatter at delta u ~ 1e-3-1e-2 per drift "
            "scale (Paper 0 flatness ledger), maps to a fractional drift "
            "inhomogeneity across the ~1e10 m loop of eta ~ "
            "delta u x L_leg/R_H ~ 1e-19, under which the canonical "
            "amplitude falls to ~1e-31 s; the eta=1 ceiling quoted is "
            "therefore a maximally generous bound, not the expected value"),
        "worldline_dip_note": (
            "a comparison link transported through the terrestrial region "
            "rather than along the screen does not rescue the amplitude: "
            "the same corpus consistency requirement that bounds the "
            "resonator channel pins the terrestrial shear response at "
            "S_Sigma <= 1e-10, suppressing any near-Earth contribution "
            "below the quoted ceiling"),
        "statement": (
            "Under the canonical disformal connection evaluated on the "
            "measured screen-triangle geometry (legs ~ %.2g m, u ~ 6e-7 "
            "galactic depth, grad u ~ galactic, udot ~ H0), the "
            "admissible-branch (B0 = +1) holonomy ceiling at order-unity "
            "drift inhomogeneity is ~%.1e s per loop -- ~%.2f orders below "
            "the narrowband-equivalent delay H_implied ~ %.1e s implied by "
            "psi = omega_eff H at omega_eff = 2pi x 1380 MHz.  Reaching the "
            "implied amplitude would require either eta ~ %.0f (drift "
            "inhomogeneity ~77 times the H0 scale across a 0.1 AU loop -- "
            "inconsistent with the landscape smoothness bounded by the "
            "redshift scatter) or B0 ~ %.0f, departing ~%.1f orders from "
            "the Paper-28 admissible normalisation and ~%.1f orders from "
            "the Paper-0 calibrated magnitude.  Under the "
            "ambient-excursion field convention selected by the Earth-loop "
            "consistency argument the ceiling is far lower still.  The "
            "coherent closure phase therefore cannot be read as the "
            "literal omega_eff x H_resid circulation on this geometry -- "
            "the amplitude reconciliation joins the parity result "
            "(orientation-even, no Stokes-flux scaling) in identifying psi "
            "as a phase-domain non-factorizability signature of the "
            "ordered leg product.  The direct time-domain holonomy "
            "constraint on this sightline is the delay-channel bound "
            "|H| <= 0.46 ns."
            % (L, abs(H_c1),
               math.log10(H_implied_w / abs(H_c1)) if H_c1 else float("nan"),
               H_implied_w,
               scenarios["absolute_depth__B0_1"]["eta_required_at_this_B0"],
               B0_req,
               math.log10(B0_req),
               math.log10(B0_req / 3.2e-3))),
    }

    outpath = RESULTS_DIR / "step_057_holonomy_amplitude_closure.json"
    with open(outpath, "w") as f:
        json.dump(res, f, indent=2, cls=NpEncoder)
    print_status(f"wrote {outpath}", "INFO")
    print(json.dumps(res["verdict"], indent=2, cls=NpEncoder))


if __name__ == "__main__":
    main()
