"""
Physics constants and fundamental units for the Temporal Equivalence Principle (TEP) framework.
Standardized for Jakarta v0.8.

References:
- Jakarta v0.8 Manuscript (2026)
- NIST CODATA 2018
- WGS84 Geodetic System
"""

from typing import Final

# Fundamental Constants
C_LIGHT: Final[float] = 299792458.0  # Speed of light in vacuum [m/s]
G_NEWTON: Final[float] = 6.67430e-11  # Newton's gravitational constant [m³ kg⁻¹ s⁻²]
H_BAR: Final[float] = 1.054571817e-34  # Reduced Planck constant [J s]

# Planck Units (Natural Units for TEP Scalar Field)
# Standardized to Reduced Planck Mass for Jakarta v0.8 consistency
M_PL_GEV: Final[float] = 2.435e18  # Reduced Planck Mass [GeV]
M_PL_KG: Final[float] = 4.341e-9    # Reduced Planck Mass [kg]

# Earth Parameters (WGS84)
R_EARTH: Final[float] = 6371000.0   # Mean Earth radius [m]
M_EARTH: Final[float] = 5.9722e24   # Earth mass [kg]
GM_EARTH: Final[float] = 3.986004418e14  # Earth gravitational parameter [m³ s⁻²]
J2_EARTH: Final[float] = 1.08263e-3  # Earth's second dynamic form factor (J2)

# Astronomy Units
AU_METERS: Final[float] = 1.495978707e11  # Astronomical Unit [m]
PC_METERS: Final[float] = 3.08567758e16   # Parsec [m]
LY_METERS: Final[float] = 9.46073047e15   # Light year [m]

# TEP Baseline Theory Parameters (Jakarta v0.8)
LAMBDA_BASELINE_GEV: Final[float] = 1.0e-5  # 10 keV scale [GeV]
BETA_BASELINE: Final[float] = 1.0           # Unit scalar coupling
N_TOPOLOGY: Final[int] = 1                   # Continuous gradient suppression index

# TEP Field Relaxation and Screening (Empirical Anchors)
LAMBDA_TEP_M: Final[float] = 4200000.0       # Unified relaxation length [m]
R_TRANSITION_M: Final[float] = 4146000.0     # PREM-derived transition radius [m]
RHO_CRITICAL_G_CM3: Final[float] = 20.0      # Universal critical density [g/cm³]
SUPPRESSION_EXPONENT: Final[float] = 0.334    # Theoretical density scaling (≈ 1/3)

# Screening Factor (Jakarta v0.8 standard)
# S_⊕ = (R_transition / R_Earth)^(1/3) ≈ (4146/6371)^(1/3) ≈ 0.349
CHARACTERISTIC_SUPPRESSION: Final[float] = 0.349

# Conversion Factors
KG_M3_TO_GEV4: Final[float] = 4.318e-21  # kg/m³ to GeV⁴ conversion (Natural Units)

def get_tep_metadata() -> dict:
    """Return theory version and metadata for data provenance."""
    return {
        "theory_version": "Jakarta v0.8",
        "paradigm": "Temporal Topology (Continuous Gradient)",
        "coupling_convention": "A(phi) = exp(2 beta phi / M_Pl)",
        "mpl_definition": "Reduced Planck Mass (2.435e18 GeV)",
        "standard_constants": {
            "Lambda_TEP": f"{LAMBDA_TEP_M/1e3} km",
            "R_transition": f"{R_TRANSITION_M/1e3} km",
            "rho_c": f"{RHO_CRITICAL_G_CM3} g/cm³",
            "S_earth": f"{CHARACTERISTIC_SUPPRESSION}"
        }
    }
