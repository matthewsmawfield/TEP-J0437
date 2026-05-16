"""
TEP-J0437 Configuration Module
================================

Centralized configuration management for the TEP pipeline.
All tunable parameters are defined here for easy modification.
"""

import json
import functools
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any

from .json_numpy import NpEncoder

# Project root directory
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# Top-level constants for easy import
RANDOM_SEED = 42
RESULTS_DIR = PROJECT_ROOT / "results"

# Unit conversion factors (CODATA 2018 / IAU 2015) - must be defined before use
PC_TO_KM = 3.085677581491367e13  # 1 parsec in km (IAU 2015 definition)
MAS_YR_TO_KM_S = 4.74047e-3  # 1 mas/yr * pc to km/s (IAU 2015)
C_LIGHT_KM_S = 299792.458  # Speed of light in km/s (CODATA 2018 exact)

# Pulsar and Screen Parameters
J0437_PSI_DEG = 128.0  # Principal scattering axis (degrees, kinematic alignment)
# NOTE: NU_REF_MHZ removed - frequency must be loaded from epoch metadata
# Hardcoding 1400 MHz would compromise chromatic vs achromatic discrimination tests

# Physical Constants for J0437-4715 (from Reardon et al. 2021, ApJ, 923, 32)
# Used in step_008 alternative explanations testing
C_LIGHT = 2.99792458e8  # Speed of light in vacuum [m/s] (CODATA 2018 exact)
LAMBDA_LBAND = 0.21  # Wavelength [m] for L-band observations (approx 1.4 GHz)
D_EFF = 37.5 * PC_TO_KM * 1e3  # Effective distance to scattering screen [m] (37.5 pc from Reardon et al. 2021)
V_EFF = 104.38e3  # Effective transverse velocity [m/s] (104.38 km/s from Reardon et al. 2021)
# Fresnel scale calculated from wavelength and effective distance
FRESNEL = np.sqrt(LAMBDA_LBAND * D_EFF / (2 * np.pi))

# Pulsar Physical Parameters (from peer-reviewed literature)
# J0437-4715 parameters (Reardon et al. 2021, ApJ, 923, 32; Deller et al. 2008, ApJ, 685, 674)
J0437_DIST_PC = 156.3  # Distance [pc] (parallax from Deller et al. 2008)
J0437_PM_RA = 121.439  # Proper motion RA [mas/yr] (Reardon et al. 2021)
J0437_PM_DEC = -71.438  # Proper motion Dec [mas/yr] (Reardon et al. 2021)
J0437_PM_MAG = 106.9  # Proper motion magnitude [mas/yr] (sqrt(pm_ra^2 + pm_dec^2))
J0437_PB_DAYS = 5.74104635  # Orbital period [days] (Reardon et al. 2024, PPTA-DR3)
# Epoch of periastron T0 [MJD] — must match the same timing solution as omega, i, Omega below
J0437_T0_MJD = 54530.1722  # Reardon et al. 2024, PPTA-DR3 (ApJL 971 L19; arXiv:2407.07132)
J0437_A1_LC = 3.36671466  # Projected semi-major axis x = (a sin i)/c [light-seconds], PPTA-DR3
J0437_S_SCREEN = 0.6  # Relative screen distance (D_s / D_p) from Reardon et al. 2021

# Binary orientation (same PPTA-DR3 solution as T0 / x above). Used for sky-plane reflex velocity.
J0437_ECC = 1.91805e-5
J0437_OM_DEG = 1.359  # Longitude of periastron ω (deg)
J0437_KOM_DEG = 208.3  # Longitude of ascending node Ω (deg)
J0437_INC_DEG = 137.506  # Orbital inclination i (deg)
# J2000 position for equatorial projection of v_orbit (Reardon et al. 2024, Table)
_J0437_RA_H = 4.0 + 37.0 / 60.0 + 15.9284042 / 3600.0
_J0437_DEC_D = -(47.0 + 15.0 / 60.0 + 9.303700 / 3600.0)
J0437_RA_RAD = float(np.deg2rad(_J0437_RA_H * 15.0))  # hours -> degrees -> rad
J0437_DEC_RAD = float(np.deg2rad(_J0437_DEC_D))

# J1603-7202 parameters (Walker et al. 2022, MNRAS, 510, 3462)
J1603_DIST_PC = 250.0  # Distance [pc] (DM-based from Walker et al. 2022)
J1603_PM_RA = -6.6  # Proper motion RA [mas/yr] (Walker et al. 2022)
J1603_PM_DEC = -25.9  # Proper motion Dec [mas/yr] (Walker et al. 2022)
J1603_PB_DAYS = 0.0  # Non-binary
J1603_T0_MJD = 0.0  # Non-binary
J1603_A1_LC = 0.0  # Non-binary
J1603_K_KMS = 0.0  # Non-binary
J1603_S_SCREEN = 0.5  # Relative screen distance (D_s / D_p) from Walker et al. 2022

# J0613-0200 parameters (Deller et al. 2009, ApJ, 701, 1243; Reardon et al. 2021)
# Binary MSP in a higher-density ISM environment than J0437; intended as a
# real-data control pulsar when PPTA DR2 epochs are available.
J0613_DIST_PC = 900.0  # Distance [pc] (parallax from Deller et al. 2009)
J0613_PM_RA = 3.0  # Proper motion RA [mas/yr] (Reardon et al. 2021)
J0613_PM_DEC = -9.0  # Proper motion Dec [mas/yr] (Reardon et al. 2021)
J0613_PB_DAYS = 1.198512595  # Orbital period [days] (Desvignes et al. 2016)
J0613_T0_MJD = 52617.112  # Epoch of periastron [MJD]
J0613_A1_LC = 1.525  # Projected semi-major axis [light-seconds]
J0613_S_SCREEN = 0.5  # Assumed relative screen distance
J0613_ENV_DENSITY = 1e-23  # ISM density [g/cm^3] (estimated ~10x J0437)
J0613_SCREEN_DIST_PC = 300.0  # Estimated screen distance [pc]

# Default configuration
DEFAULT_CONFIG = {
    "pipeline": {
        "random_seed": 42,
        "verbose": True,
        "log_dir": "logs",
        "results_dir": "results",
        "data_dir": "data"
    },
    "step_002_secondary_spectra": {
        "min_snr": 1.5,
        "max_snr": 50.0,
        "eta_threshold": 0.001,
        "tau_min_ns": 0.001,
        "tau_max_ns": 100.0,
        "fd_min_mhz": 0.001,
        "fd_max_mhz": 100.0,
        "use_direct_peak_finding": True,
        "min_arclets_for_triplets": 3,
        "max_arclets_for_triplets": 20
    },
    "step_003_closure_delays": {
        "max_closure_delay_us": 0.05,
        "min_snr": 5.0,
        "outlier_threshold_sigma": 5.0,
        "enable_bootstrap": True,
        "n_bootstrap": 10000,
        "confidence_level": 0.95
    },
    "analysis": {
        "significance_threshold": 5.0,
        "magnitude_equality_threshold": 0.1,
        "min_epochs_for_analysis": 10,
        "min_triplets_for_analysis": 100,
        "theory_status": "a_priori - TEP theory developed before this analysis",
        "methodology_status": "developed to test existing TEP theory",
        "analysis_tests": [
            {
                "name": "aggregate_mean_all_triplets",
                "type": "supporting",
                "purpose": "Verify bipolar cancellation"
            },
            {
                "name": "h_magnitude_detection",
                "type": "primary",
                "purpose": "Test a priori TEP prediction of non-zero holonomy"
            },
            {
                "name": "negative_delays_clockwise",
                "type": "supporting",
                "purpose": "Test clockwise loop component"
            },
            {
                "name": "positive_delays_counterclockwise",
                "type": "supporting",
                "purpose": "Test counter-clockwise loop component"
            }
        ],
        "multiple_comparison_correction": "bonferroni",
        "n_tests": 4,
        "alpha_uncorrected": 0.05,
        "falsification_criteria": {
            "h_magnitude_zero": "If |H| consistent with zero after correction",
            "single_sign_significant": "If only one sign shows significance",
            "magnitude_equality_fails": "|H-| significantly different from |H+|",
            "sign_distribution_skewed": "If sign distribution not ~50/50",
            "both_signs_significant": "If both signs equally significant (systematic error)"
        }
    },
    "visualization": {
        "figure_dpi": 300,
        "figure_format": "png",
        "color_scheme": "viridis",
        "show_progress": True
    }
}


class TEPConfig:
    """Configuration manager for TEP pipeline."""
    
    def __init__(self, config_path: Path = None):
        """
        Initialize configuration.
        
        Parameters
        ----------
        config_path : Path, optional
            Path to custom configuration JSON file
        """
        self.config = DEFAULT_CONFIG.copy()
        
        if config_path and config_path.exists():
            self.load_from_file(config_path)
    
    def load_from_file(self, path: Path):
        """Load configuration from JSON file."""
        with open(path, 'r') as f:
            user_config = json.load(f)
        # Merge with defaults
        self._merge_config(self.config, user_config)
    
    def _merge_config(self, default: Dict, user: Dict):
        """Recursively merge user config into default."""
        for key, value in user.items():
            if key in default and isinstance(default[key], dict) and isinstance(value, dict):
                self._merge_config(default[key], value)
            else:
                default[key] = value
    
    def get(self, section: str, key: str = None) -> Any:
        """Get configuration value."""
        if key is None:
            return self.config.get(section, {})
        return self.config.get(section, {}).get(key)
    
    def set(self, section: str, key: str, value: Any):
        """Set configuration value."""
        if section not in self.config:
            self.config[section] = {}
        self.config[section][key] = value
    
    def save(self, path: Path):
        """Save current configuration to file."""
        with open(path, 'w') as f:
            json.dump(self.config, f, indent=2, cls=NpEncoder)


# Global configuration instance
_config = None

def get_config(config_path: Path = None) -> TEPConfig:
    """Get or create global configuration instance."""
    global _config
    if _config is None:
        _config = TEPConfig(config_path)
    return _config
