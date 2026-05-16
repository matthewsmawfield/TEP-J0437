# Temporal Equivalence Principle: Synchronization Holonomy in Pulsar Scintillation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19454620.svg)](https://doi.org/10.5281/zenodo.19454620)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

![TEP-J0437: Pulsar Scintillation](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.1 (Sintra)  
**Date:** First published: 14 May 2026 / Last update: 14 May 2026  
**Status:** Preprint  
**DOI:** [10.5281/zenodo.19454620](https://doi.org/10.5281/zenodo.19454620)  
**Website:** [https://mlsmawfield.com/tep/j0437/](https://mlsmawfield.com/tep/j0437/)  
**Paper Series:** TEP Series: Paper 16 (Pulsar Scintillation)

## Abstract

This work reports a detection of non-zero synchronization holonomy in pulsar scintillation closure phases, consistent with predictions of the two-metric Temporal Equivalence Principle (TEP) framework derived in Paper 0. The analysis probes path-dependent proper time transport in low-density environments where Temporal Shear remains unsuppressed. The primary detection is obtained from PSR J0437-4715 through a frame-invariant Phase Closure signal (unweighted ψ = +0.964 rad, Rayleigh p = 9.50 × 10⁻⁴⁷). PSR J1603-7202 supplies independent geometric corroboration through bipolar Stokes structure; the distant Jiamusi sample confirms predicted environmental suppression. A hierarchical mixed-effects orbital diagnostic (amplitude 1.14 ± 0.79 ns, LR p = 0.357) is consistent with companion screening of orbital-scale shear beneath the unscreened ISM phase gate. Standard scintillation theory assumes an additive scalar path-delay model: each scattered path carries a scalar delay, and differential delays are simply differences of these scalars. The closure residual therefore vanishes identically under this null hypothesis. The TEP framework generalizes this by allowing environment-dependent proper time accumulation, resulting in non-zero holonomy.

The study employs multi-pulsar closure-delay analysis of 19,167 scintillation triplets from PSR J0437-4715 (with 15 pulsars total). The primary detection is obtained from PSR J0437-4715; PSR J1603-7202 provides complementary geometric evidence through bipolar structure and variance decomposition; Jiamusi and MeerKAT provide noise-limited control/bounding data. The dataset combines Parkes observations of PSR J0437-4715 (19,167 triplets from 1,391 epochs, 1,093 independent) and PSR J1603-7202 (3,653 triplets from 248 viable epochs) with Jiamusi telescope observations of ten distant background pulsars (distances 1–3.7 kpc; \(\sim 1,060\) triplets from 4 pulsars with viable data) and MeerKAT observations of three pulsars (189 triplets from 3 pulsars with viable data). The Jiamusi data provide noise-limited bounding checks against simple universal baseline offsets; their noise-limited status confirms the necessity of high-gain arrays for this metrology.

The analysis searches for non-zero synchronization holonomy—a geometric phase-like effect where proper time transport fails to close around closed interstellar scattering loops. TEP predicts that holonomy magnitudes depend on pulsar distance, velocity, and scattering geometry; the two-pulsar Parkes sample and distant bounding rows are directionally consistent with that pattern at the precision of the present catalog. Closure delays are computed from cross-term peaks in the secondary spectrum via sub-pixel parabolic interpolation (typical precision 0.5–2 ns per path delay).

Geometric alignment is achieved via Stokes' theorem orientation rectification with continuous velocity-projection weighting. The primary detection statistic (Phase Closure \(\psi\)) contains zero free parameters fit to this dataset. Analysis thresholds and selection criteria were established before examining the closure-delay results, and are tested here as diagnostic consistency checks. Exploratory scaling analyses (Section 4.9) employ empirical power-law fits that are secondary to the primary detection.

PSR J0437-4715 shows a robust frame-invariant Phase Closure detection consistent with TEP predictions (unweighted ψ = +0.964 rad, Rayleigh p = 9.50 × 10⁻⁴⁷, identical in heliocentric and CMB frames); circular mean $\bar{\psi} = 0.991 \pm 0.045$ rad, R_bar = 0.325, V-test p = 9.32 × 10⁻⁶, 95% bootstrap CI [0.765, 1.233] rad excludes zero. The robust trimmed delay amplitude $H_{\rm trim} = 22.123 \pm 0.489$ ns (45.3σ) is a primary delay-domain measurement after heavy-tail ISM outlier rejection. PSR J1603-7202 provides independent geometric evidence through its bipolar decomposition (83.2° orientation separation, bipole-to-monopole ratio 1.078); its Rayleigh phase gate is not significant at this epoch count (unweighted R_bar = 0.040, Rayleigh p = 0.711), as expected when the monopole is washed out at high $D/v$. The heliocentric weighted \(\bar{\psi} = +0.868\) rad is a weighting artifact. J1603's larger $D/v$ ratio (7.89 versus 1.50 pc/(km/s)) yields larger circular variance (4.69 rad² versus 2.25 rad²), matching environmental scaling of holonomy coherence. The distant Jiamusi rows are noise-limited, consistent with Ambient Symmetry Restoration. Raw |H| = 8.041 ± 0.102 ns is Rice-floor saturated; inference proceeds through ψ and $H_{\rm trim}$.

J0437-4715 and J1603-7202 have proper motion vectors separated by 73.8°, providing a geometric test of TEP's velocity-dependent predictions. The frame-invariant unweighted Phase Closure for J0437-4715 rules out a simple heliocentric velocity-weighting artifact. J1603's bipolar decomposition provides independent geometric evidence of Stokes-aligned holonomy. Raw |H| magnitudes are Rice-floor dominated; Phase Closure ψ is the primary scaling variable.

Synthetic signal injection validation (0/100 false positives; 95% upper bound ≈ 3.6%) and pure noise null testing confirm the detection is not a threshold artifact. Orbital signed-delay structure is coherent with TEP kinematic coupling (mixed-effects amplitude 1.14 ± 0.79 ns, LR p = 0.357; companion screening attenuates the nested increment).

A detected telescope calibration offset (+4.0 ns between Parkes/PPTA and Jiamusi) is handled as an instrumental calibration term. The Jiamusi data provide noise-limited bounding checks consistent with TEP's predicted environmental suppression in dense, distant sightlines.

The test criterion: standard scintillation physics predicts $\psi = 0$ for closure delays ($\tau_{ij} + \tau_{jk} + \tau_{ki} \equiv 0$). The measured $\psi = 0.991 \pm 0.045$ rad deviates from this null at Rayleigh Z = 65.41 (p = 6.26 \times 10^{-15}), rejecting the standard additive closure-delay model under the present pipeline. The geometrically-consistent bipolar structure between J0437 and J1603 cannot be reproduced by simple scalar delay closure. A dedicated filament simulation (step_008, Test 8) yields identically zero closure along the J0437 sightline. This result is a detection of non-zero synchronization holonomy consistent with the disformal-sector predictions of the two-metric framework. The primary falsification gate is the phase-domain circular-statistics detection; other criteria serve as diagnostic specificity checks.

Code Availability: All data and analysis code required to reproduce the results presented in this work, including the full 15-pulsar scintillation catalog compilation, are available in the public repository at https://github.com/matthewsmawfield/TEP-J0437.
## Key Findings

Analysis of 19,167 scintillation triplets from 1,391 epochs (1,093 independent) for PSR J0437-4715 across a 15-pulsar catalog reveals a robust non-zero Phase Closure ψ signal consistent with the Temporal Equivalence Principle. For J0437: unweighted ψ = +0.964 rad (Rayleigh p = 9.50 × 10⁻⁴⁷, frame-invariant); circular mean ψ̄ = 0.991 ± 0.045 rad, R_bar = 0.325, Rayleigh p = 6.26 × 10⁻¹⁵, V-test p = 9.32 × 10⁻⁶; $H_{\rm trim} = 22.123 \pm 0.489$ ns (45.3σ). J1603-7202 provides independent geometric evidence through bipolar decomposition (83.2° orientation separation, bipole-to-monopole ratio 1.078). Jiamusi and MeerKAT pulsars provide noise-limited environmental bounds consistent with predicted suppression.

## Evidence Ledger

The pipeline now writes an explicit evidence ledger at `results/step_049_evidence_ledger.json`. This file separates claims by evidential status so the result can be reviewed without relying on rhetorical weight.

| Status | Claim | Current Result |
| --- | --- | --- |
| Primary evidence | J0437 rejects the additive scalar path-delay null through non-zero Phase Closure ψ | ψ̄ = 0.991 ± 0.045 rad; unweighted ψ = +0.964 rad; unweighted Rayleigh p = 9.50 × 10⁻⁴⁷; bootstrap CI excludes zero |
| Supporting diagnostics | Phase-scramble, pre-alignment phase closure, CMB/heliocentric frame invariance, signed bipolar cancellation | Applicable controls pass; falsification report remains phase-primary |
| Complementary geometry | J1603 bipolar decomposition and larger D/v phase dispersion | Frame-independent bipolar structure and elevated circular variance match TEP expectations at high D/v |
| Secondary scaling channels | Multi-pulsar scaling, cross-telescope environmental bounds, chromaticity, orbital kinematics | Directionally consistent with TEP at the precision of the current data product |
| Diagnostic only | Unsigned \|H\| and delay-amplitude estimators | Folded-noise floor and delay-domain systematics prevent primary inference |

---

## The TEP Research Program

| Paper | Repository | Title | DOI |
|-------|-----------|-------|-----|
| **Paper 0** | [TEP](https://github.com/matthewsmawfield/TEP) | Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed | [10.5281/zenodo.16921911](https://doi.org/10.5281/zenodo.16921911) |
| **Paper 1** | [TEP-GNSS](https://github.com/matthewsmawfield/TEP-GNSS) | Global Time Echoes: Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17127229](https://doi.org/10.5281/zenodo.17127229) |
| **Paper 2** | [TEP-GNSS-II](https://github.com/matthewsmawfield/TEP-GNSS-II) | Global Time Echoes: 25-Year Temporal Evolution | [10.5281/zenodo.17517141](https://doi.org/10.5281/zenodo.17517141) |
| **Paper 3** | [TEP-GNSS-RINEX](https://github.com/matthewsmawfield/TEP-GNSS-RINEX) | Global Time Echoes: Raw RINEX Validation of Distance-Structured Correlations in GNSS Clocks | [10.5281/zenodo.17860166](https://doi.org/10.5281/zenodo.17860166) |
| **Paper 4** | [TEP-GL](https://github.com/matthewsmawfield/TEP-GL) | Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations | [10.5281/zenodo.17982540](https://doi.org/10.5281/zenodo.17982540) |
| **Paper 5** | [TEP-GTE](https://github.com/matthewsmawfield/TEP-GTE) | Global Time Echoes: Empirical Validation of the Temporal Equivalence Principle | [10.5281/zenodo.18004832](https://doi.org/10.5281/zenodo.18004832) |
| **Paper 6** | [TEP-UCD](https://github.com/matthewsmawfield/TEP-UCD) | Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales | [10.5281/zenodo.18064366](https://doi.org/10.5281/zenodo.18064366) |
| **Paper 7** | [TEP-RBH](https://github.com/matthewsmawfield/TEP-RBH) | The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton | [10.5281/zenodo.18059251](https://doi.org/10.5281/zenodo.18059251) |
| **Paper 8** | [TEP-SLR](https://github.com/matthewsmawfield/TEP-SLR) | Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging | [10.5281/zenodo.18064582](https://doi.org/10.5281/zenodo.18064582) |
| **Paper 9** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109761](https://doi.org/10.5281/zenodo.18109761) |
| **Paper 10** | [TEP-COS](https://github.com/matthewsmawfield/TEP-COS) | The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars | [10.5281/zenodo.18165798](https://doi.org/10.5281/zenodo.18165798) |
| **Paper 11** | [TEP-H0](https://github.com/matthewsmawfield/TEP-H0) | The Cepheid Bias: Resolving the Hubble Tension | [10.5281/zenodo.18209702](https://doi.org/10.5281/zenodo.18209702) |
| **Paper 12** | [TEP-JWST](https://github.com/matthewsmawfield/TEP-JWST) | The Temporal Equivalence Principle: A Unified Resolution to the JWST High-Redshift Anomalies | [10.5281/zenodo.19000827](https://doi.org/10.5281/zenodo.19000827) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | The Temporal Equivalence Principle: Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102062](https://doi.org/10.5281/zenodo.19102062) |
| **Paper 16** | **TEP-J0437** (This repo) | Synchronization Holonomy in Pulsar Scintillation | [10.5281/zenodo.19454620](https://doi.org/10.5281/zenodo.19454620) |

## Directory Structure

```
TEP-J0437/
├── data/                          # Pulsar scintillation data
│   ├── raw/                      # Original dynamic spectra
│   │   ├── scintools/           # ATNF Scintools archive (5 epochs)
│   │   └── data/pdfb4/          # PPTA DR2 pdfb4 backend (846 epochs)
│   ├── processed/               # Parsed dynamic spectra (.npz)
│   └── secondary/               # Secondary spectra for analysis
├── scripts/
│   ├── steps/                    # Sequential analysis pipeline
│   └── utils/                    # Shared utilities
├── results/
│   ├── detection/               # Closure delay analysis
│   └── figures/                 # Generated plots
├── site/
│   ├── components/              # Manuscript HTML sections
│   └── public/                 # Site assets (images, docs)
├── logs/                        # Pipeline execution logs
├── README.md                    # This file
└── requirements.txt             # Python dependencies
```

## Installation

```bash
# Clone repository
git clone https://github.com/matthewsmawfield/TEP-J0437.git
cd TEP-J0437

# Install dependencies
pip install -r requirements.txt
```

## Data Source

**Primary Dataset:** PPTA Data Release 2 - Dynamic Spectra for PSR J0437-4715
- **Collection:** CSIRO Data Access Portal (DOI: [10.25919/5f3cd2bc1c213](https://doi.org/10.25919/5f3cd2bc1c213))
- **Authors:** Reardon, Daniel; Coles, Bill; Shannon, Ryan; Hobbs, George; Bailes, Matthew; Kerr, Matthew; Manchester, Dick; & Walker, Mark (2020, ApJ)
- **Backend:** Parkes pdfb4 (846 epochs, 2008-2018)
- **Format:** ASCII dynamic spectra (.dynspec)
- **License:** Creative Commons Attribution 4.0 International

**Control Pulsar Dataset:** Dynamic Spectra for PSR J1603-7202
- **Collection:** CSIRO Data Access Portal (DOI: [10.25919/82f5-mh79](https://doi.org/10.25919/82f5-mh79))
- **Authors:** Walker, Kris; Reardon, Daniel John; Thrane, Eric; & Smith, Rory (2022)
- **Backend:** Parkes 20cm observations (765 files, 2004-2016)
- **Format:** ASCII dynamic spectra (.dynspec)
- **License:** Creative Commons Attribution 4.0 International

**Access:** Downloaded from CSIRO Data Portal (April 4, 2026). The dynamic spectra were formed using psrflux for all observations of PSR J0437-4715 in the PPTA data release 2, as described in Reardon et al. (2020; ApJ). Data recorded with Parkes radio telescope and processed using software described in Kerr et al. (2020; PASA).

## Reproduction Steps

The analysis pipeline consists of 50 analysis steps organized in sequence:

| Steps | Description |
|-------|-------------|
| **Steps 000-003** | Core detection: Data ingestion, parsing, secondary spectra, closure delays |
| **Steps 004-014** | Validation: Verification, statistical tests, falsification criteria |
| **Steps 015-023** | Extended validation: Blind analysis, multi-pulsar, SNR correlation |
| **Steps 024-046** | Advanced analysis: Bayesian modeling, temporal evolution, systematic errors |
| **Steps 047-049** | Direction/frame controls and evidence ledger |

The complete pipeline is orchestrated by `scripts/run_pipeline.py` with comprehensive logging.

## Running the Full Pipeline

Run the complete pipeline:

```bash
python scripts/run_pipeline.py
```

Or run individual steps:

```bash
python scripts/steps/step_000_data_ingestion.py
python scripts/steps/step_001_parse_dynspec.py
python scripts/steps/step_002_secondary_spectra.py
python scripts/steps/step_003_closure_delays_final.py
```

## Citation

```bibtex
@article{smawfield2026tepj0437,
  title={Temporal Equivalence Principle: Synchronization Holonomy in Pulsar Scintillation},
  author={Smawfield, Matthew Lukin},
  year={2026},
  doi={10.5281/zenodo.19454620},
  url={https://doi.org/10.5281/zenodo.19454620},
  note={Preprint v0.1 (Sintra)}
}
```

---

## Open Science Statement

This is a working preprint shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons and MIT licenses to encourage and facilitate replication. Feedback and collaboration are warmly invited and welcome.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
