# Temporal Equivalence Principle: Synchronization Holonomy in Pulsar Scintillation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19454620.svg)](https://doi.org/10.5281/zenodo.19454620)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

![TEP-J0437: Pulsar Scintillation](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.2 (Sintra)  
**First published:** 17 May 2026 · **Last updated:** 17 May 2026  
**Status:** Preprint (Open for Collaboration)  
**DOI:** [10.5281/zenodo.19454620](https://doi.org/10.5281/zenodo.19454620)  
**Website:** [https://mlsmawfield.com/tep/j0437/](https://mlsmawfield.com/tep/j0437/)  
**Paper Series:** TEP Series: Paper 16 (Pulsar Scintillation)

## Abstract

Standard scintillation theory treats each scattered ray as carrying a scalar delay, so differential delays around a closed triplet cancel identically: \(\tau_{ij}+\tau_{jk}+\tau_{ki}=0\). The Temporal Equivalence Principle (TEP) instead predicts that proper-time transport is path-dependent in low-density, unscreened environments, producing a non-zero synchronization holonomy. This paper reports the rejection of the scalar-delay null hypothesis in pulsar scintillation using the phase-domain closure statistic \(\psi\), a zero-centered circular observable that separates geometric phase from folded-noise bias.

The primary target is PSR J0437-4715, analyzed with 19,167 scintillation triplets from 1,391 closure-capable Parkes/PPTA epochs; 1,093 epochs form the independent sample. The broader 15-pulsar catalog includes PSR J1603-7202, ten Jiamusi pulsars, and three MeerKAT pulsars. While J0437-4715 provides the primary phase-domain rejection of the \(\psi=0\) null, PSR J1603-7202 contributes complementary bipolar geometric structure from a distinct line of sight at different velocity geometry. Furthermore, the distant Jiamusi and MeerKAT samples are noise-limited, consistent with TEP's predicted environmental suppression in dense environments.

J0437-4715 shows a non-zero Phase Closure signal. The weighted circular mean is \(\bar{\psi}=0.984\pm0.046\) rad with \(R_{\rm bar}=0.308\); the directional V-test shows significant directional concentration relative to the pre-specified null direction \(\mu_0 = 0\) at \(p=2.04\times10^{-5}\), and the 95% bootstrap confidence interval \([0.737,1.235]\) rad excludes zero. The distribution is non-uniform (Rayleigh \(p=1.34\times10^{-44}\)), identical in heliocentric and CMB-frame analyses, confirming that the rejection is not a frame-dependent artifact.

PSR J1603-7202 has a 73.8° proper-motion separation from J0437-4715 and exhibits a frame-independent bipolar geometric structure, matching TEP predictions for high-dispersion sightlines where the monopole is washed out. The Jiamusi and MeerKAT samples are consistent with the expected environmental suppression at large distances. The raw unsigned delay magnitude for J0437-4715, \(|H|=8.100\pm0.102\) ns, operates at the expected folded-normal noise floor \(E[|H|]=6.810\) ns (see Section 2.1.1). A robust trimmed amplitude \(H_{\rm trim}=21.991\pm0.483\) ns (45.5σ) serves as a secondary robustness diagnostic; it does not carry primary inferential weight.

Multiple independent checks support the phase-domain result. Phase-scramble and pre-alignment controls pass, unweighted \(\psi\) is strictly frame-invariant, signed-delay cancellation behaves as expected for a bipolar signal, and rigorous synthetic noise tests confirm zero false positives. A signed-delay orbital diagnostic shows phase-locked structure directionally consistent with TEP kinematic coupling, though the hierarchical mixed-effects amplitude is not independently significant (\(p = 0.372\), 2 df), as expected for a partially screened orbital channel. Multi-pulsar scaling, chromaticity, cross-telescope environmental bounds, and orbital structure provide multiscale consistency checks on the TEP framework.

The principal empirical result is the detection of non-zero Phase Closure on the J0437 sightline. Standard scalar-delay ISM models—including thin-screen Kolmogorov, multi-screen, refractive wandering, chromatic plasma, and Doppler-delay covariance models—predict \(\psi = 0\) and are rejected by this observation. The geometric structure from J1603-7202 and the environmental-suppression consistency in the distant pulsars are directionally consistent with the Temporal Equivalence Principle's non-integrable time transport.

Code Availability: All data and analysis code required to reproduce the results presented in this work, including the full 15-pulsar scintillation catalog compilation, are available in the public repository at https://github.com/matthewsmawfield/TEP-J0437.

## Key Findings

Analysis of 19,167 scintillation triplets from 1,391 epochs (1,093 independent) for PSR J0437-4715 across a 15-pulsar catalog reveals a robust non-zero Phase Closure ψ signal consistent with the Temporal Equivalence Principle. For J0437: unweighted ψ = +1.120 rad (Rayleigh p = 1.34 × 10⁻⁴⁴, frame-invariant); circular mean ψ̄ = 0.984 ± 0.046 rad, R_bar = 0.308, Rayleigh p = 1.39 × 10⁻¹³, V-test p = 2.04 × 10⁻⁵; $H_{\rm trim} = 21.991 \pm 0.483$ ns (45.5σ). J1603-7202 provides independent geometric evidence through bipolar decomposition (83.2° orientation separation, bipole-to-monopole ratio 1.078). Jiamusi and MeerKAT pulsars provide noise-limited environmental bounds consistent with predicted suppression.

## Evidence Ledger

The pipeline now writes an explicit evidence ledger at `results/step_049_evidence_ledger.json`. This file separates claims by evidential status so the result can be reviewed without relying on rhetorical weight.

| Status | Claim | Current Result |
| --- | --- | --- |
| Primary evidence | J0437 rejects the additive scalar path-delay null through non-zero Phase Closure ψ | ψ̄ = 0.984 ± 0.046 rad; unweighted ψ = +1.120 rad; unweighted Rayleigh p = 1.34 × 10⁻⁴⁴; bootstrap CI [0.737, 1.235] rad excludes zero |
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
| **Paper 15** | [TEP-EFA](https://github.com/matthewsmawfield/TEP-EFA) | Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly | [10.5281/zenodo.19454863](https://doi.org/10.5281/zenodo.19454863) |
| **Paper 16** | **TEP-J0437** (This repo) | Synchronization Holonomy in Pulsar Scintillation | [10.5281/zenodo.19454620](https://doi.org/10.5281/zenodo.19454620) |
| **Paper 17** | [TEP-LLR](https://github.com/matthewsmawfield/TEP-LLR) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446029](https://doi.org/10.5281/zenodo.19446029) |

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
  note={Preprint v0.2 (Sintra)}
}
```

---

## Open Science Statement

This is a working preprint shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons and MIT licenses to encourage and facilitate replication. Feedback and collaboration are warmly invited and welcome.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
