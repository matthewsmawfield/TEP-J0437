# Temporal Equivalence Principle: Synchronization Holonomy in Pulsar Scintillation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19454619.svg)](https://doi.org/10.5281/zenodo.19454619)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
![TEP-J0437: Pulsar Scintillation](site/public/image.webp)

**Author:** Matthew Lukin Smawfield  
**Version:** v0.6 (Sintra)  
**First published:** 17 May 2026 · **Last updated:** 24 October 2026  
**Status:** Preprint  
**DOI:** [10.5281/zenodo.19454619](https://doi.org/10.5281/zenodo.19454619)  
**Website:** [https://mlsmawfield.com/tep/j0437/](https://mlsmawfield.com/tep/j0437/)  
**Paper Series:** TEP Series · Group E (V. Pulsars, Astrophysical Clocks & Distance-Ladder Tests) · Paper 16 (`E16`)

## Abstract

<!-- Generated from site/components/0_abstract.html; update the manuscript source, then run tools/sync_current_corpus_readmes.py. -->

The scalar-factorizable thin-screen delay model tested here assigns each scattered ray an additive scalar delay, so differential delays around a closed triplet cancel identically: \(\tau_{ij}+\tau_{jk}+\tau_{ki}=0\). The Temporal Equivalence Principle (TEP) instead predicts that proper-time transport is path-dependent in dilute, unscreened environments, producing a non-zero synchronization holonomy. This paper reports the rejection of the scalar-delay null hypothesis in pulsar scintillation using the phase-domain closure statistic \(\psi\), a zero-centered circular observable that separates geometric phase from folded-noise bias.

The primary target is PSR J0437-4715, analyzed with 19,167 scintillation triplets from 1,391 closure-capable Parkes/PPTA epochs; 1,093 epochs form the independent sample. The broader 15-pulsar catalog includes PSR J1603-7202, ten Jiamusi pulsars, and three MeerKAT pulsars. While J0437-4715 provides the primary phase-domain rejection of the \(\psi=0\) null, PSR J1603-7202 contributes complementary bipolar geometric structure from a distinct line of sight at different velocity geometry. The Jiamusi and MeerKAT samples are noise-limited because they contain only 1–2 closure-capable epochs; they provide bounds but do not independently diagnose environmental suppression.

J0437-4715 shows a non-zero Phase Closure signal. Because the Rayleigh and V-test statistics measure circular concentration against uniformity (and are rotation-invariant with respect to the mean angle), circular concentration and zero-mean rejection are reported separately: the unweighted circular mean is \(\psi=1.120\) rad (95% bootstrap CI \([+0.988,+1.254]\) rad excluding zero at \(5.7\sigma\)–\(12.7\sigma\) across 15–240-epoch contiguous MJD blocks and \(16.8\sigma\) i.i.d.; circular concentration Rayleigh \(p=1.34\times10^{-44}\)), and the weighted circular mean is \(\bar{\psi}=0.984\pm0.046\) rad with 95% circular-bootstrap confidence interval \([0.737,1.235]\) rad excluding zero at \(5.2\sigma\)–\(6.6\sigma\) across 15–240-epoch MJD blocks (\(7.9\sigma\) i.i.d. at effective sample size \(n_{\rm eff}=313\), \(\text{SE}_{\rm boot}=0.125\) rad; \(21.2\sigma\) analytical circular SE, Step 059; corroborated in the delay domain by the MJD-ordered 15–240-epoch block-bootstrap noise-subtracted excess at \(3.8\sigma\)–\(10.8\sigma\), Steps 052 and 060; circular concentration \(\bar{R}=0.308\), Rayleigh \(p=1.39\times10^{-13}\), V-test \(p=1.02\times10^{-5}\) relative to \(\mu_0 = 0\)). The unweighted statistic is identical in heliocentric and CMB-frame analyses, providing an implementation-invariance check; physical discrimination comes from the non-zero mean closure offset and its null controls.

PSR J1603-7202 has a 73.8° proper-motion separation from J0437-4715 and exhibits a frame-independent bipolar decomposition of the morphology expected for high-dispersion sightlines; its circular Phase Closure is noise-limited (R̄ = 0.019, Rayleigh p = 0.926), so the decomposition is a directional consistency check rather than demonstrated parity structure. Delay-domain unsigned amplitude estimators (\(M_{\rm cl}\)) and their folded-normal Rice noise-floor bounds are documented in Appendix B.

Multiple independent checks support the phase-domain result. Phase-scramble and pre-alignment controls pass, unweighted \(\psi\) is strictly frame-invariant, signed-delay cancellation behaves as expected for a bipolar signal, and rigorous synthetic noise tests yielded zero false positives in 100 synthetic null trials. In the 5.74-day binary orbital channel, while the naive signed-delay mean is attenuated by within-epoch bipolar cancellation (\(1.11 \pm 0.79\) ns, \(p = 0.372\), 2 df), evaluating the identical REML mixed-effects model on the unsigned delay magnitude \(M_{\rm cl}\) yields \(A = 5.71 \pm 1.10\) ns (\(\text{LR} = 26.57\), \(p = 1.70 \times 10^{-6}\); \(6.42 \pm 1.19\) ns, \(\text{LR} = 28.70\), \(p = 5.86 \times 10^{-7}\) on the 1,093 independent epochs) and on the paired noise-subtracted excess \(M_{\rm cl} - E[M_{\rm cl}]\) halves the inter-epoch random scatter (\(\sigma_\alpha: 20.95 \to 10.60\) ns) to resolve a statistically significant modulation of \(A = 1.14 \pm 0.45\) ns (\(\text{LR} = 6.32\), \(p = 0.0425\), Step 046b), corroborated in the phase domain by the 1-parameter Thiele–Innes binary reflex velocity projection onto the \(128^\circ\) screen axis (\(2\Delta\ln L = 5.70\), 1 df \(p = 0.0170\), \(\Delta\mathrm{AIC} = -3.70\), Step 051). Multi-pulsar scaling, chromaticity, cross-telescope environmental bounds, and orbital structure provide multiscale consistency checks on the TEP framework.

The observation rejects both scalar-factorizable and non-factorizable asymmetric sub-image-blended thin-screen nulls. Unlike spatial VLBI visibility closure (which samples a single sky image at \(\mathbf{u}_{12}+\mathbf{u}_{23}+\mathbf{u}_{31}=0\)), secondary-spectrum cross-term closure isolates localized image pairs \((\tau_b-\tau_a, f_{D,b}-f_{D,a})\); even when \(M_{\rm sub}=5\) asymmetric anisotropic sub-images with refractive phase wedges blend non-factorizably inside the \(3\times 3\) pixel patch, the Hermitian overlap kernel centers the closure distribution at \(\bar{\psi}_{\rm sim}=+0.015\) rad (\(1.38\sigma\) from zero, Step 061). Across the full \(4.45\times\) wideband Parkes DR2 frequency span from 50 cm (\(710\) MHz, 216 epochs), 20 cm (\(1{,}376\) MHz, 717 epochs, \(Z=193.4\)), to 10 cm (\(3{,}165\) MHz, 458 epochs, \(Z=189.9\)), the 20-cm/10-cm phase-delay ratio is \(2.441 \pm 0.042\) (50-cm/10-cm: \(5.020 \pm 0.118\)), matching phase-fixed holonomy (\(\nu_{10}/\nu_{20}=2.299\), \(\nu_{10}/\nu_{50}=4.455\)) while excluding cold-plasma dispersive scaling (\(\propto\nu^{-2} \implies 5.29\) and \(19.85\)) at \(>65\sigma\) and refractive wedge scaling (\(\propto\nu^{-4} \implies 27.9\) and \(394\)), with the 10-cm paired delay excess (\(+8.70\pm 0.70\) ns, \(12.4\sigma\)) exceeding 20 cm (\(+4.57\pm 0.24\) ns, \(19.3\sigma\); Step 061). Within the TEP measurement model, the non-zero Phase Closure is a synchronization-holonomy candidate sourced by the non-exact disformal sector: on the measured screen-triangle geometry (Step 057), the canonical disformal connection yields \(\approx 1.5\times10^{-12}\) s as a ceiling at order-unity drift inhomogeneity, with the \(B_0=+1\) unit benchmark (\(1.3\%\) of the narrowband-equivalent \(H_{\rm implied}\approx 1.13\times10^{-10}\) s on the median leg, \(5.0\%\) on the 95th-percentile leg) and reaches \(1.14\times10^{-10}\) s (\(100.5\%\) on the median leg; \(B_0^{\rm req}=20.1\) on the \(p_{95}\) leg) at the corpus's GW170817 galactic-profile envelope \(B_0\lesssim 77.7\) (Paper 0, Step 50; while terrestrial optical-resonator bounds constrain \(B(\phi)\Sigma²\) at Earth's screened surface where \(\mathcal{S}_\Sigma² \ll 1\), decoupled from the unscreened interstellar \(B(\phi)\Sigma\,\dot{u}\) loop integral, Paper 0 §10.1), while remaining within the delay-domain bound (\(\lesssim 0.46\) ns). The geometric structure from J1603-7202 is directionally consistent with the Temporal Equivalence Principle's non-integrable time transport.

Code Availability: All data and analysis code required to reproduce the results presented in this work, including the full 15-pulsar scintillation catalog compilation, are available in the public repository at https://github.com/matthewsmawfield/TEP-J0437.

Keywords: Temporal Equivalence Principle, pulsar scintillation, synchronization holonomy, closure delays, modified gravity, interstellar medium, multi-pulsar analysis, PSR J0437-4715, PSR J1603-7202, arc curvature scaling, scattering strength

Current manuscript source: [0_abstract.html](site/components/0_abstract.html).

## Evidence Ledger

The pipeline now writes an explicit evidence ledger at `results/step_049_evidence_ledger.json`. This file separates claims by evidential status so the result can be reviewed without relying on rhetorical weight.

| Status | Claim | Current Result |
| --- | --- | --- |
| Primary evidence | J0437 rejects the additive scalar path-delay null through non-zero Phase Closure ψ | ψ̄ = 0.984 ± 0.046 rad; unweighted ψ = +1.120 rad; unweighted Rayleigh p = 1.34 × 10⁻⁴⁴; bootstrap CI [0.737, 1.235] rad excludes zero |
| Secondary delay-domain diagnostic | J0437 symmetric 10% trimmed amplitude | M_{cl,sym10} = 21.991 ± 0.75 ns; secondary robustness diagnostic, does not carry primary inferential weight |
| ISM systematic rejection | step_008 scalar-factorizable null simulations | Seven scalar-factorizable ISM null implementations; all returned zero mean Phase Closure, none reproduced the observed coherent orientation-dependent signal |
| Supporting diagnostics | Phase-scramble, pre-alignment phase closure, CMB/heliocentric frame invariance, signed bipolar cancellation | Applicable controls pass; falsification report remains phase-primary |
| Complementary geometry | J1603 bipolar decomposition and larger D/v phase dispersion | Frame-independent bipolar structure and elevated circular variance are complementary checks; phase closure remains noise-limited |
| Secondary scaling channels | Multi-pulsar scaling, cross-telescope environmental bounds, chromaticity, orbital kinematics | Directionally consistent with TEP at the precision of the current data product |
| Diagnostic only | Unsigned M_cl and delay-amplitude estimators | Folded-noise floor and delay-domain systematics prevent primary inference |

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
| **Paper 6** | [TEP-UCD](https://github.com/matthewsmawfield/TEP-UCD) | Universal Critical Density: Unifying Atomic, Galactic, and Compact Object Scales | [10.5281/zenodo.18064365](https://doi.org/10.5281/zenodo.18064365) |
| **Paper 7** | [TEP-RBH](https://github.com/matthewsmawfield/TEP-RBH) | The Soliton Wake: A Runaway Black Hole as a Gravitational Soliton | [10.5281/zenodo.18059250](https://doi.org/10.5281/zenodo.18059250) |
| **Paper 8** | [TEP-SLR](https://github.com/matthewsmawfield/TEP-SLR) | Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging | [10.5281/zenodo.18064581](https://doi.org/10.5281/zenodo.18064581) |
| **Paper 9** | [TEP-EXP](https://github.com/matthewsmawfield/TEP-EXP) | What Do Precision Tests of General Relativity Actually Measure? | [10.5281/zenodo.18109760](https://doi.org/10.5281/zenodo.18109760) |
| **Paper 10** | [TEP-COS](https://github.com/matthewsmawfield/TEP-COS) | The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars | [10.5281/zenodo.18165798](https://doi.org/10.5281/zenodo.18165798) |
| **Paper 11** | [TEP-H0](https://github.com/matthewsmawfield/TEP-H0) | The Cepheid Bias: Resolving the Hubble Tension | [10.5281/zenodo.18209702](https://doi.org/10.5281/zenodo.18209702) |
| **Paper 12** | [TEP-JWST](https://github.com/matthewsmawfield/TEP-JWST) | The Temporal Equivalence Principle: A Unified Resolution to the JWST High-Redshift Anomalies | [10.5281/zenodo.19000827](https://doi.org/10.5281/zenodo.19000827) |
| **Paper 13** | [TEP-WB](https://github.com/matthewsmawfield/TEP-WB) | The Temporal Equivalence Principle: Temporal Shear Recovery in Gaia DR3 Wide Binaries | [10.5281/zenodo.19102061](https://doi.org/10.5281/zenodo.19102061) |
| **Paper 15** | [TEP-EFA](https://github.com/matthewsmawfield/TEP-EFA) | Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly | [10.5281/zenodo.19454862](https://doi.org/10.5281/zenodo.19454862) |
| **Paper 16** | **TEP-J0437** (This repo) | Synchronization Holonomy in Pulsar Scintillation | [10.5281/zenodo.19454619](https://doi.org/10.5281/zenodo.19454619) |
| **Paper 17** | [TEP-LLR](https://github.com/matthewsmawfield/TEP-LLR) | Lunar Laser Ranging and the Nordtvedt Effect | [10.5281/zenodo.19446028](https://doi.org/10.5281/zenodo.19446028) |

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
  doi={10.5281/zenodo.19454619},
  url={https://doi.org/10.5281/zenodo.19454619},
  note={Preprint v0.6 (Sintra)}
}
```

---

## Open Science Statement

This is a working preprint shared in the spirit of open science—all manuscripts, analysis code, and data products are openly available under Creative Commons and MIT licenses to encourage and facilitate replication. Feedback and collaboration are warmly invited and welcome.

---

**Contact:** matthew@mlsmawfield.com  
**ORCID:** [0009-0003-8219-3159](https://orcid.org/0009-0003-8219-3159)
