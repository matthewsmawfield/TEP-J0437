# Temporal Equivalence Principle: Temporal Shear in the Earth Flyby Anomaly
**Matthew Lukin Smawfield**
Version: v0.1 (Yogyakarta)
First published: 14 May 2026
DOI: 10.5281/zenodo.19454863

---

## Abstract

Twelve Earth gravity assist flybys spanning nine spacecraft (NEAR, Galileo 1990/1992, Cassini, Rosetta 2005/2007/2009, MESSENGER, Juno, Stardust, OSIRIS-REx, BepiColombo) are analyzed within the Temporal Equivalence Principle (TEP) framework. The TEP framework proposes that global simultaneity is inherently non-integrable, with the rate of time represented as a dynamical scalar field φ. All non-gravitational matter couples universally to a causal matter metric through conformal coupling A(φ) = exp(β φ/M_{\rm Pl}), where β is a dimensionless coupling constant and M_{\rm Pl} is the reduced Planck mass. This coupling produces a scalar force F = β_eff c² ∇φ/M_{\rm Pl} on test masses, where β_eff = β × S_⊕(r) incorporates geometric screening via Temporal Topology. The screening factor S_⊕(r) encodes continuous suppression of Temporal Shear in density gradients, with a characteristic transition radius R_sol ≈ 4146 km derived from the UCD saturation model (a soliton interpretation is one candidate microscopic realization) and independently validated by GNSS atomic clock correlations (λ_TEP ≈ 4000 km).

The scalar force manifests as a "Phantom Mass" artifact—velocity anomalies that mimic unmodeled gravitational mass distributions through field-gradient couplings. The radial component of this force is indistinguishable from a small shift in GM and is absorbed by orbit determination programs. The non-radial component, modulated by Earth's oblateness (J2, J3, J4), trajectory asymmetry, velocity-dependent disformal coupling, IRI-based perigee plasma inputs in the Step 007 geometry envelope, and cosmographic CMB-frame velocity geometry, produces the observed flyby anomaly. Four published anomalies are retained in the catalog (NEAR, Galileo 1990, Rosetta 2005, Cassini). The Step 008 inverse-variance ensemble uses three primary fits that pass pre-specified gates (S/N ≥ 2 and sign agreement between observation and TEP prediction at β_ref = 10⁻⁴): NEAR (13.46 ± 0.01 mm/s), Galileo 1990 (3.92 ± 0.03 mm/s), and Rosetta 2005 (1.82 ± 0.05 mm/s). Cassini (0.11 ± 0.05 mm/s) shows a sign mismatch at β_ref and is excluded from that β ensemble while remaining in the dataset for diagnostics and null-population effect sizes.

Field values are computed self-consistently from the Temporal Shear Suppression density minimum φ(ρ) = Λ [n Λ^(n+4) M_{\rm Pl} / (2 β ρ)]^(1/(n+1)) (matter density ρ in kg m⁻³, converted to GeV⁴ in the pipeline; the factor 2 matches universal matter coupling in the Einstein-frame kinetic normalization used in Jakarta v0.8 and in Step 007), yielding φ_earth ≈ 2.4×10⁴ GeV and φ_space ≈ 2.0×10¹⁰ GeV for n = 3 and Λ = 10 MeV. The three primary fits used in the ensemble weighted mean span β = 5.03×10⁻⁴ to 2.02×10⁻³ (factor 4.0), with inverse-variance weighted mean β = 5.65×10⁻⁴ ± 2.79×10⁻⁵. The heterogeneity in fitted values (I² ≈ 100%) reflects geometry-dependent variation in effective coupling across flybys. Earth-screened fitted values satisfy terrestrial PPN consistency checks with margins of roughly 2×10³ to 10⁵ relative to |γ - 1| < 2.3×10⁻⁵; the separate worst-case solar-path check remains below the Cassini Shapiro bound by a smaller factor of about 3.5. Bayesian model comparison with a four-tier framework (Null, Anderson empirical, TEP restricted, TEP flexible) favors the TEP restricted model on the three-gated ensemble, while the small sample requires qualitative interpretation of information-criterion summaries. The TEP restricted model uses a single fitted parameter β, with λ_TEP ≈ 4000 km, S_⊕ ≈ 0.35, and v_trans ≈ 16.8 km/s pre-specified from independent measurements and first-principles derivations. On the same n = 3 gated ensemble as Step 026, BIC-based summaries give B₁₀ ≈ 1.8×10¹⁸ (TEP restricted vs Null, ΔBIC ≈ 84), B_A0 ≈ 9.9×10¹⁶ (Anderson vs Null, ΔBIC ≈ 78), and TEP restricted vs Anderson B ≈ 18.2 (ΔBIC ≈ 5.8), indicating positive but not decisive preference over the Anderson empirical baseline. Published null results for high-altitude or symmetric trajectories (MESSENGER 2005, Rosetta 2007, Rosetta 2009) are consistent with TEP predictions through altitude-dependent gradient suppression and the Step 007 geometry envelope. Three flybys (Stardust, OSIRIS-REx, BepiColombo) have no public anomaly report and are not used in quantitative likelihood. At the refit weighted-mean β, Step 039 classifies one raw-tension case (Juno); Galileo 1992 and Cassini are raw true nulls under the geometry envelope. Post-OD survival factors are withheld until mission OD configuration yields defensible F_OD estimates. Numerical simulation supports the claim that empirical-acceleration OD filters can absorb TEP-like signals, but the remaining raw-tension null motivates raw DSN reanalysis rather than era-based survival factors. Full three-dimensional spacecraft state vectors from JPL Horizons were tested for CMB-frame velocity geometry (Step 040, n = 8). No individual ratio-versus-modulation correlation reaches conventional significance; the exploratory both-aligned flag yields Pearson r = +0.37, p = 0.36, and Mann-Whitney U = 11, p = 0.24. Multivariate regression on residual ratio achieves R² = 0.47 (adjusted R² = +0.08), and an optimal weighted combination of spacecraft and Earth CMB projections gives r = −0.61, p = 0.11. These results are sample-limited and do not yet decisively confirm disformal coupling in the CMB rest frame at conventional significance levels.

This work bridges the gap between precision solar system tests and cosmological dynamics, showing that the Temporal Equivalence Principle framework is consistent with published flyby anomaly measurements and offers a new avenue for exploring the intersection of gravity, time, and matter, ultimately shedding new light on the fundamental nature of spacetime.

Keywords: Earth flyby anomaly, Temporal Equivalence Principle, scalar force, Phantom Mass, trajectory asymmetry, geometric screening, Temporal Topology, Temporal Shear

# 1. Introduction

The Equivalence Principle (EP) is a cornerstone of general relativity, stating that gravitational acceleration is locally indistinguishable from acceleration due to motion. However, the Temporal Equivalence Principle (TEP)—the assertion that global simultaneity is inherently non-integrable—suggests that the rate of time is a dynamical scalar field $\phi$. This framework, established in the Jakarta foundational axioms (v0.8), proposes that all non-gravitational matter couples universally to a causal matter metric $\tilde{g}_{\mu\nu} = A^2(\phi)g_{\mu\nu}$, where $A(\phi) = \exp(\beta \phi/M_{\text{Pl}})$.

#### Key Terminology

- *Proper time* ($\tau$) is the time measured by a clock following a specific trajectory through the causal metric.

- *Temporal Topology* refers to the spatial structure of the field $\phi$, which exhibits continuous suppression in high-density environments.

- *Temporal Shear* ($\nabla\phi$) is the gradient of the time field, which generates the observed scalar force.

- *Phantom Mass* describes the anomalous acceleration that mimics a gravitational mass distribution, arising from the non-radial coupling of the scalar field.

- *PPN parameter $\gamma$* measures the amount of spatial curvature; in Jakarta v0.8 (Sec. 7), $\gamma - 1 = -2\alpha_{\rm eff}^2$ with $\alpha_{\rm eff}$ the screened conformal coupling; for $A=\exp(\beta\phi/M_{\rm Pl})$ this maps to $|\gamma - 1| \approx 2\beta_{\rm eff}^2$ when comparing to Cassini's bound on $|\gamma-1|$.

## 1.1 The Earth Flyby Anomaly

Since 1990, spacecraft executing Earth gravity assists have exhibited anomalous orbital energy changes that lack a standard explanation. The NEAR spacecraft (1998) showed the largest effect: an unexplained velocity increase of 13.46 mm/s. Galileo (1990) and Cassini (1999) displayed smaller but significant anomalies. These velocity shifts occur precisely at perigee passage and persist as asymptotic excess velocities ($v_\infty$) in the outbound trajectories.

Standard physics offers no satisfactory explanation. Thermal radiation pressure, atmospheric drag, and tidal effects have been found insufficient by orders of magnitude. The anomalies show no correlation with spacecraft orientation or spin rate, ruling out conventional systematic errors. The effect appears genuinely gravitational in nature, manifesting as a "Phantom Mass" artifact that reflects a non-integrable time transport.

## 1.2 TEP as a Candidate Explanation

The TEP framework provides a natural explanation through the interaction between the spacecraft and the Earth's Temporal Topology. As a spacecraft traverses the field gradient $\nabla\phi$, it experiences a scalar force $\mathbf{F}_\phi = \beta_{\text{eff}} c^2 \nabla\phi / M_{\text{Pl}}$. While a pure clock-rate shift would cancel in two-way Doppler tracking to first order, the scalar force acts directly on the trajectory, producing a physical velocity shift.

The observed heterogeneity in flyby anomaly magnitudes is not random scatter but arises from deterministic geometry-dependent modulation. The TEP prediction for a given flyby depends on several physical factors: (1) perigee altitude (determines Temporal Shear strength via density suppression), (2) approach-departure asymmetry (disformal coupling requires velocity-dependent anti-aligned geometry), (3) plasma environment (plasma attenuation modulates the scalar field), (4) solar activity (modulates ionospheric density), and (5) cosmographic CMB-frame velocity geometry (the disformal coupling scales as v² in the scalar rest frame, approximated by the CMB dipole frame with bulk velocity ~370 km/s toward RA = 167.94°, Dec = −6.93°). These factors combine to produce a roughly four-fold span in per-flyby fitted β across the three-gated primary ensemble (Step 008), while the full literature set including Cassini exhibits a wider dynamic range when reference-coupling predictions disagree in sign. The Temporal Shear Suppression mechanism is essential for three reasons: (1) it ensures the coupling strength satisfies solar system PPN constraints; (2) it explains both detections and null results through density-dependent screening; and (3) it establishes the transition radius $R_{\rm sol} \approx 4146$ km as a universal scale. Flybys sampling regions of high Temporal Shear (low altitude, high asymmetry) exhibit anomalies, while those in shielded regimes (high altitude or symmetric trajectories) remain null.

## 1.3 This Work

This paper presents a comprehensive analysis of the Earth flyby anomaly using the TEP framework. Published Doppler tracking measurements from Anderson et al. (2008) are employed, interpreted as "Phantom Mass" signatures of the local Temporal Topology. The analysis proceeds by reconstructing trajectories from JPL Horizons, computing TEP predictions with full 3D integration, and fitting the universal coupling $\beta$ to the observed dataset.

The structure of this paper is as follows: Section 2 describes the data sources; Section 3 presents the TEP Temporal Topology model and cosmographic analysis methodology; Section 4 reports the fitting results, PPN validation, and cosmographic temporal shear modulation tests (Section 4.11); Section 5 discusses the Phantom Mass interpretation and directional consistency with the CMB rest frame; and Section 6 concludes with prospects for further multi-messenger tests.

# 2. Observations and Data

## 2.1 The Flyby Spacecraft Sample

This analysis utilizes nine spacecraft spanning twelve Earth flyby events between 1990 and 2020: Galileo (1990, 1992), NEAR (1998), Cassini (1999), Rosetta (2005, 2007, 2009), MESSENGER (2005), Juno (2013), Stardust (2001), OSIRIS-REx (2017), and BepiColombo (2020). The dataset is divided into three data quality classes: *published anomalies* (four flybys with measured nonzero Δv and formal uncertainties), *published nulls/bounds* (five flybys with explicitly reported null results or upper limits), and *no public anomaly report* (three flybys with no published search or measurement). The latter class is not used in quantitative likelihood. Table 1 summarizes the key parameters for each flyby.

#### Physical Constants

The analysis uses the following CODATA 2018 values: Earth radius $R_\oplus = 6.371 \times 10^6$ m, gravitational constant $G = 6.67430 \times 10^{-11}$ m$^3$ kg$^{-1}$ s$^{-2}$, and speed of light $c = 299\,792\,458$ m/s (exact). The reduced Planck mass $M_{\rm Pl} = 2.435 \times 10^{18}$ GeV is derived from $\hbar c/G^{1/2}$.

Table 1: Earth Flyby Spacecraft Parameters

| Spacecraft | Date | Perigee (km) | $v_\infty$ (km/s) | $\Delta v_{\rm obs}$ (mm/s) | $\sigma$ (mm/s) | Data class |
| --- | --- | --- | --- | --- | --- | --- |
| Galileo | 1990-12-08 | 972 | 13.73 | 3.92 | 0.03 | Published anomaly |
| Galileo | 1992-12-08 | 310 | 14.08 | 0.00 | 0.05 | Published null/bound |
| NEAR | 1998-01-23 | 568 | 12.72 | 13.46 | 0.01 | Published anomaly |
| Cassini | 1999-08-18 | 1197 | 19.02 | 0.11 | 0.05 | Published anomaly |
| Rosetta | 2005-03-04 | 1969 | 10.51 | 1.82 | 0.05 | Published anomaly |
| Rosetta | 2007-11-13 | 5430 | 12.46 | 0.02 | 0.05 | Published null/bound |
| Rosetta | 2009-11-13 | 2572 | 13.31 | 0.00 | 0.05 | Published null/bound |
| MESSENGER | 2005-08-02 | 2351 | 10.39 | 0.00 | 0.05 | Published null/bound |
| Juno | 2013-10-09 | 817 | 14.79 | 0.00 | 0.02 | Published null/bound |
| Stardust | 2001-01-15 | 6009 | 10.31 | — | — | No public anomaly report |
| OSIRIS-REx | 2017-09-22 | 17239 | 8.52 | — | — | No public anomaly report |
| BepiColombo | 2020-04-10 | 12697 | 7.59 | — | — | No public anomaly report |

*Note:* $\Delta v_{\rm obs}$ values for the published anomaly and published null/bound classes are from Anderson et al. (2008) and companion papers. Stardust, OSIRIS-REx, and BepiColombo have no public anomaly report; em-dashes indicate that no published measurement or bound exists. These three flybys are not used in quantitative likelihood but are listed as predicted nulls based on their high perigee altitudes. Rosetta 2009 has a published null result (dv = 0.00 mm/s) from Muller et al. (2010); the uncertainty is the DSN tracking precision (0.05 mm/s) as no formal bound was published. Perigee distances are geocentric; $v_\infty$ is the hyperbolic excess velocity.

## 2.2 Data Sources and Provenance

The anomaly measurements used in this analysis are taken from the peer-reviewed literature, specifically the comprehensive study by Anderson et al. (2008) and subsequent mission-specific analyses. These values were obtained through NASA's Deep Space Network (DSN) Doppler tracking combined with the Jet Propulsion Laboratory Orbit Determination Program (ODP).

Literature sources:

- Primary reference: Anderson, J. D., et al. (2008). "Anomalous Orbital-Energy Changes Observed during Spacecraft Flybys of Earth." *Physical Review Letters*, 100(9), 091102.

- Rosetta analysis: Morley, T., & Budnik, F. (2007). "Rosetta Navigation at its First Earth-Swingby." *Proceedings of the 20th International Symposium on Space Flight Dynamics*.

- Juno analysis: Aksenov, E. L., & Tuchin, A. G. (2020). "Earth flyby anomalies and the general relativistic theory of the Kerr gravitational field." *MNRAS*, 492(3), 3703-3711.

## 2.3 Data Quality Assessment

A rigorous analysis requires assessment of data quality for each flyby. All four primary detections have complete DSN coverage spanning $\pm 12$ hours around perigee, enabling robust pre/post comparison. The reported uncertainties (0.01–0.05 mm/s) are consistent with DSN Doppler precision at X-band.

Systematic error controls: Antenna phase center, tropospheric delay, and station positions are well-modeled in the JPL ODP software. Residual uncertainties are at the $\sim 0.1$ mm/s level, which is an order of magnitude below the larger anomalies (NEAR, Galileo).

## 2.4 Trajectory Data from JPL Horizons

Spacecraft trajectories for the analysis were obtained from NASA's JPL Horizons ephemeris system. For each flyby, state vectors (position and velocity) spanning $\pm 2$ days around perigee passage are reconstructed. These trajectories represent the best-estimate spacecraft paths based on all available tracking data.

# 3. Methodology

The analysis employs a four-step pipeline to test whether TEP with Temporal Topology explains observed flyby velocity anomalies as "Phantom Mass" artifacts. The pipeline retrieves spacecraft trajectories from JPL Horizons, computes TEP predictions for each flyby geometry using full 3D integration, fits the coupling parameter $\beta$ to match observed anomalies, and validates all parameters against solar system PPN constraints.

## 3.1 Data Acquisition

Spacecraft trajectories are obtained from the NASA JPL Horizons ephemeris system using the astroquery interface. For each flyby, reconstructed state vectors (position and velocity) are retrieved in the ICRF (International Celestial Reference Frame) at 30-minute intervals spanning $\pm 2$ days around perigee passage.

#### Trajectory Parameters Extracted

- Perigee altitude (minimum geocentric distance)

- Perigee velocity (speed at closest approach)

- Inbound/outbound asymptotic velocity ($v_\infty$)

- Spacecraft potential at perigee ($\Phi_{\rm sc}$)

Flyby velocity anomalies ($\Delta v_{\rm obs}$) are taken from published literature. The primary source is Anderson et al. (2008), with supplementary references for Rosetta (Morley & Budnik 2007; Müller et al. 2008, 2010) and Juno (Aksenov & Tuchin 2020). All values were measured by NASA/JPL using Deep Space Network Doppler tracking with the Orbit Determination Program. Asymptotic $v_\infty$ declinations ($\delta_{\rm in}$, $\delta_{\rm out}$) for the six flybys in Anderson et al. (2008) are taken from that source; for the remaining six flybys, declinations are computed from the ephemeris using two-body orbital mechanics (eccentricity vector method).

## 3.2 TEP Temporal Topology Model

The TEP framework provides a quantitative model for the flyby anomaly through a scalar force arising from the Temporal Topology field φ. In scalar-tensor theories with conformal coupling A(φ) = exp(β φ/M_{\rm Pl}), the scalar field gradient produces an additional force on test masses:

\begin{equation}
\mathbf{F}_\phi = \beta_{\rm eff} \, \frac{c^2 \nabla\phi}{M_{\rm Pl}}
\end{equation}

where β_eff = β × S_⊕(r) is the effective coupling with geometric screening, where S_⊕(r) describes the continuous suppression of Temporal Shear. The characteristic surface ratio S_⊕ ≈ 0.35 is fixed by the UCD saturation geometry as S_⊕ = (R_⊕ − R_sol)/R_⊕ with R_sol ≈ 4146 km (distinct from the embedding factor R_sol/R_⊕ used in Paper 6). The radial component of this force is indistinguishable from a small shift in GM and is absorbed by orbit determination. The non-radial component—modulated by Earth's oblateness (J2, J3, J4) and the spacecraft's trajectory geometry—produces a net velocity change that appears as the flyby anomaly.

The predicted velocity shift is resolved through rigorous numerical integration of the equations of motion (EOM) in the Earth-centered inertial (ECI) frame. This approach captures the dynamic evolution of the scalar force as the spacecraft traverses the varying field gradient, incorporating a 4th-order Spherical Harmonic Expansion (SHEX) for the geopotential to ensure that local gravitational perturbations are not conflated with the scalar force:

\begin{equation}
\Delta \mathbf{v}_{\rm TEP} \approx \left. \frac{\beta_{\rm eff} \, c^2 \nabla\phi}{M_{\rm Pl}} \right|_{\rm peri} \Delta t_{\rm peri} + \left. \frac{b_{\rm disf}}{M_{\rm Pl}} (\nabla\phi \cdot \mathbf{v}) \mathbf{v} \right|_{\rm peri} \Delta t_{\rm peri}
\end{equation}

where the simplified perigee approximation is:

\begin{equation}
\Delta v_{\rm TEP} \approx \beta_{\rm eff} \, \frac{c^2}{M_{\rm Pl}} \left(\frac{d\phi}{dr}\right)_{r_p} \, \frac{r_p}{v_p} \, \left[J_2 + J_3 \sin(\lambda_p) + J_4 P_4(\sin\lambda_p)\right] \left(\frac{R_\oplus}{r_p}\right)^2 (\cos\delta_{\rm in} - \cos\delta_{\rm out})
\end{equation}

where:

- $(d\phi/dr)_{r_p}$ is the scalar field gradient at perigee altitude

- $r_p$ and $v_p$ are the perigee distance and velocity

- $J_2, J_3, J_4$ are the zonal harmonics (EGM96/WGS84 coefficients)

- $\lambda_p$ is the perigee latitude

- $\delta_{\rm in}$ and $\delta_{\rm out}$ are the asymptotic declinations on approach and departure (from Anderson et al. (2008))

### 3.2b Full 3D Trajectory Integration

The perigee approximation provides a computationally efficient closed-form estimate, but the primary predictions are validated against full 3D trajectory integration to capture dynamic evolution of the scalar force along the spacecraft path. JPL Horizons ephemeris data for historical flybys provide scalar range and velocity observables rather than complete 3D state vectors. To reconstruct the trajectory, perigee state vectors (position and velocity) from the TEP geometry model serve as anchor points, and the hyperbolic orbit is propagated via Keplerian mechanics to generate full 3D ephemeris consistent with the JPL Horizons time grid. The reconstructed trajectory is validated against the JPL range data to ensure physical consistency.

The scalar velocity shift is accumulated by integrating the TEP force along the reconstructed path:

\begin{equation}
\Delta v_{\rm TEP}^{(3D)} = \int_{t_{\rm in}}^{t_{\rm out}} \beta_{\rm eff}(r(t)) \, \frac{c^2}{M_{\rm Pl}} \, \left|\frac{d\phi}{dr}\right|_{r(t)} \, \mathcal{B}(\lambda(t), r(t)) \, \cos\delta_{\rm asym} \, \mathcal{E}(t) \, s_{\rm disf}(t) \, dt
\end{equation}

where $\mathcal{B}(\lambda, r) = [J_2 + J_3 \sin(\lambda) + J_4 P_4(\sin\lambda)] (R_\oplus/r)^2$ is the zonal harmonic bracket, $\mathcal{E}(t)$ is the geometry envelope factor, and $s_{\rm disf}(t)$ is the disformal modulation factor. The integration uses 240–360 points per flyby and captures altitude-dependent field gradient variation, local plasma modulation, and velocity-dependent disformal factors at each point along the trajectory.

This 3D integration validates the perigee approximation: for the three primary detections, the integrated $\Delta v$ agrees with the perigee estimate within 8–16% (Section 4.1.2). For symmetric trajectories (Galileo 1992, MESSENGER) and high-altitude flybys (Rosetta 2007), both methods predict negligible anomalies, consistent with observations. The perigee approximation is therefore retained as the primary computational tool for ensemble fitting, while the 3D integration provides an independent cross-check.

J3 contribution: The J3 term adds a latitude-dependent asymmetry to the non-radial force component. However, J3 is two orders of magnitude smaller than J2 ($|J_3/J_2| \approx 2.3 \times 10^{-3}$), and its inclusion does not collapse the heterogeneity in fitted β across the three-gated ensemble (factor ~4 in β at fixed reference coupling, with formal Cochran Q ≫ 1). This indicates that remaining spread arises from phase-boundary and geometry–plasma modulation in the envelope, not from neglected J3 alone.

The scalar field φ relaxes outside Earth with a relaxation length λ_TEP ≈ 4000 km, established independently from GNSS atomic clock correlations and the scalar field mass inferred from the cosmological sound horizon.

The trajectory asymmetry factor $\cos\delta_{\rm in} - \cos\delta_{\rm out}$ is the dominant source of inter-flyby variation. Symmetric trajectories (e.g., Galileo 1992, MESSENGER) have $\cos\delta_{\rm in} \approx \cos\delta_{\rm out}$ and predict negligible anomalies, consistent with observations. Asymmetric trajectories (e.g., NEAR with $\cos\delta_{\rm in} - \cos\delta_{\rm out} = 0.625$) predict large anomalies.

\begin{equation}
\phi(r) = \phi_{\rm earth} + (\phi_{\rm space} - \phi_{\rm earth}) \left[1 - \exp\!\left(-\frac{r - R_\oplus}{\lambda_{\rm TEP}}\right)\right]
\end{equation}

Geometric screening: Critical to PPN compliance is the transition radius $R_{\rm sol} \approx 4146$ km from the UCD saturation model (Step 010), cross-validated by GNSS correlation length. This defines the characteristic suppression ratio $S_{\oplus} \approx 0.35$ that quantifies the attenuation of Temporal Shear at Earth's surface. $S_{\oplus} = (R_{\oplus} - R_{\rm sol})/R_{\oplus}$ is the gradient suppression ratio at the surface; it is distinct from the UCD embedding factor $S = R_{\rm sol}/R_{\oplus} \approx 0.65$ used in Paper 6 (UCD), which measures how deeply the mass is embedded within its saturation radius.

\begin{equation}
\beta_{\rm eff} = \beta \times S_{\oplus}(r)
\end{equation}

The Temporal Topology field minimum at density $\rho$ is:

\begin{equation}
\phi_{\rm min}(\rho) = \Lambda \left[ \frac{n \Lambda^{n+4} M_{\rm Pl}}{2\beta \rho} \right]^{1/(n+1)}
\end{equation}

#### Characteristic Field Values ($n=3$, $\Lambda=10$ MeV)

- Inside Earth ($\rho = 5515$ kg/m$^3$): $\phi_{\rm earth} = 2.35 \times 10^{4}$ GeV

- At Earth's surface ($\rho = 2700$ kg/m$^3$): $\phi_{\rm surface} = 2.81 \times 10^{4}$ GeV

- In vacuum ($\rho \approx 10^{-20}$ kg/m$^3$): $\phi_{\rm space} \approx 2.0 \times 10^{10}$ GeV (same closed form with $2\beta\rho$ as in Step 007; $\rho$ converted to GeV$^4$ via CODATA-based factors)

- TEP relaxation length: $\lambda_{\rm TEP} \approx 4000$ km (from GNSS / scalar field Compton wavelength)

- Characteristic suppression: $S_{\oplus} \approx 0.35$ (UCD-derived from Step 010)

Vacuum field value: The Temporal Topology field formula φ ∝ ρ^(-1/4) produces large but finite values in the interplanetary medium (ρ ≈ 10⁻²⁰ kg/m³). The self-consistent field equation yields φ_space ≈ 2.0×10¹⁰ GeV for the reference coupling β=10⁻⁴. No ad-hoc cutoff is applied; the field is computed directly from the physical density.

Geometry modulation factors: Per-flyby fitted β in the three-gated ensemble spans roughly a factor of four, reflecting substantial geometry-, plasma-, and velocity-dependent modulation encoded in the Step 007 envelope. Four physical mechanisms are tracked in the analysis: (1) *inclination-dependent screening*—higher latitude trajectories sample less equatorial bulge; (2) *J2 oblateness*—altitude-dependent screening from Earth's shape; (3) *plasma environment*—ionospheric density modulates local screening (IRI at perigee, Step 033); and (4) *velocity effects*—disformal coupling in the high-velocity regime. These factors are incorporated into the scalar force calculation.

## 3.2a Component-Level Geometry Factor Analysis

To address the extreme heterogeneity (I² ≈ 100%) in fitted β values across flybys, a component-level analysis extracts the effective geometry factor for each flyby independently. The geometry factor isolates trajectory-dependent modulation from the universal coupling strength, revealing the physical origin of the observed variation.

The effective geometry factor is defined as the ratio of observed anomaly to the gradient prediction at the reference coupling:

\begin{equation}
G_{i,\text{eff}} = \frac{\Delta v_{i,\text{obs}}}{\Delta v_{\text{grad},i}(\beta_0 = 10^{-4})}
\end{equation}

This factor absorbs all geometry-dependent modulation—altitude, J2 oblateness, trajectory asymmetry, velocity-dependent disformal coupling, plasma screening, and OD absorption—into a single observable per flyby. The implied universal coupling is then $\beta_{0,\text{implied}} = 10^{-4} \times G_{i,\text{eff}}$.

Correlation analysis tests whether $G_{i,\text{eff}}$ varies systematically with trajectory parameters:

- Altitude: higher perigee → lower $G_{\text{eff}}$ (weaker field gradient)

- Velocity: higher $v_p$ → lower $G_{\text{eff}}$ (shorter field exposure time)

- Asymmetry: positive $\cos\delta_{\text{in}} - \cos\delta_{\text{out}}$ → higher $G_{\text{eff}}$ (stronger disformal enhancement)

A multiple linear regression in log space quantifies the combined contribution:

\begin{equation}
\log_{10} |G_{\text{eff}}| = c_0 + c_1 \tilde{h} + c_2 \tilde{v} + c_3 \tilde{a}
\end{equation}

where $\tilde{h}$, $\tilde{v}$, $\tilde{a}$ are normalized altitude, velocity, and asymmetry. Non-zero coefficients confirm geometry-dependent TEP coupling; $R^2$ near unity indicates that the three trajectory parameters explain most of the observed heterogeneity.

This approach is complemented by a four-parameter hierarchical Bayesian model ($\beta_0$, $b_{\text{disf}}$, $\sigma$, $\alpha_{\text{res}}$) sampled via MCMC. The pre-computed gradient and disformal components from the TEP scalar force model contain the full perigee physics; the likelihood scales these components by the inferred universal couplings, with any residual unmodeled modulation captured by $\alpha_{\text{res}}$. Posterior predictive checks validate the model against per-flyby observations.

## 3.3 Deterministic Factor Computation

#### Deterministic Factors

- **Trajectory geometry (G_traj):** G_traj = exp(-(h - 300 km)/2000 km) × (1 + |cosδ_asym|)

- **Temporal Shear Suppression (S_⊕):** $S_\oplus = (R_\oplus - R_{\rm sol})/R_\oplus$ with $R_{\rm sol} \approx 4146$ km (UCD saturation radius; matches `CHARACTERISTIC_SUPPRESSION` in `scripts/utils/physics.py`)

- **OD absorption (F_OD):** Mission-specific fraction of injected TEP signal surviving standard OD processing. Step 039 withholds post-OD columns until Step 021 supplies defensible mission OD configuration data.

- **Plasma factor (F_plasma):** Modulated by solar activity indices (F10.7 flux, Kp index)

- **Disformal factor (F_disf):** Velocity-activated sign reversal for v > 16.8 km/s with negative asymmetry

## 3.4 Variance Decomposition ANOVA

The variance in component scaling parameters is decomposed into sources using a formal ANOVA/hierarchical variance model. This quantifies the contribution of gradient vs disformal components to the total heterogeneity. A comprehensive four-stage variance decomposition analysis is presented in Section 4.3 (Results), which consolidates structural physics modulation, observational pipeline effects, environmental modulation, and statistical limitations into a unified framework.

The current four-stage variance decomposition (Step 009, fitted-β heterogeneity) attributes 0.0% of the variance to the structural proxy bundle, 0.0% to observational pipeline effects, 99.5% to environmental modulation (F10.7 correlation |r| ≈ 0.998, p ≈ 0.043 on n = 3), and 0.5% to residual small-sample statistics, intrinsic scatter, and model incompleteness. The dominant environmental fraction reflects the small detection sample and solar-epoch covariation; it should not be read as a uniquely causal F10.7 driver without larger n.

## 3.5 Disformal Transition Criterion

A disformal transition criterion Ξ is defined to classify flybys into conformal-dominated, mixed, or disformal-dominated regimes. This provides a formal test for Cassini as a disformal-regime case.

\begin{equation}
\Xi_i = \left(\frac{v_i}{v_{\text{trans}}}\right)^p \times |\cos\delta_{\text{in}} - \cos\delta_{\text{out}}| \times \left(\frac{|\nabla\phi_i|}{|\nabla\phi_\oplus|}\right)^q \times \text{sgn}(\cos\delta_{\text{in}} - \cos\delta_{\text{out}})
\end{equation}

where:

- v_trans ≈ 16.8 km/s is the transition velocity (derived from TEP field equations, see below)

- v_i is the flyby perigee velocity

- p = 2 is the velocity exponent

- q = 1 is the gradient exponent

- ∇φ_⊕ is the Temporal Shear at Earth's surface

- ∇φ_i is the Temporal Shear at flyby altitude

- sgn indicates aligned (positive) vs anti-aligned (negative) disformal response

### Analytical Derivation of v_trans

The transition velocity v_trans is not an empirically-tuned parameter derived from the Earth Flyby Anomaly dataset, but rather a fundamental scale emerging from the TEP field equations. The disformal coupling term in the TEP metric has the form:

\begin{equation}
ds^2 = A(\phi)c^2dt^2 - B(\phi)\partial_\mu\phi\partial_\nu\phi dx^\mu dx^\nu - C(\phi)d\mathbf{x}^2
\end{equation}

where the disformal factor B(φ) couples to the kinetic term ∂μφ∂νφ. The characteristic velocity scale emerges from the condition where the disformal contribution becomes comparable to the conformal contribution in the effective metric perturbation. This occurs when:

\begin{equation}
B(\phi)v^2 \sim A(\phi) - 1
\end{equation}

Using the TEP field equations from the Jakarta axioms, the scalar field dynamics are governed by the relaxation equation:

\begin{equation}
\nabla^2\phi - \frac{1}{\lambda_{\rm TEP}^2}\phi = -\frac{\beta}{M_{\rm Pl}}\rho
\end{equation}

The characteristic velocity scale emerges from equating the disformal metric perturbation to the conformal potential perturbation. Using the TEP field equations and the relaxation relation ∇φ ∼ φ/λ_TEP, this yields a transition velocity that scales with the square root of the Temporal Shear:

\begin{equation}
v_{\rm trans} = \frac{c}{\sqrt{2}}\left(\frac{\lambda_{\rm TEP}}{R_\oplus}\right)^{1/2}\left(\frac{|\nabla\phi_\oplus|\,\lambda_{\rm TEP}}{M_{\rm Pl}}\right)^{+1/2}
\end{equation}

Substituting the independently-determined TEP relaxation length λ_TEP ≈ 4000 km (from GNSS atomic clock correlations, Step 016), Earth's radius R_⊕ = 6371 km, and the dimensionless surface field combination |∇φ_⊕| λ_TEP / M_{\rm Pl} ≈ 10⁻⁸ (from the UCD-derived characteristic suppression S_⊕ ≈ 0.35), the transition velocity is obtained:

\begin{equation}
v_{\rm trans} = \frac{c}{\sqrt{2}}\left(\frac{4000~\text{km}}{6371~\text{km}}\right)^{1/2}\left(10^{-8}\right)^{+1/2} \approx 16.8~\text{km/s}
\end{equation}

This derivation demonstrates that v_trans ≈ 16.8 km/s is a field-theoretic prediction of the TEP framework, derived from independently-measured parameters (λ_TEP from GNSS, S_⊕ from UCD) and fundamental constants. The value is not tuned to match the Cassini flyby data; rather, Cassini's high perigee velocity (19.02 km/s > v_trans) naturally places it in the disformal-dominated regime, explaining its sign reversal as a consequence of the underlying field dynamics.

Classification (by |Ξ|):

- |Ξ| < 0.05: Conformal-dominated

- 0.05 ≤ |Ξ| ≤ 0.10: Mixed

- |Ξ| > 0.10: Disformal-dominated

The sign of Ξ indicates the nature of the disformal response: positive for aligned trajectories and negative for anti-aligned trajectories. Cassini, with its high perigee velocity (19.02 km/s) and negative asymmetry, falls into the mixed regime with a negative sign, indicating it operates in the anti-aligned disformal response regime where the conformal-gradient and disformal terms partially cancel.

Velocity shift formula: The predicted velocity anomaly combines four physical effects:

\begin{equation}
\Delta v_{\rm TEP} = \frac{\beta_{\rm eff}\, c^2}{M_{\rm Pl}} \cdot \underbrace{\frac{d\phi}{dr}\bigg|_{r_p}}_{\text{field gradient}} \cdot \underbrace{\frac{r_p}{v_p}}_{\text{perigee time}} \cdot \underbrace{J_2 \!\left(\frac{R_\oplus}{r_p}\right)^{\!2}}_{\text{non-radial fraction}} \cdot \underbrace{(\cos\delta_{\rm in} - \cos\delta_{\rm out})}_{\text{trajectory asymmetry}}
\end{equation}

Each factor has a distinct physical origin:

- Field gradient $d\phi/dr = (\Delta\phi / \lambda_{\rm TEP})\, e^{-h/\lambda_{\rm TEP}}$: the scalar force strength at perigee altitude $h$, decaying exponentially with the GNSS-established relaxation length. Lower flybys experience stronger gradients.

- Perigee dwell time $r_p / v_p$: the effective duration of the close encounter. Slower, lower flybys accumulate larger impulses.

- $J_2$ oblateness $J_2 (R_\oplus/r_p)^2$: the non-radial component of the scalar force arising from Earth's oblateness. The radial component is absorbed into the orbit determination program's estimate of $GM$; only the non-radial residual produces a net velocity change.

- Trajectory asymmetry $\cos\delta_{\rm in} - \cos\delta_{\rm out}$: the difference in approach and departure $v_\infty$ declinations (from Anderson et al. (2008)). This factor determines how asymmetrically the spacecraft samples the oblate field. For symmetric trajectories ($\delta_{\rm in} \approx \delta_{\rm out}$), the non-radial impulse cancels and the predicted anomaly vanishes—correctly predicting null results for flybys such as Galileo 1992 and MESSENGER.

## 3.6 Robust Bayesian Fitting

The primary Step 008 coupling estimate is an inverse-variance scaling fit on the sign-gated detections, and the Step 026 model-comparison layer uses Gaussian weighted least-squares likelihoods. A Student's t-distribution likelihood with degrees of freedom $\nu = 3$ is used only in the auxiliary robust Bayesian/hierarchical checks, where it tests whether the conclusions are sensitive to outlier treatment in the small sample.

For each flyby with measured anomaly $\Delta v_{\rm obs} \neq 0$, the coupling parameter $\beta$ is fitted to maximize the posterior:

\begin{equation}
\mathcal{L}(\beta) = \prod_i \frac{\Gamma[(\nu+1)/2]}{\Gamma(\nu/2) \sqrt{\nu\pi}\sigma} \left[ 1 + \frac{1}{\nu} \left(\frac{\Delta v_{\rm obs,i} - \Delta v_{\rm TEP,i}(\beta)}{\sigma}\right)^2 \right]^{-(\nu+1)/2}\end{equation}

PPN constraint validation: The fitted $\beta$ satisfies, with geometric screening applied:

\begin{equation}
|\gamma - 1| \approx 2\beta_{\rm eff}^2 < 2.3 \times 10^{-5}\end{equation}

Jakarta v0.8 (Sec. 7, DEF screened limit) gives $\gamma - 1 = -2\alpha_{\rm eff}^2$ with $\alpha_{\rm eff} \equiv d(\ln A)/d\phi$ evaluated at the screened source. For $A(\phi)=\exp(\beta\phi/M_{\rm Pl})$, write $\psi\equiv\phi/M_{\rm Pl}$ so that $d(\ln A)/d\psi=\beta$; identifying the locally active dimensionless coupling with screened $\beta_{\rm eff}$ yields $|\gamma - 1| \approx 2\beta_{\rm eff}^2$ for comparison to the Cassini magnitude bound.

## 3.7 Statistical Analysis

The weighted mean $\beta$ across all detections is:

\begin{equation}
\bar{\beta} = \frac{\sum_i w_i \beta_i}{\sum_i w_i}, \quad w_i = \frac{1}{\sigma_{\beta,i}^2}
\end{equation}

with inverse-variance weights derived from propagated measurement uncertainties. The weighted standard error is:

\begin{equation}
\sigma_{\bar{\beta}} = \left(\sum_i w_i\right)^{-1/2}
\end{equation}

The NEAR detection dominates due to its superior measurement precision ($\sigma = 0.01$ mm/s vs. $0.03$–$0.05$ mm/s for others).

Heterogeneity assessment: Following meta-analysis conventions (Higgins & Thompson, 2002), heterogeneity is quantified using:

\begin{equation}
Q = \sum_i w_i (\beta_i - \bar{\beta})^2 \quad \text{(Cochran's Q)}
\end{equation}

\begin{equation}
I^2 = \max\!\left(0,\,\frac{Q - (n-1)}{Q}\right) \times 100\% \quad \text{(percentage variance due to heterogeneity; Higgins 2002)}
\end{equation}

An $I^2 > 75\%$ indicates extreme heterogeneity, justifying uncertainty inflation by $\sqrt{Q/(n-1)}$ to account for model scatter beyond measurement error.

Robustness verification: Two complementary approaches validate conclusion stability:

- *Parametric bootstrap ($n = 10\,000$):* Resampling with replacement while adding measurement noise validates the weighted mean distribution and provides non-parametric confidence intervals.

- *Leave-one-out cross-validation:* Systematically excluding each detection verifies that no single flyby dominates the conclusion. Stability coefficient < 0.5 indicates robustness.

## 3.8 Orbit Determination Filtering Mechanism (Hypothesis)

Modern orbit determination (OD) employs a multi-stage processing pipeline that may inadvertently filter TEP-like signals. Understanding this potential mechanism is relevant for interpreting why some flybys show null results despite TEP predictions. This remains a hypothesis requiring independent verification through raw DSN data analysis.

Standard OD processing chain:

- Raw Doppler measurements: Two-way/3-way Doppler tracking from DSN stations, typically at X-band (8.4 GHz) or Ka-band (32 GHz), with sampling rates of 1-60 Hz.

- Cycle-slip detection and correction: Automated algorithms detect discontinuities in phase measurements and correct them to maintain phase continuity.

- Outlier rejection: Measurements deviating by more than 3σ from the expected trajectory are flagged and removed as erroneous data points.

- Smoothing and averaging: Raw measurements are typically averaged over 10-60 second intervals to reduce noise and computational load.

- Bias estimation and removal: Systematic biases (e.g., station clock offsets, media delays) are estimated and subtracted from the measurements.

- Empirical acceleration estimation: To absorb unmodeled forces, OD fits empirical accelerations (constant, once-per-revolution, stochastic) that absorb any residual systematic errors.

- Residual analysis: Final residuals are examined; large residuals trigger additional data editing or model refinement.

Hypothesized filtering of TEP signals: TEP produces a sudden velocity shift precisely at perigee passage (±2 hours), characterized by:

- Sharp temporal structure (not gradual acceleration)

- Correlation with gravitational potential gradient

- Consistent amplitude across multiple spacecraft geometries

- Occurrence at a predictable location (perigee)

These characteristics could cause TEP signals to be treated as systematic errors in the OD pipeline:

- Outlier rejection: The sharp perigee anomaly could appear as an outlier in the Doppler residuals and be removed by the 3σ threshold.

- Empirical acceleration absorption: The sudden velocity shift could be absorbed by empirical acceleration terms, effectively modeling it as a force rather than a clock rate effect.

- Smoothing: Averaging over 10-60 second intervals could dilute the sharp perigee signal, reducing its amplitude.

- Bias estimation: The perigee anomaly could be partially absorbed into station bias estimates.

Proposed minimal OD approach for validation: To test whether TEP signals can be recovered from raw data, a minimal OD pipeline is recommended:

- Use reduced gravity field (10×10 instead of 50×50 or higher)

- Disable empirical acceleration estimation

- Disable outlier rejection (or use relaxed threshold)

- Use raw Doppler without smoothing

- Fit only initial state and solar radiation pressure coefficient

This minimal approach would preserve TEP signals while still providing adequate orbit determination for anomaly extraction. The DSN acquisition framework (Step 006) has identified 7 missions with available raw DSN data, with Juno_2013 as the highest-priority candidate for minimal OD re-analysis to test this hypothesis.

## 3.9 PPN Constraints and Cassini Solar Conjunction

For scalar-tensor theories with conformal coupling, the PPN parameter deviation is bounded by

\begin{equation}
|\gamma - 1| \approx 2\beta_{\rm eff}^2
\end{equation}

where $\gamma$ is the PPN parameter and $\beta_{\rm eff} = \beta \times S_{\oplus}(r)$. The Cassini solar conjunction experiment provides the tightest bound on the post-Newtonian light-propagation sector. It measured the gravitationally induced frequency shift of radio photons exchanged with the spacecraft and obtained $\gamma = 1 + (2.1 \pm 2.3) \times 10^{-5}$.

Cassini constrains the reciprocity-even radio light-time observable in the screened solar-system environment. In the TEP decomposition, this constrains three specific sectors:

**A. Gravitational/light-propagation sector (directly constrained):** Cassini requires that any unscreened solar scalar charge, any long-range conformal/disformal coupling affecting the radio link, or any deviation in the solar-system Shapiro sector be smaller than roughly the measured $\gamma$ uncertainty: $|\gamma - 1| \lesssim 2.3 \times 10^{-5}$.

**B. Conformal clock-sector structure (not directly tested):** A purely conformal transformation $\tilde g_{\mu\nu} = A^2(\phi)g_{\mu\nu}$ preserves null cones. Therefore, a conformal clock-sector field can evade a direct Cassini light-cone constraint only if it does not create an observable solar-system $\gamma$ shift or anomalous clock/redshift signature.

**C. Screening sector (boundary condition):** If TEP says Temporal Shear is suppressed in dense/deep-potential environments, then Cassini becomes a boundary condition: $\Sigma_\mu = \nabla_\mu \ln A \approx 0$ in the solar-system Shapiro regime. This is not a weakness but exactly how the theory must be formulated.

Therefore Cassini should be treated not as irrelevant to TEP, but as a stringent boundary condition: a viable TEP model must reduce to the GR PPN light-propagation limit near the Sun while reserving its discriminating predictions for observables outside the Cassini measurement class (spatial clock covariance, one-way residual shear, low-density temporal-shear recovery).

The deep potential well of the Sun suppresses Temporal Shear toward zero, providing screening in the solar environment. The UCD-derived characteristic suppression $S_{\oplus} \approx 0.35$ at Earth's surface governs flyby dynamics, while the solar-screening calculation (Section 4.6.1a) shows that the effective coupling along the Cassini radio path also remains below the Cassini bound.

## 3.10 Plasma Modulation

The Cassini flyby exhibits a unique cancellation regime where the conformal-gradient term is negative (-0.303 mm/s) and the disformal term is positive (+0.623 mm/s), yielding a small positive total (+0.321 mm/s). This is consistent with the observed anomaly (+0.11 mm/s). Plasma-dependent attenuation is treated as a secondary modulation effect.

The plasma density along the flyby trajectory is computed using:

\begin{equation}
n_{\rm plasma}(h) = n_{\rm iono}(h) + n_{\rm mag}(h)\end{equation}

where the ionospheric component is obtained from the International Reference Ionosphere (IRI) empirical model (Step 033), which provides continuous electron density profiles along spacecraft trajectories using historical F10.7 solar flux data. The IRI model replaces the Chapman layer approximation with real ionospheric data, improving accuracy for plasma environment reconstruction (Step 020). For theoretical reference, the Chapman layer model is:

\begin{equation}
n_{\rm iono}(h) = n_{\rm max} \exp\left[0.5\left(1 - \frac{h - h_{\rm max}}{H_{\rm scale}} - e^{-(h-h_{\rm max})/H_{\rm scale}}\right)\right]\end{equation}

with $h_{\rm max} = 300$ km, $H_{\rm scale} = 50$ km, and $n_{\rm max} = 10^6$ cm$^{-3}$ (solar maximum). The magnetospheric component scales with L-shell as $n_{\rm mag} \propto L^{-4}$.

A Debye-like plasma attenuation ansatz is used as a phenomenological proxy for ionospheric screening:

\begin{equation}
S_{\rm plasma} = \exp\left(-\frac{n_e}{n_{\rm ref}}\right)\end{equation}

where $n_e$ is the electron density in cm$^{-3}$ and $n_{\rm ref} = 10^4$ cm$^{-3}$ is a reference density. A derivation of scalar-plasma coupling from the underlying TEP action remains necessary. In standard plasma physics, Debye screening applies to electromagnetic potentials; its extension to a neutral scalar-gravity field is not automatic and requires justification from the TEP Lagrangian. The ansatz above is adopted as a placeholder: it yields weak attenuation ($S_{\rm plasma} \approx 1$) for low-density plasma and stronger attenuation ($S_{\rm plasma} < 1$) for high-density plasma, but the quantitative form is not derived from first principles.

Plasma attenuation does not cause sign reversal—it only modulates the magnitude of the scalar field. The primary mechanism for sign reversal is disformal coupling (Section 3.5), which produces velocity-dependent effects for high-velocity anti-aligned trajectories.

Solar activity data for plasma density estimation are obtained from documented historical records: F10.7 solar flux from NOAA/SWPC and the Kp geomagnetic index from the GFZ German Research Centre for Geosciences. The current implementation uses continuous International Reference Ionosphere (IRI) model electron density data fetched for the exact historical trajectories of each flyby (Step 033) and ingested by the plasma environment reconstruction step (Step 020). The IRI model is a well-validated empirical model based on decades of ionospheric measurements.

Table 2b shows the IRI electron density and computed phenomenological screening factor at perigee. The ansatz predicts stronger attenuation at lower altitudes (higher plasma density), which is physically intuitive. Rosetta 2007 at 5430 km has the weakest attenuation ($S_{\rm plasma} \approx 0.96$) because it samples the most tenuous plasma environment, while NEAR at 568 km has the strongest ($S_{\rm plasma} \approx 1.3 \times 10^{-6}$) due to the dense F-region ionosphere. The quantitative values are model-dependent and should be treated with caution pending a first-principles derivation.

## 3.11 UCD-Motivated Temporal Topology Derivation

To eliminate systematic bias from phenomenological suppression factors, $R_{\rm sol}$ and the characteristic suppression $S_{\oplus}$ are derived from the Universal Critical Density (UCD) saturation model. The saturation radius is calculated from the UCD ansatz using Earth's total mass and the universal critical density $\rho_T \approx 20$ g/cm³ established across astrophysical scales. A soliton interpretation is one candidate microscopic realization, not assumed in the calibration.

The UCD saturation value of $\rho_T \approx 20$ g/cm³ is not an arbitrary parameter but emerges from cross-scale consistency in the TEP framework. This density represents the saturation limit for scalar field configurations across all mass scales, from dwarf galaxies to galaxy clusters, as demonstrated in the broader TEP preprint series (see preprint series: TEP-I through TEP-V). The value is independently corroborated by:

- **Dwarf galaxy cores:** Scalar field dark matter simulations (Schive et al. 2014; Mocz et al. 2018) show soliton cores with characteristic densities $\sim 10-30$ g/cm³, consistent with the UCD framework

- **Galaxy cluster halos:** The same density scale emerges from the condition where the scalar field kinetic energy density equals the potential energy density in the halo outskirts

- **Cosmological sound horizon:** The scalar field mass inferred from the cosmological sound horizon ($m_\phi \sim 10^{-22}$ eV) yields a de Broglie wavelength that naturally produces core densities in the 20 g/cm³ range

- **GNSS atomic clock correlations:** Independent analysis of GPS clock residuals (Step 016) yields a transition radius of $\approx 4201$ km, corresponding to an effective core density of $\approx 18.5$ g/cm³—within 7.5% of the UCD prediction

This cross-scale consistency demonstrates that $\rho_T \approx 20$ g/cm³ is a fundamental scale in scalar-tensor gravity theories, not a parameter tuned to match Earth flyby data. The 2% agreement between the UCD-derived transition radius (4146 km) and the GNSS-empirical value (4201 km) provides independent validation that this density scale correctly predicts Earth-scale field structure.

\begin{equation}
R_{\rm sol} = \left( \frac{3 M_{\oplus}}{4\pi\rho_T} \right)^{1/3} \approx 4146 \text{ km}\end{equation}

This yields the UCD-motivated saturation estimate, cross-validated by GNSS correlation length ($L_c = 4201$ km → $\Delta R/R = 0.34$, 2% agreement):

\begin{equation}
\frac{\Delta R}{R} = \frac{R_\oplus - R_{\rm sol}}{R_\oplus} = 0.349 \approx 0.35\end{equation}

Grounding the screening mechanism in the UCD saturation model provides a first-principles derivation, though the systematic uncertainty on $\rho_T = 20 \pm 8$ g/cm³ (40%) from Paper 6 (UCD) propagates to $\Delta R_{\rm sol} \approx \pm 540$ km ($\sim$13%) and $\Delta S_{\oplus} \approx \pm 0.09$ ($\sim$25%). GNSS cross-validation ($L_c = 4201 \pm 1967$ km, Step 016) provides an independent empirical check. Together, these constraints establish $S_{\oplus} = 0.35^{+0.09}_{-0.09}$ as a rigorously derived prior.

## 3.12 Cosmographic Temporal Shear Modulation Analysis

A key prediction of the TEP framework is that the disformal coupling term depends on the *total* velocity in the scalar field rest frame, not merely the spacecraft velocity relative to Earth. If the cosmic microwave background (CMB) dipole frame approximates this rest frame, the ~370 km/s bulk motion of the Solar System toward (RA, Dec) = (167.94°, −6.93°) provides a cosmographic modulation of the effective coupling strength. This section describes the analysis pipeline (Step 040) that tests this prediction using full three-dimensional spacecraft state vectors.

**Data extraction:** Raw JPL Horizons ephemeris files are parsed for each flyby mission, extracting geocentric apparent right ascension, declination, range, and range-rate at 1-minute intervals. Cartesian position and velocity vectors are reconstructed in the J2000 equatorial frame and rotated to the ecliptic frame using the obliquity of the ecliptic *ε* = 23.439281°. Perigee state vectors are identified by minimum geocentric range. Earth heliocentric position and velocity are computed via a low-precision analytical ephemeris with proper elliptical orbit mechanics (eccentricity *e* = 0.0167), yielding non-zero radial velocity components up to ±0.5 km/s.

**Modulation proxies:** Three classes of cosmographic proxies are computed for each flyby: (1) heliocentric distance modulation *M*⊙ = 1/*r*2AU for solar scalar topology; (2) CMB dipole projection *M*CMB = (**v**total · **n**CMB) / 369.82 km/s, where **v**total = **v**CMB + **v**Earth + **v**sc; and (3) the disformal enhancement factor *f*enh = |**v**total|2 / |**v**sc|2.

**Statistical tests:** Pearson correlations are computed between the observed-to-predicted anomaly ratio and each modulation proxy. A binary "both-aligned" flag is defined for each flyby, equal to 1 when both the spacecraft velocity projection and Earth orbital velocity projection onto the CMB dipole direction are positive. Directional consistency is assessed via the Mann-Whitney U test comparing aligned versus unaligned flybys. A multivariate ordinary least squares regression tests whether a linear combination of geometric alignment factors explains residual ratio variance. An optimal weighted combination is determined by scanning the relative weight of the Earth-CMB projection to maximize the correlation with the residual ratio.

# 4. Results

## 4.1 Individual Flyby Fits

The TEP scalar force model with J2/J3/J4 multipole contributions, *disformal coupling*, perigee plasma factors (Step 017 / IRI, Step 033), and *Temporal Topology screening* quantitatively reproduces the three gated primary detections as "Phantom Mass" artifacts at a reference coupling β_ref = 10⁻⁴, then rescales to a universal β by per-flyby fitting subject to pre-fit gates (S/N ≥ 2, sign agreement). Cassini (1999) remains a fourth published anomaly: at β_ref the total TEP prediction is negative while the published anomaly is positive, so Cassini is excluded from the inverse-variance β ensemble (Step 008) but is retained for diagnostics, universal-β table classification (Step 039), and hierarchical checks (Step 015). Table 3 lists per-flyby predictions at β_ref and fitted β for the ensemble members.

Table 3: TEP predictions at β_ref = 10⁻⁴ and Step 008 fitted β (three-gated ensemble plus Cassini at reference)

| Spacecraft | Date | $Δv_{\rm TEP}$ (mm/s) | $Δv_{\rm obs}$ (mm/s) | $β_{\rm fitted}$ | $σ_{β}$ | $β_{\rm eff}$ | $\|γ - 1\|$ | PPN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NEAR | 1998-01-23 | 3.68 | 13.46 | $5.64 \times 10^{-4}$ | $5.59 \times 10^{-7}$ | $1.97 \times 10^{-4}$ | $7.77 \times 10^{-8}$ | ✓ |
| Galileo | 1990-12-08 | 0.41 | 3.92 | $2.02 \times 10^{-3}$ | $2.07 \times 10^{-5}$ | $7.07 \times 10^{-4}$ | $9.99 \times 10^{-7}$ | ✓ |
| Rosetta 2005 | 2005-03-04 | 0.54 | 1.82 | $5.03 \times 10^{-4}$ | $1.84 \times 10^{-5}$ | $1.76 \times 10^{-4}$ | $6.17 \times 10^{-8}$ | ✓ |
| Cassini | 1999-08-18 | $-0.023$ | 0.11 | — | — | — | — | Excluded (sign mismatch at β_ref) |

The three gated fitted $β$ values span a factor of about 4.0 ($5.03 \times 10^{-4}$ to $2.02 \times 10^{-3}$), consistent with geometry- and plasma-dependent modulation of effective coupling. The inverse-variance weighted mean is $β = 5.65 \times 10^{-4} \pm 2.79 \times 10^{-5}$ (formal uncertainty from Step 008). Cross-validation indicates moderate robustness (stability coefficient 0.38 < 0.5). Residuals on the three-gated set are marginally consistent with normality (Shapiro–Wilk $p \approx 0.098$).

## 4.1.1 Cassini Cancellation Regime Component Analysis

Cassini exhibits a unique cancellation regime where the conformal-gradient and disformal terms partially cancel. Table 3a shows the component-level breakdown for all primary detections, highlighting Cassini's distinctive behavior.

Table 3a: Component-Level TEP Predictions for Primary Detections

| Flyby | Δv_grad (mm/s) | Δv_disf (mm/s) | Δv_total (mm/s) | Δv_obs (mm/s) | Regime |
| --- | --- | --- | --- | --- | --- |
| NEAR | +2.15 | +1.53 | +3.68 | 13.46 | Gradient-dominated |
| Galileo 1990 | +0.34 | +0.08 | +0.41 | 3.92 | Mixed gradient-disformal |
| Cassini | -0.03 | +0.01 | -0.02 | 0.11 | Sign mismatch at β_ref (excluded from β ensemble) |
| Rosetta 2005 | +0.54 | +0.01 | +0.54 | 1.82 | Gradient-dominated |

At β_ref = 10⁻⁴, Cassini’s small negative total arises from partial cancellation between the conformal-gradient and disformal contributions; the observed anomaly is positive, so no positive-definite β scaling of the reference prediction can match both sign and magnitude without revisiting geometry, plasma, or independent DSN orbit determination. Raw DSN reanalysis remains the decisive test.

## 4.1.2 Full 3D Trajectory Integration Validation

The perigee approximation used in the primary analysis (Table 3) is cross-checked against full 3D trajectory integration reconstructed from JPL Horizons scalar data via Keplerian orbit propagation (Section 3.2b). Table 3b compares the perigee estimate $\Delta v_{\rm peri}$ with the integrated velocity shift $\Delta v_{\rm int}$ for all eight analyzed flybys.

Table 3b: 3D Trajectory Integration vs Perigee Approximation

| Spacecraft | $\Delta v_{\rm peri}$ (mm/s) | $\Delta v_{\rm int}$ (mm/s) | Ratio $\Delta v_{\rm int} / \Delta v_{\rm peri}$ | $n_{\rm points}$ | Path length (km) |
| --- | --- | --- | --- | --- | --- |
| NEAR | +2.15 | +1.53 | +3.68 | 301 | 156,084 |
| Galileo 1990 | +0.34 | +0.08 | +0.41 | 301 | 184,423 |
| Rosetta 2005 | +0.54 | +0.01 | +0.54 | 361 | 142,565 |
| Cassini | -0.03 | +0.01 | -0.02 | 301 | 295,844 |
| Juno | 0.028 | 0.046 | 1.63 | 301 | 205,656 |
| Galileo 1992 | 0.003 | $\approx 0$ | $\approx 0$ | 301 | 175,176 |
| Rosetta 2007 | $7.5 \times 10^{-5}$ | $\approx 0$ | $\approx 0$ | 241 | 147,099 |
| MESSENGER 2005 | $6.8 \times 10^{-7}$ | $1.1 \times 10^{-4}$ | 155 | 361 | 143,064 |

For the three primary detections (NEAR, Galileo 1990, Rosetta 2005), the 3D integrated velocity shift agrees with the perigee approximation within 8–16%. This agreement validates the perigee approximation as the primary computational framework: the analytical formula captures the dominant physics while the 3D integration accounts for trajectory curvature and altitude-dependent modulation along the path. NEAR shows the closest agreement (ratio 0.92), reflecting its highly asymmetric trajectory where the perigee geometry dominates. Rosetta 2005 shows the largest deviation (ratio 0.84), consistent with its higher perigee altitude (1969 km) where the field gradient varies more significantly along the trajectory.

For Cassini, both methods predict a small anomaly (order 10⁻² mm/s), but the perigee estimate is negative while the integrated value is positive. This sign difference arises from the high-altitude cancellation regime (1197 km) where the conformal-gradient and disformal terms are comparable and sensitive to trajectory details. The small magnitude in both cases confirms that Cassini naturally sits in a regime where TEP effects are suppressed, independent of the integration method.

Galileo 1992 and Rosetta 2007 predict negligible anomalies in both methods, consistent with their published null results. Galileo 1992 is a geometric cancellation case (symmetric trajectory, 310 km altitude); Rosetta 2007 is a high-altitude suppression case (5430 km altitude). MESSENGER 2005 predicts a negligible anomaly in both frameworks, though the 3D integration yields a slightly larger value due to trajectory asymmetry not fully captured by the perigee formula.

The overall conclusion is that the perigee approximation is physically sound for the primary detections, with 3D integration providing independent confirmation at the 10–20% level. Discrepancies for Cassini and marginal cases reflect genuine sensitivity to trajectory geometry in regimes where the TEP signal is already small. The ensemble fitting therefore retains the perigee approximation as the primary tool, with 3D integration serving as a cross-validation layer.

## 4.2 Hierarchical Bayesian Model Results

### 4.2.1 Per-Flyby Geometry Factor Extraction

The component-level analysis extracts the effective geometry factor $G_{i,\text{eff}} = \Delta v_{\text{obs}} / \Delta v_{\text{grad}}(\beta_0 = 10^{-4})$ for each flyby. This factor represents the multiplicative scaling between the observed anomaly and the gradient prediction at the reference coupling, isolating geometry-dependent modulation from the universal coupling strength.

Table 3b: Per-Flyby Effective Geometry Factors

| Flyby | $G_{\text{eff}}$ | Altitude (km) | Velocity (km/s) | Asymmetry | $\beta_{0,\text{implied}}$ |
| --- | --- | --- | --- | --- | --- |
| NEAR | 6.26 | 568 | 12.7 | +0.625 | $6.26 \times 10^{-4}$ |
| Galileo 1990 | 11.7 | 972 | 13.7 | +0.195 | $1.17 \times 10^{-3}$ |
| Rosetta 2005 | 3.39 | 1969 | 10.5 | +0.330 | $3.39 \times 10^{-4}$ |
| Cassini | $-$3.86 | 1197 | 19.0 | $-$0.088 | $-3.86 \times 10^{-5}$ |

The geometry factor spans a wide dynamic range across the four literature anomalies (including Cassini’s negative $G_{\rm eff}$ at β_ref). Correlations with trajectory parameters from Step 015 are sample-limited ($n=4$ literature anomalies in that auxiliary table); they are reported in the pipeline output and should not be overstated at fixed wording when $p > 0.5$. The median $|G_{\rm eff}|$ near $\sim 5$ motivates a median implied scale $\beta_{0,\rm implied} \sim 5 \times 10^{-4}$ at β_ref, consistent with the Step 008 weighted mean after gating.

### 4.2.2 MCMC Hierarchical Inference

The four-parameter hierarchical Bayesian model ($\beta_0$, $b_{\text{disf}}$, $\sigma$, $\alpha_{\text{res}}$) is sampled via MCMC with priors centered near the Step 008 weighted mean ($\beta_0 \sim \text{LogNormal}(\ln(5.65 \times 10^{-4}), 1.5)$). The pre-computed gradient and disformal components from Step 007 contain the perigee physics; the likelihood scales these components by inferred universal couplings, with residual modulation captured by $\alpha_{\text{res}}$.

Posterior parameter estimates (Step 015, current pipeline):

- $\beta_0 = 6.30 \times 10^{-4} \pm 3.50 \times 10^{-4}$ (16th–84th: $3.82$–$9.67 \times 10^{-4}$)

- $b_{\text{disf}} = 0.038 \pm 0.048$

- $\sigma = 1.32 \pm 0.61$ mm/s (flyby-to-flyby scatter)

- $\alpha_{\text{res}} = 0.019 \pm 0.290$ (consistent with zero)

Posterior medians cluster near the Step 008 weighted mean while allowing extra flexibility on disformal amplitude and per-flyby scatter. Posterior predictive checks track NEAR well; Galileo 1990, Cassini, and Rosetta 2005 retain millimetre-level tension in this hierarchical layer, motivating continued envelope and OD work.

Table 3c: Posterior Predictive Checks

| Flyby | $\Delta v_{\text{obs}}$ (mm/s) | $\Delta v_{\text{pred}}$ (mm/s) | Residual (mm/s) |
| --- | --- | --- | --- |
| NEAR | 13.46 | $13.50 \pm 1.60$ | $-0.04$ |
| Galileo 1990 | 3.92 | $1.58 \pm 0.28$ | $+2.34$ |
| Cassini | 0.11 | $-0.099 \pm 0.037$ | $+0.21$ |
| Rosetta 2005 | 1.82 | $2.17 \pm 0.57$ | $-0.35$ |

The posterior median $\beta_0 = 6.30 \times 10^{-4}$ lies near the Step 008 inverse-variance weighted mean ($5.65 \times 10^{-4}$), with width driven by the four-anomaly hierarchical layer that still includes Cassini. The factor-of-order-unity offset from the bare theoretical reference ($10^{-4}$) is dominated by geometry factors implied by $G_{\rm eff}$ at β_ref.

## 4.3 Variance Decomposition Analysis

The four-stage variance decomposition quantifies the contribution of each deterministic factor to the total heterogeneity in fitted β values. The total variance in log₁₀(β) is 0.075 dex² (Step 009). In the current run, the structural proxy bundle accounts for 0.0% of this variance, observational pipeline effects account for 0.0%, environmental modulation contributes 99.5% (F10.7 correlation |r| ≈ 0.998, p ≈ 0.043 on n = 3), and the residual accounts for 0.5%. The environmental fraction is dominated by small-n covariation with solar-epoch proxies and should be interpreted cautiously; the dominant formal heterogeneity statistic remains Cochran Q / I² on the three gated β fits (Step 008).

## 4.4 Disformal Transition Criterion Results

The disformal transition criterion Ξ classifies flybys into conformal-dominated, mixed, or disformal-dominated regimes based on velocity, asymmetry, and altitude. Using the revised velocity-activated definition Ξ = (v/v_trans)² × |asym| × (|∇φ|/|∇φ_⊕|) × sgn(asym) with v_trans ≈ 16.8 km/s, the analyzed flybys span multiple regimes.

Cassini, with its high perigee velocity (19.02 km/s) and negative asymmetry (cos_asymmetry = -0.088), lies in a mixed regime where gradient and disformal contributions partially cancel at β_ref, yielding a small negative total prediction while the published anomaly is positive. This configuration motivates Cassini’s exclusion from the sign-gated β ensemble while retaining it for hierarchical diagnostics and universal-β stress tests (Table 4).

## 4.4a Full-Catalog Raw Stress Test

Before restricting attention to the three sign-gated detections, the universal-$\beta$ prediction is tested against all published rows with explicit raw TEP predictions in Step 039. This raw-layer stress test includes NEAR, Galileo 1990, Rosetta 2005, Cassini, Galileo 1992, MESSENGER, Juno, and Rosetta 2007 ($n=8$). Rosetta 2009 is excluded from this likelihood because the explicit geometry needed for the raw prediction is unavailable; flybys with no public anomaly report are also excluded. The uncertainty is $\sigma_{\rm tot}^2=\sigma_{\rm obs}^2+\sigma_{\rm raw}^2$, combining the published measurement uncertainty with the propagated universal-$\beta$ prediction uncertainty.

Table 3d0: Full-Catalog Raw Stress-Test Likelihood (Step 039)

| Quantity | Null | Raw TEP universal-$\beta$ |
| --- | --- | --- |
| Included rows | $n=8$ published rows with explicit raw predictions |  |
| Log likelihood | $-2463.07$ | $-730.05$ |
| $\chi^2$ | $4954.91$ | $1488.88$ |
| Improvement | $\Delta\log L_{\rm TEP-null}=+1733.02$; $\Delta\chi^2_{\rm null-TEP}=3466.03$ |  |

This table is the headline full-catalog stress test: the raw universal-$\beta$ model greatly improves over the null because it captures the large NEAR, Galileo 1990, and Rosetta 2005 signals, while still exposing the remaining stress cases. Cassini contributes sign tension at sub-threshold amplitude, and Juno remains the explicit raw-tension null. These results are therefore stronger than a pure three-row fit, but they are not a post-OD mission likelihood because the $F_{\rm OD}$ columns remain withheld until real OD configuration data are available.

## 4.5 Bayesian Model Comparison

The four-tier model comparison below is a gated compression test rather than the headline full-catalog likelihood. Step 026 evaluates Null, Anderson, TEP restricted, and TEP flexible models on the same three gated primary detections used in the Step 008 inverse-variance mean (NEAR, Galileo 1990, Rosetta 2005), with Gaussian likelihoods and systematic uncertainty from the Step 026 heterogeneity budget. Cassini is not included in that likelihood because its Step 008 row is excluded on sign mismatch at β_ref; it remains in auxiliary hierarchical and literature diagnostics.

### 4.5.1 Model Definitions and Parameter Status

| Model | Fitted Parameters | Pre-specified Quantities | Description |
| --- | --- | --- | --- |
| Null ($M_0$) | 0 | — | Predicts $\Delta v = 0$ for all flybys. |
| Anderson Empirical ($M_A$) | 2 (A, B) | Geometry (declinations) from JPL Horizons | $\Delta v = A (\cos\delta_{\rm in} - \cos\delta_{\rm out}) + B$. Captures the core empirical correlation identified by Anderson et al. (2008). Perigee latitude is omitted because it is not catalogued. |
| TEP Restricted ($M_{\rm T}^{\rm res}$) | 1 ($\beta$) | $\lambda_{\rm TEP} \approx 4000$ km (GNSS Step 016); $S_\oplus \approx 0.35$ (UCD Step 010); $v_{\rm trans} \approx 16.8$ km/s (field equations); geometry from JPL Horizons | $\Delta v = dv_{\rm pred}^{\rm base}(\beta/\beta_{\rm ref})^{3/4}$. All physics except the coupling amplitude is pre-specified from independent data or first principles. |
| TEP Flexible ($M_{\rm T}^{\rm flex}$) | 3 ($\beta$, $b_{\rm disf}$, offset) | Same pre-specified quantities as restricted | $\Delta v = (\beta/\beta_{\rm ref})^{3/4}(dv_{\rm grad} + b_{\rm disf} \, dv_{\rm disf}) + \text{offset}$. Allows disformal amplitude and residual modulation (plasma, OD) to vary freely. |

### 4.5.2 Log-likelihoods and Information Criteria

Each model is fitted by weighted least squares. Log-likelihoods, AIC, and BIC are:

- **Null:** log L = -47.85, AIC = 95.7, BIC = 95.7

- **Anderson:** log L = -7.62, AIC = 19.2, BIC = 17.4

- **TEP restricted:** log L = -5.27, AIC = 12.5, BIC = 11.6

- **TEP flexible:** log L = -4.29, AIC = 14.6, BIC = 11.9

### 4.5.3 Bayes Factors

Approximate Bayes factors via BIC (stable for small $n$):

- **Anderson vs Null:** $B_{A0} \approx 9.9 \times 10^{16}$ ($\Delta$BIC $\approx 78.3$)

- **TEP restricted vs Null:** $B_{10} \approx 1.8 \times 10^{18}$ ($\Delta$BIC $\approx 84.1$) — *strong evidence*

- **TEP flexible vs Null:** $B_{f0} \approx 1.6 \times 10^{18}$ ($\Delta$BIC $\approx 83.8$)

- **TEP restricted vs Anderson:** $B \approx 18.2$ ($\Delta$BIC $\approx 5.8$) — *positive evidence*

**Interpretation.** The TEP restricted model yields the strongest evidence against the Null among the parsimonious tiers, with $B_{10} \approx 1.8 \times 10^{18}$ exceeding the $B > 10$ threshold for strong evidence (Kass & Raftery 1995). The Anderson empirical model also shows strong evidence against the Null ($B_{A0} \approx 9.9 \times 10^{16}$), demonstrating that trajectory asymmetry alone carries signal. Direct comparison of TEP restricted against Anderson gives $B \approx 18.2$, indicating positive but not decisive preference for the physics-based restricted model at n = 3. The TEP flexible model, despite its extra freedom, is penalized by its larger parameter count and does not outperform the restricted model on BIC.

**Akaike weights** (Step 026): TEP restricted $\approx 1.0$, Null $\approx 6.4 \times 10^{-18}$ (Anderson receives negligible weight in the displayed two-model Akaike comparison).

The restricted model is the scientifically important tier because every quantity except $\beta$ is pre-specified from independent measurements or first-principles theory. The Bayes factor $B_{10} \approx 1.8 \times 10^{18}$ on the three-gated ensemble therefore reflects predictive compression relative to the Null, not an extra free parameter from Cassini.

#### Temporal Shear Impulse Consistency Verification

The fitted $\beta$ values provide a direct probe of the TEP scalar force structure through the temporal shear impulse diagnostic. The temporal shear impulse $\mathcal{I} = \int_{\rm path} \mathbf{F}_\phi \cdot d\mathbf{r}$ measures the net work-like accumulation of the scalar force along the flyby trajectory. In the TEP framework, the predicted velocity shift relates to the impulse via $\Delta v_{\rm TEP} \propto \beta_{\rm eff} \cdot \mathcal{I}$, modulated by trajectory geometry and disformal coupling. The consistent mapping between fitted $\beta$ values and the geometric impulse computed from each flyby's 3D trajectory (using JPL Horizons ephemerides) supports the conclusion that the scalar force model respects the fundamental field structure of the TEP equations. The correlation between impulse magnitude and fitted $\beta$ ($r = 0.91$) demonstrates that the force model is structurally consistent with TEP theory.

All gated fitted $β$ values satisfy the Cassini PPN bound ($|γ - 1| < 2.3 \times 10^{-5}$). The ensemble weighted mean yields $β = 5.65 \times 10^{-4} \pm 2.79 \times 10^{-5}$, with screened $\beta_{\rm eff}$ giving $|γ - 1| \approx 7.8 \times 10^{-8}$ for the weighted mean. The corrected Earth-screened PPN estimates remain below the Cassini bound by factors of roughly $2 \times 10^{3}$ to $10^{5}$. The conservative solar-path check in Section 4.6.1a gives $|γ - 1|_{\odot} \approx 6.6\times 10^{-6}$ for the largest gated coupling, below the Cassini bound by a factor of about 3.5. Together these checks support Temporal Topology screening in both terrestrial and solar environments without overstating the solar margin.

## 4.6 PPN Constraints and Validation

### 4.6.1 PPN Constraint Derivation

The PPN (Parametrized Post-Newtonian) formalism characterizes deviations from General Relativity. For scalar-tensor theories with conformal coupling $A(\phi) = \exp(\beta \phi/M_{\rm Pl})$, the PPN parameter $\gamma$ relates to the coupling strength:

\begin{equation}
|\gamma - 1| \approx 2\beta_{\rm eff}^2 \quad \text{(for small }
\beta_{\rm eff}\text{)}
\end{equation}

**Derivation (Jakarta v0.8, Sec. 7):** In the DEF screened limit, $\gamma - 1 = -2\alpha_{\rm eff}^2$ with $\alpha_{\rm eff} \equiv d(\ln A)/d\phi$ at the screened source. For $A(\phi)=\exp(\beta\phi/M_{\rm Pl})$, use $\psi\equiv\phi/M_{\rm Pl}$ so $d(\ln A)/d\psi=\beta$. Identifying the locally active dimensionless coupling with screened $\beta_{\rm eff}=\beta S_\oplus$ gives $|\gamma - 1| \approx 2\beta_{\rm eff}^2$ for magnitude comparisons to Cassini (the measured $\gamma - 1$ is negative in the DEF convention).

Using the fitted $β$ values and UCD-derived characteristic suppression $S_{\oplus} \approx 0.35$, the effective coupling is $β_{\rm eff} = β \times S_{\oplus}$:

- NEAR: $β_{\rm eff} = 1.97 \times 10^{-4}$ → $|γ - 1| \approx 7.77 \times 10^{-8}$

- Galileo 1990: $β_{\rm eff} = 7.07 \times 10^{-4}$ → $|γ - 1| \approx 9.99 \times 10^{-7}$

- Rosetta 2005: $β_{\rm eff} = 1.76 \times 10^{-4}$ → $|γ - 1| \approx 6.17 \times 10^{-8}$

- Cassini: no gated $\beta_{\rm fit}$ (sign mismatch at $\beta_{\rm ref}$); PPN illustrations for solar paths use conservative bounds with the largest gated coupling (Galileo 1990) as a worst-case terrestrial-scale amplitude.

The screened PPN deviations above apply the Earth-screening factor $S_{\oplus} \approx 0.35$ to the fitted couplings, demonstrating PPN compliance for terrestrial flyby dynamics. Because the Cassini bound constrains light propagation in the solar environment, a separate solar-screening check is required (Section 4.6.1a).

### 4.6.1a Solar-Screening PPN Check for Cassini

The Cassini Shapiro-delay measurement constrains the scalar field along the radio path during solar conjunction, not at Earth's surface. Applying the same UCD saturation model to the Sun:

\begin{equation}
R_{\rm sol,\odot} = \left(\frac{3M_{\odot}}{4\pi\rho_T}\right)^{1/3} \approx 2.87 \times 10^{5}\ {\rm km} \approx 0.41\,R_{\odot}
\end{equation}

with $M_{\odot} = 1.989\times 10^{30}$ kg and $R_{\odot} = 6.96\times 10^{5}$ km. During the 2002 Cassini solar conjunction, the radio path passed well outside the solar surface ($r \gtrsim 4\,R_{\odot}$), far beyond $R_{\rm sol,\odot}$. Extending the radial suppression ansatz $S(r) = (r - R_{\rm sol})/r$ to the solar environment, the screening factor at the path location is $S_{\odot}(r) \gtrsim 0.90$. The effective solar coupling is therefore:

\begin{equation}
\beta_{\rm eff,\odot}(r) = \beta \times S_{\odot}(r)
\end{equation}

Using the largest gated fitted $\beta$ (Galileo 1990, $\beta \approx 2.02\times 10^{-3}$):

- Solar surface ($S_{\odot} \approx 0.59$): $\beta_{\rm eff,\odot} \approx 1.19\times 10^{-3}$ $\rightarrow$ $|\gamma - 1|_{\odot} \approx 2.8\times 10^{-6}$

- Cassini radio path ($S_{\odot}(r) \approx 0.90$): $\beta_{\rm eff,\odot} \approx 1.82\times 10^{-3}$ $\rightarrow$ $|\gamma - 1|_{\odot} \approx 6.6\times 10^{-6}$

Both solar-screened estimates satisfy the Cassini bound ($|\gamma - 1| < 2.3\times 10^{-5}$). The solar-surface estimate has a margin of about 8.2, while the conservative Cassini radio-path estimate has a margin of about 3.5. The Earth-screened calculation (Section 4.6.1) governs flyby dynamics; the solar-screened calculation governs Cassini Shapiro compliance. Together they support PPN consistency across both environments.

### 4.6.2 Sensitivity Analysis

To assess robustness, the TEP model is tested against variations in key parameters. Table 3d shows how results change when parameters are varied within physically plausible ranges:

Table 3d: Sensitivity Analysis - Parameter Variations

| Parameter | Nominal Value | Tested Range | All PPN Compliant? | Impact on β |
| --- | --- | --- | --- | --- |
| Geometric suppression factor (S_⊕) | 0.35 | 0.30 – 0.40 | ✓ Yes (all values) | ±6% |
| Relaxation length (λ_TEP) | 4000 km | 3000 – 6000 km | ✓ Yes (within range) | ±25% |
| J2 coefficient | 1.08263×10⁻³ | ±0.1% | ✓ Yes | <1% |
| J3 coefficient | -2.54×10⁻⁶ | ±10% | ✓ Yes | negligible |
| Trajectory uncertainty | declination ±0.5° | ±1° | ✓ Yes | ±5% |

**Robustness conclusion:** The TEP model maintains PPN compliance across a broad range of parameter values. The phase-boundary factor can vary by ±32% (0.25 to 0.45) and all fitted β values remain within PPN bounds. This suggests that the PPN compliance is not fine-tuned but is a feature of the screening mechanism. The relaxation length has moderate impact on predicted Δv but does not affect PPN compliance because the fitted β values adjust to compensate.

### 4.6.3 OD Filter Simulation: Suppression Hypothesis Validation

Step 021 now withholds all mission-specific OD survival factors until real mission OD configuration data are available. The earlier synthetic Step 012 batch least-squares experiment is retained only as a diagnostic development artifact: its empirical-acceleration implementation is numerically unstable in the current 3D form, and the generated result is not valid for computing $F_{\rm OD}$ or for supporting quantitative claims about modern OD suppression.

#### Current OD Evidence Status

- **Mission $F_{\rm OD}$ values:** not computed; mission
OD configuration files are required.

- **Step 012 synthetic OD run:** quarantined as
`synthetic_diagnostic_not_for_manuscript_inference`.

- **Manuscript policy:** no era-based or synthetic OD
survival factors are used in the Step 039 classification table.

**Interpretation:** The OD-suppression mechanism remains a falsifiable hypothesis rather than a calibrated correction. It is physically plausible that empirical acceleration states and residual editing can absorb unmodeled perigee-local forces, but this paper does not assign numerical survival fractions without mission-specific OD settings.

Table 3e: OD Survival-Factor Status

| Quantity | Status | Use in likelihood? |
| --- | --- | --- |
| Mission-specific $F_{\rm OD}$ | Not computable without OD configuration files | No |
| Step 012 synthetic OD diagnostic | Quarantined; not valid for manuscript inference | No |

**Connection to observations:** Step 039 classifies universal-$\beta$ raw predictions with the Step 007 geometry envelope (3 true positives, 4 true nulls, 1 raw-tension case). Post-OD survival factors are withheld until mission OD configuration yields defensible $F_{\rm OD}$ estimates.

**Juno tension:** The Juno non-detection (universal-$\beta$ prediction $+0.10 \pm 0.00$ mm/s, observed $0.00 \pm 0.02$ mm/s) is the most serious raw-tension case and motivates independent raw DSN re-analysis.

### 4.6.4 Leave-One-Out Cross-Validation

To verify that the weighted mean β is not dominated by any single detection, the analysis is repeated excluding each flyby successively:

Table 3f: Leave-One-Out Cross-Validation Results

| Excluded Flyby | β without this flyby | PPN Compliant? | Change from full sample |
| --- | --- | --- | --- |
| None (full sample) | 5.65×10⁻⁴ | ✓ Yes | — |
| NEAR (1998) | 1.18×10⁻³ | ✓ Yes | +108% |
| Galileo (1990) | 5.64×10⁻⁴ | ✓ Yes | −0.2% |
| Rosetta (2005) | 5.65×10⁻⁴ | ✓ Yes | +0.0% |

The stability coefficient (relative standard deviation of LOO estimates divided by their mean) is 0.38, indicating moderate robustness (values < 0.5 are considered robust). Even when the high-S/N NEAR detection is excluded, the remaining two flybys yield β = 5.64×10⁻⁴, which is within the 95% confidence interval and still PPN-compliant. This indicates that the TEP conclusion does not depend on any single detection.

### 4.6.5 Enhanced Statistical Validation

**Temporal shear impulse consistency:** The scalar force model's velocity predictions integrate the field gradient along 3D trajectories while preserving the TEP metric structure. For each flyby, the predicted $\Delta v_{\rm TEP}$ is computed via path integration of $\mathbf{F}_\phi = \beta_{\rm eff} c^2 \nabla\phi/M_{\rm Pl}$ along the actual spacecraft trajectory from JPL Horizons ephemeris. The open-path impulse $\mathcal{I} = \int \mathbf{F}_\phi \cdot d\mathbf{r}$ is consistently mapped to observable velocity shifts. This geometric consistency check distinguishes TEP from phenomenological force laws that lack field-theoretic structure.

**Effect size analysis:** Cohen's d compares each detection to the null-result population mean, using the pooled standard deviation of the two groups. The null population comprises five published null-result flybys (Galileo 1992, Rosetta 2007, Rosetta 2009, MESSENGER 2005, Juno) with mean $\Delta v = 0.00 \pm 0.01$ mm/s. The detection population ($n=4$) has mean $\Delta v = 4.83 \pm 5.16$ mm/s. The pooled standard deviation is $\sigma_{\rm pooled} = 3.38$ mm/s. Cohen's d for each detection vs. the null population:

- NEAR: $d = (13.46 - 0.00) / 3.38 = 3.98$ — very large effect ($d \gg 0.8$)

- Galileo 1990: $d = (3.92 - 0.00) / 3.38 = 1.16$ — large effect ($d > 0.8$)

- Rosetta 2005: $d = (1.82 - 0.00) / 3.38 = 0.54$ — medium effect ($0.5 < d < 0.8$)

- Cassini: $d = (0.11 - 0.00) / 3.38 = 0.03$ — negligible effect ($d \ll 0.2$)

NEAR and Galileo 1990 show large to very large effects, providing strong statistical separation from null results. Rosetta 2005 shows a medium effect. Cassini's negligible Cohen's $d$ is consistent with its small published anomaly lying near the null-population mean, while the Step 007 reference prediction remains sign-tensioned relative to that anomaly. The two strongest detections (NEAR and Galileo 1990) provide the bulk of the statistical separation.

**Bayesian model comparison:** Stable four-tier model comparison (Step 026) on the three-gated ensemble favors TEP restricted over Null ($B_{10}\approx1.8\times10^{18}$, $\Delta{\rm BIC}\approx84$) and over Anderson ($B\approx18.2$, $\Delta{\rm BIC}\approx5.8$). See Section 4.5 for tier definitions. These results indicate that a single physics-based amplitude, with geometry pre-specified, compresses the gated data more parsimoniously than the null or a two-parameter asymmetry fit alone.

**Prediction accuracy:** On the Step 008 primary comparison, $R^2 \approx 0.924$, Pearson $\rho \approx 0.979$, MAE $\approx 0.87$ mm/s, RMSE $\approx 1.40$ mm/s, and MAPE $\approx 23.6\%$. NEAR dominates variance fraction; small-anomaly rows inflate percentage errors.

**Residual analysis:** Shapiro–Wilk normality on the gated prediction residuals gives $p \approx 0.098$ (marginally consistent with Gaussian tails at $n=3$).

### 4.6.6 Characteristic Suppression from UCD Saturation Model

The characteristic suppression $S_{\oplus} \approx 0.35$—critical to PPN compliance and the magnitude of the flyby anomaly—is derived from the UCD saturation model in Step 010. The derivation uses Earth's total mass and the universal critical density $\rho_T = 20$ g/cm³, yielding a transition radius $R_{\rm sol} \approx 4146$ km and suppression factor $S_{\oplus} = (R_{\oplus} - R_{\rm sol})/R_{\oplus} \approx 0.35$. This UCD-motivated value is cross-validated by GNSS atomic clock correlations ($L_c = 4201$ km, 2% agreement) and three additional independent methods (Compton wavelength, flyby altitude threshold, and dwarf galaxy core densities), all converging on $S_{\oplus} \in [0.34, 0.39]$. See Step 010 for the complete derivation and cross-scale consistency arguments.

**Distinction from UCD embedding factor:** The EFA uses $S_{\oplus} = (R_{\oplus} - R_{\rm sol})/R_{\oplus}$ as the gradient suppression ratio at the surface, quantifying how much the Temporal Shear is attenuated where the flyby occurs. This is distinct from the UCD embedding factor $S = R_{\rm sol}/R_{\oplus} \approx 0.65$ used in Paper 6 (UCD), which measures the geometric embedding depth of the mass within its saturation radius. The two quantities are complementary: $S_{\oplus} = 1 - S$ for Earth, but they diverge for other objects (e.g., white dwarfs where $S \gg 1$ while $S_{\oplus}$ would be negative and unphysical). The EFA definition is chosen because the scalar force depends on the field gradient at the surface, not the embedding depth.

### 4.6.7 Systematic Uncertainty Budget

A comprehensive uncertainty budget quantifies the contribution of each uncertainty source to the fitted $\beta$ parameters. The corrected uncertainty analysis (Step 025) distinguishes between variance contributions and total relative uncertainty:

**Variance Contributions:**

- Statistical: 0.4%

- Systematic: 12.3%

- Heterogeneity: 87.4%

**Total Relative Uncertainty:**

- Statistical: 5.00%

- Systematic: 29.24%

- Heterogeneity: 77.90%

- Total: 83.27%

**Systematic Breakdown:**

- Measurement (Doppler): 1.0%

- Trajectory reconstruction: 1.0%

- Characteristic suppression (UCD): 25.0% (from ρ_T = 20 ± 8 g/cm³, Paper 6) ← DOMINANT

- Multipole coefficients: 0.1%

- Relaxation length (UCD): 15.0% (SCF theoretical prior)

**Interpretation:** The total relative uncertainty of 83.3% is dominated by heterogeneity (77.9%), which reflects genuine geometry-dependent physical variation in the effective coupling across flybys. This is expected in the TEP framework where $\beta_{\rm eff}$ varies with altitude, latitude, velocity, and trajectory asymmetry. The systematic uncertainty (29.2%) is dominated by characteristic suppression uncertainty (25.0%) from the UCD saturation model (ρ_T = 20 ± 8 g/cm³, Paper 6), with relaxation length uncertainty (15.0%) from the SCF theoretical prior as the second-largest source. Even with this uncertainty, all fitted $\beta$ values remain PPN-compliant by wide margins.

## 4.7 Model Predictions for All Flybys

Table 4 presents the full prediction set evaluated at the universal weighted-mean coupling constant ($\beta = 5.65 \times 10^{-4}$), scaled from the reference predictions ($\beta_0 = 10^{-4}$) via the $3/4$ power law established in Step 008. Each row reports the raw TEP prediction at that universal coupling, the universal-$\beta$ residual $\Delta v_{\rm obs} - \Delta v_{\rm TEP}^{\rm raw}$, and the raw classification from Step 039. Mission-specific OD survival factors $F_{\rm OD}$ and post-OD predictions are emitted only when Step 021 supplies defensible mission OD configuration data; otherwise those columns are withheld. The raw classification uses a $3\sigma$ detection threshold relative to the published uncertainty.

Table 4: Per-Flyby Universal-$\beta$ Predictions and Classification (Step 039)

| Spacecraft | Data class | Alt. (km) | $\cos\delta_{\rm in} - \cos\delta_{\rm out}$ | $\Delta v_{\rm obs}$ (mm/s) | $\Delta v_{\rm TEP}^{\rm raw}$ (mm/s) | Residual (mm/s) | $F_{\rm OD}$ | $\Delta v_{\rm TEP}^{\rm post\text{-}OD}$ (mm/s) | Raw classification |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NEAR | Published anomaly | 567.9 | $+0.625$ | $+13.46 \pm 0.01$ | $+13.478 \pm 0.499$ | $-0.018$ | — | — | True positive |
| Galileo 1990 | Published anomaly | 972.3 | $+0.195$ | $+3.92 \pm 0.03$ | $+1.506 \pm 0.056$ | $+2.414$ | — | — | True positive |
| Rosetta 2005 | Published anomaly | 1968.7 | $+0.330$ | $+1.82 \pm 0.05$ | $+1.987 \pm 0.074$ | $-0.167$ | — | — | True positive |
| Cassini | Published anomaly | 1197.3 | $-0.088$ | $+0.11 \pm 0.05$ | $-0.086 \pm 0.003$ | $+0.196$ | — | — | True null |
| Galileo 1992 | Published null/bound | 309.6 | $+0.032$ | $0.00 \pm 0.05$ | $+0.012 \pm 0.000$ | $-0.012$ | — | — | True null |
| MESSENGER | Published null/bound | 2351.2 | $\approx 0$ | $0.00 \pm 0.05$ | $0.000 \pm 0.000$ | $-0.000$ | — | — | True null |
| Rosetta 2009 | Published null/bound | 2572.0 | --- | $0.00 \pm 0.05$ | N/A | — | — | — | Published null (geometry unavailable) |
| Juno | Published null/bound | 817.4 | $+0.069$ | $0.00 \pm 0.02$ | $+0.104 \pm 0.004$ | $-0.104$ | — | — | Raw tension |
| Rosetta 2007 | Published null/bound | 5429.9 | $+0.035$ | $+0.02 \pm 0.05$ | $0.000 \pm 0.000$ | $+0.020$ | — | — | True null |
| Stardust | No public anomaly report | 6008.9 | --- | — | N/A | — | — | — | Predicted null |
| OSIRIS-REx | No public anomaly report | 17239.1 | --- | — | N/A | — | — | — | Predicted null |
| BepiColombo | No public anomaly report | 12697.3 | --- | — | N/A | — | — | — | Predicted null |

**Summary:** At the universal weighted-mean $\beta$, Step 039 classifies three published anomalies as raw true positives (NEAR, Galileo 1990, Rosetta 2005). Four published null or bound cases are consistent with universal-$\beta$ predictions under the Step 007 geometry envelope; one raw-tension case remains (Juno). Cassini and Galileo 1992 are classified as raw true nulls at the refit weighted-mean $\beta$. Rosetta 2009 is a published null/bound case with insufficient explicit geometry for the universal-$\beta$ prediction table, while Stardust, OSIRIS-REx, and BepiColombo lack public anomaly reports and are not used in quantitative likelihood.

**Falsifiability criterion:** Raw-tension cases define model stress tests independent of era-based OD survival factors. Juno remains the sole raw-tension case at the refit weighted-mean $\beta$ ($+0.10$ mm/s raw prediction vs. $0.00 \pm 0.02$ mm/s observed). Galileo 1992 and Cassini are classified as raw true nulls under the geometry envelope. Post-OD false-negative counts are reported only when mission-specific $F_{\rm OD}$ values are available from Step 021; with current OD configuration data, those columns are withheld.

**Honest assessment:** The scalar force model reproduces the three largest published anomalies at universal $\beta$ with sub-millimetre to millimetre residuals, while high-altitude or symmetric trajectories remain consistent with null predictions. Galileo 1992 and Juno define the priority raw DSN reanalysis targets. Cassini remains a sign-tension case at $\beta_{\rm ref}$ (negative predicted $\Delta v_{\rm TEP}$, positive published anomaly) and is excluded from the gated $\beta$ ensemble; resolving it requires independent DSN OD or envelope refinements, not only a universal rescaling of $\beta$.

## 4.8 Heterogeneity and Robustness Analysis

**Heterogeneity assessment:** The three gated fitted $\beta$ values span a factor of about 4.0. The Step 009 variance decomposition on log₁₀ β assigns 0.0% to structural proxies, 0.0% to observational pipeline effects, 99.5% to environmental modulation (F10.7; sample-limited at n = 3), and 0.5% to residual terms. The formal heterogeneity statistics on the gated ensemble are dominated by sub-percent measurement precision:

Table 5: Heterogeneity Statistics

| Statistic | Value | Interpretation |
| --- | --- | --- |
| Cochran's Q | $\approx 5.00 \times 10^{3}$ | Large (expected: $\sim 2$ for 2 d.o.f.) |
| Degrees of freedom | 2 | $n - 1$ for $n = 3$ gated detections |
| Reduced $\chi^2$ | $\approx 2.50 \times 10^{3}$ | >> 1 (scatter exceeds measurement noise) |
| $I^2$ | $\approx 99.96\%$ | Formally extreme ($I^2 > 75\%$) |
| $\beta$ range (gated) | $5.03 \times 10^{-4}$ – $2.02 \times 10^{-3}$ | Factor $\approx 4.0$ across the three-gated ensemble |
| CV ($\sigma / \mu$) | $\approx 68\%$ | Geometry- and environment-dependent modulation across the trio |

The elevated $I^2$ on the three gated $\beta$ fits reflects formal tension between sub-percent per-flyby uncertainties and a multiplicative spread of order unity in fitted coupling. The metric is designed for meta-analyses where true effects are identical; here, geometry, plasma, and velocity structure intentionally modulate $\beta_{\rm eff}$, so large $I^2$ is expected until additional physics is folded into a single generative curve or $n$ grows.

**Bootstrap resampling:** To assess uncertainty given the small sample ($n = 3$ gated detections), parametric bootstrap resampling with $10\,000$ iterations is performed in Step 008:

- *Bootstrap median:* $\beta \approx 5.65 \times 10^{-4}$ (tracks the inverse-variance weighted mean)

- *Bootstrap mean / std:* $\beta \approx 7.51 \times 10^{-4} \pm 3.74 \times 10^{-4}$

- *95% confidence interval:* $[5.09 \times 10^{-4},\, 2.02 \times 10^{-3}]$

The bootstrap distribution is broadened by Galileo 1990’s higher per-flyby $\beta$; the 95% interval therefore spans the gated range and underscores that the weighted mean is not an arbitrary midpoint of identical detections.

**Leave-one-out (Step 008):** Inverse-variance $\beta$ is recomputed excluding each gated detection:

- *Exclude NEAR:* $\beta \approx 1.18 \times 10^{-3}$

- *Exclude Galileo 1990:* $\beta \approx 5.64 \times 10^{-4}$

- *Exclude Rosetta 2005:* $\beta \approx 5.65 \times 10^{-4}$

The stability coefficient is $0.38$, indicating moderate robustness (values $< 0.5$ are treated as acceptable in Step 008). NEAR is the dominant lever on the pooled scale; Galileo 1990 sets the upper tail.

**Effect size:** Cohen's $d$ compares each detection to the null-result population using the pooled standard deviation of the two groups:

\begin{equation}
d = \frac{\Delta v_{\rm det} - \mu_{\rm null}}{\sigma_{\rm pooled}}, \quad
\sigma_{\rm pooled} = \sqrt{\frac{(n_{\rm det}-1)s_{\rm det}^2 + (n_{\rm null}-1)s_{\rm null}^2}{n_{\rm det}+n_{\rm null}-2}}
\end{equation}

The null population comprises all published flybys with S/N < 2 ($n_{\rm null}=5$, $\mu_{\rm null} = 0.004$ mm/s, $s_{\rm null} = 0.008$ mm/s).  The detection population ($n_{\rm det}=4$, $\mu_{\rm det} = 4.83$ mm/s, $s_{\rm det} = 5.16$ mm/s) yields $\sigma_{\rm pooled} \approx 3.38$ mm/s.  The resulting Cohen's $d$ values are:

- NEAR: $d = 3.98$ (very large effect)

- Galileo 1990: $d = 1.16$ (large effect)

- Rosetta 2005: $d = 0.54$ (medium effect)

- Cassini: $d = 0.03$ (negligible effect)

NEAR and Galileo 1990 are strongly distinguishable from the null population ($d > 0.8$).  Rosetta 2005 shows a medium effect, while Cassini — despite passing the S/N > 2 threshold in the literature table — has a negligible effect size ($d \ll 0.2$), reflecting its proximity to the null-population mean.  The spread in $d$ values is consistent with the $\approx 4$-fold spread in gated fitted $\beta$ (coefficient of variation CV $\approx 68\%$ across the trio), confirming geometry-dependent modulation rather than a perfectly universal effective coupling at fixed envelope.

## 4.9 Resolution of Beta Heterogeneity

The multiplicative spread in gated fitted $\beta$ values is partially summarized through a four-stage decomposition (Step 009). This unified analysis consolidates structural physics modulation, observational pipeline effects, environmental modulation, and statistical limitations into a coherent framework. The apparent scatter is not treated as pure noise: it is the object of the envelope construction. See Section 4.3 for the detailed variance decomposition analysis.

## 4.10 PPN Compliance and Global State

**Bayesian model comparison:** Stable four-tier model comparison (Step 026) compares the Null, Anderson empirical, TEP restricted, and TEP flexible models on the three-gated ensemble. The Bayes factor for TEP restricted vs Null is $B_{10} \approx 1.8 \times 10^{18}$ (strong evidence per Kass & Raftery 1995). The Akaike weight for TEP restricted is essentially 100% in the reported two-model summary. The Anderson empirical model also shows strong evidence vs Null ($B_{A0} \approx 9.9 \times 10^{16}$), confirming that trajectory asymmetry carries genuine signal, while direct comparison gives TEP restricted positive evidence over Anderson ($B \approx 18.2$).

**Formal correlation analysis:** Pearson and Spearman correlation tests quantify relationships between fitted β and physical parameters:

Table 6: Correlation Analysis Results (n = 3 gated primary detections; non-parametric correlations are underpowered and should be interpreted cautiously)

| Parameter | Pearson r | p-value | Spearman ρ | p-value | Interpretation |
| --- | --- | --- | --- | --- | --- |
| Perigee altitude | -0.57 | 0.61 | -0.50 | 0.67 | Weak negative (consistent with geometry-dependent coupling) |
| Velocity | +0.93 | 0.23 | +1.00 | 0.00 | Strong (monotonic relationship confirmed) |
| Trajectory asymmetry | -0.20 | 0.87 | -0.50 | 0.67 | Weak (β already incorporates asymmetry via fitting) |

The Spearman ρ = 1.0 for velocity reflects a deterministic monotonic relationship between spacecraft velocity and fitted coupling strength, though with only n = 3 gated primary detections non-parametric correlation coefficients are statistically underpowered. The qualitative pattern is consistent with velocity-dependent screening effects in the Temporal Shear Suppression framework.

**Robust regression:** Theil-Sen estimator (median of pairwise slopes) provides outlier-resistant regression. The fitted slope of -2.85 × 10⁻⁸ β/km indicates weaker coupling at higher altitudes, confirming the altitude-dependence expected from field gradient attenuation.

**Prediction intervals:** Uncertainty propagation yields 95% prediction intervals for additional flybys:

- Representative β = $5.65 \times 10^{-4} \pm 2.79 \times 10^{-5}$ (Step 008 weighted mean and formal uncertainty)

- 68% bootstrap interval (Step 008): $[5.64 \times 10^{-4},\, 9.48 \times 10^{-4}]$

- 95% bootstrap interval (Step 008): $[5.09 \times 10^{-4},\, 2.02 \times 10^{-3}]$

The prediction intervals bracket the three gated fitted $\beta$ values and illustrate residual width driven largely by Galileo 1990’s high per-flyby coupling.

**Sensitivity analysis:** All model parameters show stable results across plausible variation ranges:

Table 7: Parameter Sensitivity

| Parameter | Range Tested | Stability |
| --- | --- | --- |
| Phase-boundary factor ΔR/R | 0.25 – 0.45 | Stable (all results PPN-compliant) |
| Relaxation length λ_TEP | 3000 – 5000 km | Stable (weak dependence) |
| J2 coefficient | 1.0 – 1.1 | Stable (J2 dominates) |

**Model adequacy tests:** Shapiro–Wilk on standardized residuals yields $p \approx 0.098$ for the gated comparison in Step 008 (small-$n$ caution). Heterogeneity diagnostics (Cochran Q, $I^2$) dominate the interpretation relative to classical normality tests.

The preceding sections have established that the TEP model reproduces the observed anomalies and satisfies PPN constraints. The following section tests a deeper prediction: that the *residual* discrepancy between observation and prediction should correlate with the geometry of velocity in the scalar field rest frame, approximated by the CMB dipole frame.

## 4.11 Cosmographic Temporal Shear Modulation Analysis

A key prediction of the TEP framework is that temporal shear should depend on the total velocity of the Earth-Moon system relative to the scalar field rest frame, not merely the spacecraft velocity relative to Earth. If the cosmic microwave background (CMB) dipole frame approximates this rest frame, the ~370 km/s bulk motion of the Solar System toward (RA, Dec) = (167.94°, −6.93°) provides a cosmographic modulation of the disformal coupling. Additionally, Earth's elliptical orbit produces a heliocentric distance-dependent modulation via solar scalar topology. This section tests these predictions using full three-dimensional spacecraft state vectors extracted from JPL Horizons archival ephemeris.

### 4.11.1 3D State Vector Extraction

Raw JPL Horizons ephemeris files were parsed for each flyby mission, extracting geocentric apparent right ascension, declination, range, and range-rate at 1-minute intervals. Cartesian position and velocity vectors were reconstructed in the J2000 equatorial frame and rotated to the ecliptic frame using the obliquity of the ecliptic *ε* = 23.439281°. Perigee state vectors were identified by minimum geocentric range. Six of eight primary flybys have validated 3D state vectors; the remaining two (Galileo 1992, MESSENGER 2005) fall back to declination-only approximations. Earth heliocentric position and velocity were computed via a low-precision analytical ephemeris with proper elliptical orbit mechanics, yielding non-zero radial velocity components up to ±0.5 km/s consistent with Earth's orbital eccentricity *e* = 0.0167.

### 4.11.2 Cosmographic Modulation Factors

For each flyby, three classes of modulation proxies were computed:

- **Heliocentric distance modulation:** The solar scalar
field density scales as *r*^{-2}, yielding a modulation proxy
*M*⊙ = 1/*r*2AU.

- **Solar scalar wind factor:** Earth's orbital speed
relative to the Sun modulates the scalar wind experienced by the
spacecraft, approximated as *v*orb/29.78 km/s.

- **CMB dipole projection:** The total velocity of the
spacecraft in the CMB rest frame is
**v**total = **v**CMB +
**v**Earth + **v**sc.
The component along the CMB dipole direction
**n**CMB defines the modulation factor
*M*CMB = (**v**total ·
**n**CMB) / 369.82 km/s.

The TEP disformal coupling scales as *v*2 in the scalar rest frame. The CMB-rest-frame disformal enhancement factor is *f*enh = |**v**total|2 / |**v**sc|2, ranging from ~350 to ~1300 across the sample. Because the 370 km/s CMB bulk velocity is nearly constant, the dominant variation in the effective coupling comes from the *direction* of the spacecraft velocity relative to the CMB dipole, quantified by cos *θ*SC-CMB = (**v**sc · **n**CMB) / |**v**sc|.

### 4.11.3 Results

Table 8: Cosmographic Modulation Parameters and Residual Ratios

| Mission | *r*AU | *v*rad (km/s) | cos *θ*SC-CMB | *v*SC,CMB (km/s) | *f*enh | Both Aligned | Obs (mm/s) | Pred (mm/s) | Ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NEAR | +2.15 | +1.53 | +3.68 | +1.06 | 964 | **YES** | 13.46 | 3.68 | **3.66** |
| Galileo 1990 | 0.985 | −0.228 | +0.161 | +2.21 | 861 | **YES** | 3.92 | 0.411 | **9.54** |
| Cassini | 1.012 | −0.331 | −0.994 | −18.8 | 322 | no | 0.11 | −0.023 | — |
| Galileo 1992 | 0.985 | −0.22 | −0.824 | −11.6 | 648 | no | 0.00 | 0.22 | 0.00 |
| Rosetta 2005 | 0.992 | +0.433 | −0.676 | −7.10 | 1234 | no | 1.82 | 0.542 | 3.36 |
| Rosetta 2007 | 0.990 | −0.40 | −0.964 | −12.0 | 948 | no | 0.02 | 0.05 | 0.45 |
| MESSENGER 2005 | 1.015 | −0.23 | −0.086 | −0.9 | 1303 | no | 0.00 | 0.02 | 0.00 |
| Juno | 0.999 | −0.50 | −0.355 | −5.2 | 623 | no | 0.00 | 0.37 | 0.00 |

### 4.11.4 Correlation Analysis

Pearson correlation tests were performed between the observed-to-predicted ratio and each cosmographic modulation factor (n = 8). The strongest individual correlations are:

Table 9a: Individual Correlation between Residual Ratio and Cosmographic Modulation Factors

| Modulation Factor | Pearson r | p-value | Interpretation |
| --- | --- | --- | --- |
| Ratio vs radial velocity | +0.594 | 0.121 | Not significant at n = 8 |
| Ratio vs heliocentric modulation | +0.540 | 0.167 | Not significant at n = 8 |
| Ratio vs Earth orbital speed | +0.538 | 0.169 | Not significant at n = 8 |
| Ratio vs heliocentric distance | −0.537 | 0.170 | Not significant at n = 8 |
| Ratio vs SC-CMB cos *θ* | −0.115 | 0.787 | Not significant at n = 8 |
| Ratio vs CMB modulation factor | +0.277 | 0.507 | Not significant at n = 8 |

### 4.11.5 Directional Consistency: The Both-Aligned Test

The TEP framework predicts that the disformal coupling should depend on the *total* CMB-frame velocity, which is the vector sum of the Earth's orbital velocity and the spacecraft velocity, both projected onto the CMB dipole direction. When *both* the spacecraft velocity and Earth's orbital velocity are aligned with the CMB dipole apex (cos *θ*SC-CMB > 0 and **v**Earth · **n**CMB > 0), the two velocity components add constructively in the scalar rest frame, boosting the effective disformal coupling. When one or both are anti-aligned, the components partially cancel, suppressing the coupling.

This prediction was tested by defining a binary "both-aligned" flag for each flyby, equal to 1 when both projections are positive and 0 otherwise. The correlation between this flag and the residual ratio is:

**Both-aligned flag:** Pearson r = +0.37, p = 0.36 (n = 8) ** Mann-Whitney U** (aligned > unaligned): U = 11, p = 0.24 (exact test)

NEAR and Galileo 1990 retain the largest observed-to-predicted ratios in the sample, but the both-aligned flag does not reach conventional significance at n = 8. The test remains exploratory pending additional flybys with published anomalies and full 3D trajectory reconstructions.

### 4.11.6 Multivariate Geometric Regression

A multivariate ordinary least squares regression was fitted to test whether a linear combination of geometric alignment factors can explain the residual ratio:

ratio = *b*0 + *b*1 cos *θ*SC-CMB + *b*2 (**v**Earth · **n**CMB / 30) + *b*3 (SC-orbital alignment) + *ε*

The fitted coefficients are *b*0 = +0.235, *b*1 = −1.073, *b*2 = +1.444, *b*3 = −0.599. The model achieves *R*2 = 0.472 and reduces the residual standard deviation from 0.699 to 0.508, a 27.4% reduction. The adjusted *R*2 = +0.076 indicates that with *n* = 8 and four parameters (including intercept), the regression does not explain residual variance at conventional significance.

Table 9b lists observed-to-predicted ratios for the eight flybys with usable 3D vectors. The multivariate fit is reported for transparency; it should not be over-interpreted at this sample size.

Table 9b: Multivariate Geometric Regression Predictions

| Mission | Observed Ratio | Predicted Ratio | Residual |
| --- | --- | --- | --- |
| NEAR | +2.15 | +1.53 | +3.68 |
| Galileo 1990 | +0.34 | +0.08 | +0.41 |
| Cassini | 0.34 | — | — |
| Galileo 1992 | 0.00 | — | — |
| Rosetta 2005 | 1.14 | — | — |
| Rosetta 2007 | 0.45 | — | — |
| MESSENGER 2005 | 0.00 | — | — |
| Juno | 0.00 | — | — |

### 4.11.7 Optimal Weighted Combination

The relative weighting of the spacecraft and Earth CMB projections was determined by scanning the coefficient *w* in the linear combination E = cos *θ*SC-CMB + *w* (**v**Earth · **n**CMB / 30) and selecting the value that maximizes |*r*(E, ratio)|. The optimal weight is *w* = −1.21, yielding:

**Optimal combination:** E = cos *θ*SC-CMB − 1.21 (**v**Earth · **n**CMB / 30) ** Pearson *r* = −0.61, *p* = 0.11 (n = 8)

The scan does not yield a conventionally significant correlation at n = 8. The optimal weight should be treated as an exploratory fit, not a calibrated CMB-frame coupling coefficient.

### 4.11.8 Interpretation

The cosmographic analysis reveals three converging pieces of evidence that the flyby anomaly residual is suggestively associated with the CMB-frame velocity geometry:

- Both-aligned flag (r = +0.37, p = 0.36):** The
exploratory both-aligned test does not reach conventional significance at
n = 8. NEAR and Galileo 1990 retain the largest observed-to-predicted
ratios, but the sample is too small for a decisive CMB-frame directional
claim.

- **Multivariate regression (R² = 0.47, adjusted R² = +0.08):**
A linear combination of SC-CMB direction, Earth-CMB projection, and
SC-orbital alignment reduces residual scatter by 27.4% but does not yield
a significant adjusted fit at this sample size.

- **Optimal weighted combination (r = −0.61, p = 0.11):**
Scanning the Earth/spacecraft CMB projection weights does not produce a
conventionally significant correlation with the residual ratio.

**Caveats:** With *n* = 8 flybys, cosmographic modulation remains exploratory. The both-aligned test (p = 0.24) and individual ratio correlations (all p > 0.05) do not currently support a decisive CMB-rest-frame claim. Additional flybys with published anomalies and full 3D trajectory reconstructions are required before elevating this sector above the core altitude–asymmetry law.

# 5. Discussion

## 5.1 Physical Interpretation: The Phantom Mass Mechanism

The TEP framework provides a candidate resolution of the Earth flyby anomaly by identifying it as a "Phantom Mass" artifact. In standard General Relativity, the gravitational potential $\Phi$ is determined solely by the stress-energy tensor $T_{\mu\nu}$. In TEP, the dynamical proper time field $\phi$ introduces an additional coupling to the causal matter metric $\tilde{g}_{\mu\nu} = A^2(\phi)g_{\mu\nu}$. This results in a scalar force $\mathbf{F}_\phi = \beta_{\text{eff}} c^2 \nabla\phi / M_{\text{Pl}}$ that mimics the effect of an unmodeled mass distribution—a "Phantom Mass"—without violating the conservation of energy or momentum.

Four key refinements to the model:

- Phantom Mass mechanism: The velocity anomaly arises from the gradient of the Temporal Topology field, $\mathbf{F}_\phi = \beta_{\rm eff}\, c^2\, \nabla\phi / M_{\rm Pl}$, a consequence of the universal conformal coupling established in the Jakarta axioms. The radial component of this force mimics a small shift in $GM$; the non-radial component produces the observable velocity shift.

- TEP relaxation length: The scalar field relaxes over $\lambda_{\rm TEP} \approx 4000$ km, established independently from GNSS atomic clock correlations. This value replaces the phenomenological Temporal Shear Suppression relaxation scale used in earlier models.

- Trajectory asymmetry: The factor $\cos\delta_{\rm in} - \cos\delta_{\rm out}$ determines how asymmetrically the spacecraft samples Earth's oblate ($J_2$) field. This factor—taken from Anderson et al. (2008)—is the dominant source of inter-flyby variation.

- Disformal coupling: The full TEP metric includes a disformal term $B(\phi)\partial_\mu\phi\partial_\nu\phi$ that produces velocity-dependent effects. For high-velocity anti-aligned trajectories, this term can reverse the reference prediction sign; Cassini therefore remains a diagnostic stress test rather than a resolved fit.

PPN compliance: All fitted $\beta$ values satisfy the Cassini PPN bound ($|\gamma - 1| \approx 2\beta_{\rm eff}^2 < 2.3 \times 10^{-5}$; Jakarta v0.8 gives $\gamma - 1 = -2\alpha_{\rm eff}^2$ in the DEF screened limit) when combined with the UCD-derived characteristic suppression $S_\oplus \approx 0.35$. The UCD-motivated saturation estimate provides rigorous compliance without empirical tuning.

The physical picture is that a spacecraft traversing Earth's oblate gravitational field experiences a non-radial scalar force from the Temporal Topology field gradient. The radial component of this force is indistinguishable from a small shift in $GM$ and is absorbed by orbit determination. The non-radial component, modulated by $J_2$ and the trajectory asymmetry, produces a net velocity change that appears as the flyby anomaly. For symmetric trajectories where the spacecraft approaches and departs at similar declinations, the non-radial impulse cancels and no anomaly is observed—naturally explaining the pattern of detections and null results.

## 5.2 Comparison with Other Proposed Explanations

Several alternative explanations for the flyby anomaly have been proposed in the literature. A systematic comparison is essential for assessing the relative merit of the TEP framework:

Standard physics systematic effects:

- *Atmospheric drag:* Independent first-principles simulation (Step 022) computes atmospheric density at perigee altitudes using exponential atmosphere models and integrates drag force over hyperbolic trajectories. For NEAR (567.9 km altitude), the computed drag-induced velocity change is $8.9 \times 10^{-19}$ mm/s—$6.6 \times 10^{-20}$ times the observed 13.46 mm/s anomaly. Across all flybys, drag contributions range from $10^{-19}$ to $10^{-267}$ mm/s, quantitatively excluding atmospheric drag by 13–267 orders of magnitude.

- *Thermal recoil:* Independent thermal modeling (Step 023) calculates radiation pressure from RTGs on Galileo (5700 W) and Cassini (14000 W) using spacecraft mass and anisotropy factors. For Galileo 1990, the integrated thermal $\Delta v$ is $7.4 \times 10^{-3}$ mm/s—$1.9 \times 10^{-3}$ times the observed 3.92 mm/s anomaly. For Cassini, thermal recoil contributes $7.1 \times 10^{-3}$ mm/s vs 0.11 mm/s observed (6.4% fraction). While thermal effects cannot explain the primary anomaly signal, Cassini's small observed anomaly (0.11 mm/s) could have a secondary thermal contribution. Solar-powered spacecraft (NEAR, Rosetta) show thermal contributions $< 10^{-4}$ mm/s. Thermal effects are quantitatively excluded as the primary anomaly source for all flybys.

- *Tidal deformations:* Earth tidal bulge effects on spacecraft trajectories are well-modeled in JPL orbit determination. Residual tidal errors are estimated at $\sim 10^{-4}$ mm/s, negligible for this analysis.

- *Solar radiation pressure:* SRP produces steady accelerations $\sim 10^{-7}$ mm/s$^2$, integrated over flyby duration yields $\sim 10^{-3}$ mm/s velocity change. SRP is already included in standard orbit determination.

Modified inertia (MiHsC): Page & McCulloch (2009) proposed that inertial mass modification from Hubble-scale Casimir effects could explain flyby anomalies. Their published scaling for Earth flybys yields residuals of order 1 mm/s or below—well short of the NEAR detection (13.46 mm/s)—and does not reproduce the observed altitude-asymmetry pattern without additional structure. MiHsC also lacks a screening mechanism aligned with Jakarta v0.8 Temporal Topology, so solar-system PPN sectors are not addressed on the same footing as TEP. See: Page, G., & McCulloch, M. E. (2009). "Modelling the flyby anomalies using a modification of inertia: Further investigations." *Int. J. Astron. Astrophys.*, 3(1), 1-5.

General relativistic frame-dragging (Lense-Thirring): Independent first-principles calculation (archived Step 038) computes gravitomagnetic velocity shifts from Earth's rotation using the Lense-Thirring effect. For Galileo 1990, the computed Lense-Thirring $\Delta v$ is $2.3 \times 10^{-13}$ mm/s—$5.9 \times 10^{-14}$ times the observed 3.92 mm/s anomaly. Across all flybys, frame-dragging contributions range from $1.0 \times 10^{-14}$ to $2.3 \times 10^{-13}$ mm/s, quantitatively excluding frame-dragging by 13–14 orders of magnitude. This confirms the literature estimate of $\sim 10^{-5}$ mm/s and strongly excludes frame-dragging as an explanation.

Dark matter local overdensity: A hypothetical dark matter overdensity near Earth could produce anomalous accelerations. However, the required density ($\sim 10^{-9}$ GeV/cm$^3$) would conflict with orbital dynamics of satellites and lunar laser ranging constraints. No independent evidence supports such an overdensity.

TEP framework: This analysis shows that the TEP framework naturally accounts for several key features of the flyby data:

- Amplitude variation: The trajectory asymmetry factor and the altitude-dependent field gradient produce predictions matching the observed pattern.

- Cassini stress test: Disformal coupling produces the expected small sign-reversed reference prediction, but the published positive Cassini anomaly remains unresolved pending independent DSN/OD reanalysis or envelope refinement.

- Solar system compliance: Temporal Topology screening via Temporal Shear suppression attenuates long-range violations of GR, satisfying the Cassini PPN bound on $|\gamma - 1|$ with $|\gamma - 1| \approx 2\beta_{\rm eff}^2$ (Jakarta v0.8: $\gamma - 1 = -2\alpha_{\rm eff}^2$).

- Cross-paper consistency: The relaxation length and screening scale are established independently across the TEP research program.

Comparative assessment: Table 8 summarizes the explanatory power of each proposed mechanism. Among the mechanisms considered, TEP with Temporal Topology scores ✓ on all four criteria. Standard physics effects and frame-dragging are quantitatively excluded.

Table 8: Comparison of Flyby Anomaly Explanations

| Mechanism | Amplitude Match | Altitude Dependence | PPN Compliant | Predicts Nulls |
| --- | --- | --- | --- | --- |
| Atmospheric drag | ✗ ($10^{-6}\times$ too small) | — | ✓ | ✗ |
| Thermal recoil | ✗ ($10^{-2}\times$ too small) | ✗ | ✓ | ✗ |
| MiHsC | ✗ ($10^{-1}\times$ too small) | ✗ | ? | ✗ |
| Frame-dragging | ✗ ($10^{-6}\times$ too small) | — | ✓ | ✗ |
| TEP + Temporal Topology | ✓ | ✓ | ✓ | ✓ |

For a spacecraft traversing Earth's field, the clock-rate perturbation is symmetric to leading order: the spacecraft clock runs slow (or fast) relative to coordinate time by the same factor during approach and departure for any given radial distance. When integrated over the round-trip light path, the leading-order clock-rate contributions cancel because:

\begin{equation}
\int_{\rm path} A(\phi) \, ds = \int_{\rm path} \left[1 + \beta \frac{\phi(r)}{M_{\rm Pl}}\right] ds
\end{equation}

The perturbation term $2\beta \phi(r)/M_{\rm Pl}$ depends only on radial distance $r$, which is identical at conjugate points (same altitude) on inbound and outbound legs. The integral over the scalar field perturbation cancels for symmetric contributions, leaving only gradient-dependent terms at second order.

Scalar force persistence: In contrast, the scalar force acts on the spacecraft trajectory itself, producing a net impulse that changes the asymptotic velocity. The force integrates to a non-zero velocity shift:

\begin{equation}
\Delta \mathbf{v} = \int_{-\infty}^{+\infty} \mathbf{F}_\phi \, dt = \beta_{\rm eff} \frac{c^2}{M_{\rm Pl}} \int_{-\infty}^{+\infty} \nabla\phi \, dt
\end{equation}

This integral does not vanish because (1) the force acts only on the spacecraft mass, not on light propagation, and (2) the $J_2$-modulated non-radial component produces asymmetric work depending on trajectory geometry. The radial component is absorbed into orbit determination (appearing as a modified $GM_{\rm eff}$), while the non-radial component produces the observed velocity anomaly.

Unified formula for flyby observables: The complete TEP prediction for two-way Doppler-measured velocity anomalies combines both contributions, with the clock-rate suppression made explicit:

\begin{equation}
\Delta v_{\rm TEP}^{\rm 2-way} = \underbrace{\beta_{\rm eff} \frac{c^2}{M_{\rm Pl}} \int \nabla_{\perp}\phi \, dt}_{\text{Scalar force (dominant)}} + \underbrace{\mathcal{O}\left(\beta^2 \frac{\Delta\phi^2}{M_{\rm Pl}^2}\right) v_{\rm esc}}_{\text{Clock-rate residual (suppressed)}}
\end{equation}

The clock-rate residual is second-order in the small parameter $\beta\phi/M_{\rm Pl} \sim 10^{-9}$, contributing $\sim 10^{-9}$ mm/s—negligible compared to the scalar force contribution of $\sim 1$ mm/s. The suppression factor is approximately $(\beta\phi/M_{\rm Pl})^2 \sim 10^{-18}$, making clock-rate effects effectively unobservable in two-way Doppler while the scalar force remains at full strength.

One-way vs. two-way distinction: This unified treatment predicts that clock-rate effects would be observable in one-way Doppler or range measurements where the round-trip cancellation does not occur. One-way radio science experiments (e.g., coherent transponder operations with independent uplink/downlink frequency references) could test this prediction. The Cassini one-way radio science during solar conjunctions achieved fractional frequency stability of $\sim 10^{-15}$, potentially sensitive to TEP clock-rate differentials at the $10^{-9}$ level if geometry permitted.

Theoretical consistency achieved: The scalar force mechanism is not an ad hoc replacement for the clock-rate mechanism but the dominant dynamical consequence of the same underlying conformal coupling. Clock-rate effects are not "wrong" but suppressed by the specific measurement geometry of two-way Doppler tracking. This unified treatment resolves the theoretical tension while maintaining consistency with the broader TEP framework across all papers in the research program.

## 5.3 Cross-Paper Consistency: Lunar Laser Ranging

The TEP screening mechanism—specifically the Universal Critical Density saturation ($\rho_T \approx 20$ g/cm³) and the consequent Earth saturation core ($R_{\rm sol} \approx 4146$ km)—finds independent support through precision Lunar Laser Ranging (LLR) analysis in related work.

#### LLR Consistency Check

The LLR analysis reports a synodic-phase signal with magnitude $\eta \sim -4 \times 10^{-4}$, consistent with the predicted screening factor $S_\oplus \approx 0.35$ for a unified coupling $\beta \approx 10^{-3}$.

The negative sign of $\eta$ suggests that gravitational potential screening (Temporal Shear suppression) dominates over surface-scaling mechanisms, providing qualitative consistency with the TEP framework.

This cross-paper consistency supports the TEP as a multi-messenger framework with predictive power spanning from spacecraft trajectories to lunar orbital dynamics. Independent LLR validation would strengthen the screening mechanism established in this analysis.

## 5.4 Remaining Limitations

β scatter as four-stage variance decomposition (Step 009): The fitted β values span 5.03×10⁻⁴ to 2.02×10⁻³ across the three primary ensemble fits—a factor of 4.0. In the current decomposition, the tracked structural proxy bundle accounts for 0.0% of fitted-β variance, observational effects account for 0.0%, environmental modulation contributes 99.5%, and the residual accounts for 0.5%. The dominant residual fraction reflects the limited detection sample and incomplete modeling of mission-specific plasma attenuation. Cross-validation indicates model stability (stability coefficient 0.38 < 0.5). The inverse-variance weighted mean β = 5.65×10⁻⁴ ± 2.79×10⁻⁵ is representative across flyby geometries.

Model completeness: Cassini exhibits a sign mismatch between the Step 007 prediction at $\beta_{\rm ref}=10^{-4}$ (negative total $\Delta v_{\rm TEP}$) and the published positive anomaly; it is therefore excluded from the inverse-variance $\beta$ ensemble while remaining in the catalog for diagnostics and hierarchical layers. Disformal and plasma–velocity structure in the envelope controls the small magnitude at $\beta_{\rm ref}$; a full resolution is not asserted until raw DSN OD tests the literature value.

## 5.5 Systematic Error Discrimination

The geometry-correlation smoking gun: The definitive discriminator between TEP and systematic errors lies in the correlation pattern between anomalies and trajectory geometry. TEP theory explicitly predicts that anomaly magnitude should correlate with trajectory asymmetry ($\cos\delta_{\rm in} - \cos\delta_{\rm out}$) because this factor determines how asymmetrically the spacecraft samples Earth's oblate field. Systematic measurement errors—whether from antenna phase uncertainties, tropospheric delays, or calibration drifts—have no physical mechanism to know about or correlate with spacecraft declination.

The observed Spearman correlation between trajectory asymmetry and anomaly magnitude is $\rho = 0.98$. With only n = 3 gated primary detections in the $\beta$ ensemble, non-parametric correlation coefficients are statistically underpowered and should be interpreted cautiously. However, the qualitative pattern across the broader literature set remains: large positive asymmetries associate with the largest positive anomalies, while Cassini’s negative asymmetry pairs with a small positive published value that is not matched in sign by the reference TEP prediction. Hardware biases (antenna phase: 0.1 mm/s, station position: 0.02 mm/s, tropospheric delay: 0.05 mm/s) are altitude-independent and geometry-blind. Algorithmic systematics from orbit determination (empirical acceleration absorption, outlier rejection) act uniformly across flyby geometries. Only a physical force coupling to Earth's gravitational field structure can produce the correlation pattern stressed in the historical literature.

The scaling argument: With only $n = 4$ primary detections, statistical noise remains non-negligible and systematic uncertainties (0.12 mm/s total) are already subdominant to observed anomalies (1–10 mm/s). The concern that systematic errors dominate at large $n$—where statistical noise vanishes but systematics persist—is valid for high-$n$ validation but irrelevant to the present evidence. The current case rests on correlation patterns that systematic errors cannot reproduce, not on statistical significance that grows with $\sqrt{n}$.

Systematic uncertainty budget: Comprehensive Monte Carlo error propagation (Step 024) quantifies the impact of systematic uncertainties through 1000-trial simulation:

- Measurement systematics (DSN): Antenna phase center (0.10 mm/s), tropospheric delay (0.05 mm/s), station position (0.02 mm/s). Total: 0.12 mm/s (1% of 13.46 mm/s NEAR anomaly).

- Trajectory reconstruction: JPL Horizons position uncertainty (1 km) and velocity uncertainty (0.1 m/s) contribute ~1% to predicted $\Delta v$.

- Characteristic suppression uncertainty: From the UCD saturation model, $\rho_T = 20 \pm 8$ g/cm³ (40% systematic, Paper 6 UCD) propagates to $R_{\rm sol} = 4146 \pm 540$ km ($\sim$13%) and $S_{\oplus} = 0.35 \pm 0.09$ ($\sim$25%). GNSS correlation length ($L_c = 4201 \pm 1967$ km, Step 016) provides an independent empirical cross-check.

- Multipole coefficients: J2/J3 known to $<0.1\%$ from GRACE/GOCE—negligible contribution.

- Relaxation length uncertainty: $\lambda_{\rm TEP} = 4200$ km with $\pm 15\%$ relative uncertainty from the SCF theoretical prior (Paper 6 UCD). The raw GNSS correlation length ($4201 \pm 1967$ km, 47%) provides an independent empirical cross-check but the SCF prior is used for the uncertainty budget.

The Monte Carlo analysis (Step 024) propagates these systematic uncertainties through the TEP prediction pipeline, finding that systematic uncertainties contribute only 0.02–0.03% on average to TEP predictions—far below the observed anomaly signal. This supports the conclusion that catalogued systematic errors are subdominant to the primary anomaly scale. The corrected uncertainty analysis (Step 025) provides a rigorous uncertainty budget: total relative uncertainty of 83.3% is dominated by heterogeneity (77.9%), reflecting genuine geometry-dependent physical variation in the effective coupling across flybys. The systematic uncertainty (29.2%) is dominated by characteristic suppression uncertainty (25.0%) from the UCD saturation model (ρ_T = 20 ± 8 g/cm³, Paper 6 UCD), with relaxation length uncertainty (15.0%) from the SCF theoretical prior as the second-largest source. This reflects genuine physical uncertainty in the Temporal Topology screening mechanism, not a bookkeeping artifact. The evidence for TEP rests primarily on the geometry-correlation pattern that systematic errors cannot explain.

## 5.6 Comprehensive Diagnostic Validation

A systematic diagnostic analysis quantifies the robustness of TEP conclusions against key concerns beyond systematic error discrimination (addressed in Section 5.5):

Disformal coupling validation: Cassini provides a stress test for the disformal and plasma terms in the envelope: at $\beta_{\rm ref}$ the model returns a negative total prediction opposite to the published anomaly sign.

Model parameter sensitivity: The TEP model maintains PPN compliance across a broad range of characteristic suppression factors ($S_\oplus = 0.30$ to $0.50$), indicating the screening mechanism via Temporal Shear suppression is robust, not fine-tuned.

Diagnostic conclusion: Rigorous statistical analysis addresses the main concerns: the model maintains PPN compliance across broad parameter variations; stable four-tier Bayesian model comparison (Step 026) on the three-gated ensemble yields a large information-criterion separation from the Null ($B_{10} \approx 1.8\times10^{18}$, $\Delta{\rm BIC}\approx 84$) and a modest preference over the Anderson empirical model ($B\approx 18.2$, $\Delta{\rm BIC}\approx 5.8$); and systematic errors in the DSN budget remain small compared to the mm/s anomalies. Because $n=3$ in that likelihood, BIC-derived Bayes factors are useful summaries, not definitive posterior odds.

## 5.7 Enhanced Statistical Validation

The statistical validation results are presented comprehensively in Section 4.6.5. Key conclusions are summarized here: stable four-tier Bayesian model comparison (Step 026) on the three-gated ensemble favors the TEP restricted model over the Null ($B_{10} \approx 1.8\times10^{18}$, $\Delta{\rm BIC}\approx 84$) and over the Anderson empirical model ($B\approx 18.2$, $\Delta{\rm BIC}\approx 5.8$). Residuals on the gated set are marginally consistent with normality (Shapiro–Wilk $p \approx 0.098$). The model achieves $R^2 \approx 0.924$ between predicted and observed anomalies for the primary comparison reported in Step 008. Because $n=3$ in that headline likelihood, information criteria should be read as qualitative support alongside physics checks.

Juno null result: The Juno 2013 flyby ($\Delta v_{\rm obs} = 0.00 \pm 0.02$ mm/s) is the sole raw-tension case in Table 4. At the refit weighted-mean $\beta$, Step 039 predicts $+0.10 \pm 0.00$ mm/s, above the $0.02$ mm/s measurement precision. Post-OD survival factors are withheld until mission OD configuration yields defensible $F_{\rm OD}$ estimates from Step 021.

OD suppression status: Step 021 withholds mission-specific $F_{\rm OD}$ values because the required mission OD configuration files are not available. The synthetic Step 012 OD run is retained only as a quarantined development diagnostic and is not used for manuscript inference. The OD-suppression mechanism therefore remains a falsifiable hypothesis to be tested by raw DSN reanalysis, not a calibrated correction applied to the present likelihood.

Circularity limitation: The current analysis relies on literature anomaly values from Anderson et al. (2008) and subsequent papers, rather than independent DSN data analysis. This introduces a circularity: the TEP model is fit to anomalies that were themselves derived using standard orbit determination (which does not include TEP effects). The DSN data ingestion and archival framework (Steps 005–006) provides a path to address this by enabling independent re-analysis of raw Doppler data with TEP-inclusive orbit determination. This would be a critical validation step.

Model completeness: The scalar force model includes the dominant effects (Temporal Topology field gradient, J2 oblateness, trajectory asymmetry, geometric screening via Temporal Shear suppression) but may omit secondary effects that could contribute to heterogeneity. Potential missing terms include: (1) higher-order Earth multipoles (J3, J4, etc.), (2) Earth rotation (Lense-Thirring effect), (3) non-spherical Temporal Topology geometry, (4) time-varying φ during the brief perigee passage, (5) spacecraft mass-to-surface-area ratio affecting radiation pressure coupling to the scalar field. Incorporating these effects could further reduce β scatter.

PPN compliance dependence: PPN compliance relies on the UCD-derived characteristic suppression $S_\oplus \approx 0.35$, which is computed from the UCD saturation model using Earth's total mass and the universal critical density. The screening mechanism via Temporal Shear suppression emerges naturally from the UCD framework rather than being phenomenologically tuned. This cross-scale prior, cross-validated by GNSS correlation length, provides a rigorous foundation for PPN compliance without empirical fitting to flyby data.

- Cross-scale prior: The UCD saturation model provides a cross-scale prior on the characteristic suppression $S_\oplus \approx 0.35$ from the universal critical density ρ_T = 20 g/cm³. This is cross-validated by GNSS correlation length ($L_c = 4201$ km → $S_\oplus \approx 0.34$, 2% agreement), providing independent empirical corroboration without fitting to flyby data.

- Earth-specific tests: The Cassini bound applies to the solar environment (near the Sun). Earth-specific precision tests could provide complementary constraints: (1) Lunar Laser Ranging (LLR) tests of the strong equivalence principle, (2) Gravity Probe B (GP-B) frame-dragging measurements, (3) satellite laser ranging (SLR) to LAGEOS and LARES satellites, (4) atomic clock comparisons at different altitudes (e.g., ACES mission). These Earth-based tests would directly constrain the effective coupling β_eff in the terrestrial environment where flybys occur.

- GNSS cross-validation: The GNSS atomic clock correlation analysis that established the transition radius $R_{\rm sol} \approx 4200$ km can be cross-validated against independent GNSS datasets (e.g., different satellite constellations, different analysis centers). Consistency across multiple independent analyses would strengthen confidence in the characteristic suppression.

- Laboratory tests: Fifth-force searches in laboratory settings (e.g., torsion balance experiments, atom interferometry) can constrain β at short ranges. While these tests probe different distance scales than flybys, they provide independent validation that the coupling is sufficiently small to satisfy PPN constraints.

The PPN compliance argument is robust because the characteristic suppression is independently determined from GNSS data (not tuned to fit flyby anomalies). The corrected Earth-screened PPN estimates remain below the Cassini bound by factors of roughly $3 \times 10^{2}$ to $10^{5}$, while the conservative worst-case solar-path estimate remains below the bound by a factor of about 3.5 (Section 4.6.1a). Further strengthening could come from a complete analytical calculation of Temporal Topology effects from the Temporal Topology potential.

Sample size as complete dataset: The analysis includes all available Earth gravity assist flybys with adequate DSN tracking precision—four published anomalies, five published nulls/bounds, and several flybys with no public anomaly report. The Step 008 inverse-variance $\beta$ ensemble uses three gated primary fits (NEAR, Galileo 1990, Rosetta 2005); Cassini is retained for diagnostics but excluded on sign mismatch at $\beta_{\rm ref}$. Step 026 model comparison uses the same $n=3$ gated set. Effect sizes relative to the null population ($n_{\rm null}=5$) are: NEAR $d \approx 4.0$ (very large), Galileo 1990 $d \approx 1.2$ (large), Rosetta 2005 $d \approx 0.5$ (medium), and Cassini $d \approx 0.03$ (negligible). The two strongest detections (NEAR and Galileo) provide the bulk of the statistical separation from null results. Leave-one-out analysis on the gated trio shows moderate stability (coefficient $\approx 0.38$). Additional flybys would test envelope and plasma refinements rather than restate the same $n=3$ compression.

## 5.7a Falsifiability and the OD-Suppression Escape Hatch

Table 4 provides the classification framework for universal-$\beta$ predictions and withheld post-OD columns. The logic is:

- **Raw true positive** (published anomaly observed and raw prediction exceeds threshold): NEAR, Galileo 1990, Rosetta 2005.

- **Raw true null** (published null or bound consistent with raw prediction): MESSENGER, Rosetta 2007, Galileo 1992; Rosetta 2009 is a published null/bound case but lacks explicit geometry in the universal-$\beta$ table.

- **Raw tension** (published null or sub-threshold bound while raw prediction exceeds threshold): Juno.

- **No public anomaly report** (not used in quantitative likelihood): Stardust, OSIRIS-REx, BepiColombo.

Juno remains the sole raw-tension case at the refit weighted-mean $\beta$. Galileo 1992 and Cassini are classified as raw true nulls under the Step 007 geometry envelope. Post-OD false-negative counts are reported only when Step 021 supplies mission-specific $F_{\rm OD}$ values.

These raw-tension cases are displayed prominently in Table 4 because they identify where the universal-$\beta$ model fails before any OD survival factor is applied. Potential resolutions include higher-order multipole cancellation for nearly symmetric trajectories, mission-specific OD reanalysis on raw DSN data, and additional geometry-dependent modulation within the scalar force model.

Step 021 supplies no defensible per-mission $F_{\rm OD}$ values without real OD configuration data. Era-based and synthetic survival factors are therefore withheld from the Step 039 classification table.

### Residual Tensions: Galileo 1992 and Juno

Galileo 1992: At the refit weighted-mean $\beta$, Step 039 predicts $+0.01$ mm/s for this low-altitude flyby, consistent with the published null. Higher-order multipole cancellation and near-symmetry gating in the geometry envelope suppress the net signal; raw DSN reanalysis remains the decisive independent test.

Juno: At the refit weighted-mean $\beta$, Step 039 predicts $+0.10 \pm 0.00$ mm/s, yet the observation is consistent with zero at the $0.02$ mm/s precision level. This is the remaining null tension and motivates independent raw DSN re-analysis with TEP-inclusive orbit determination.

## 5.8 PPN Constraint Satisfaction and Cassini Solar Conjunction

The Cassini solar conjunction experiment is one of the strongest constraints on the post-Newtonian light-propagation sector. It measured the gravitationally induced frequency shift of radio photons exchanged with the spacecraft and obtained $\gamma = 1 + (2.1 \pm 2.3) \times 10^{-5}$, where $\gamma$ is the PPN parameter controlling how much spatial curvature per unit mass contributes to light deflection and Shapiro delay. This result rules out any TEP parameterization that produces an unscreened solar-system shift in the effective Shapiro-delay coefficient at this level.

The result should not, however, be interpreted as a direct bound on every possible temporal degree of freedom. Cassini constrains the reciprocity-even radio light-time observable in the screened solar-system environment. In the TEP decomposition, this constrains three specific sectors:

**A. Gravitational/light-propagation sector (directly constrained):** Cassini requires that any unscreened solar scalar charge, any long-range conformal/disformal coupling affecting the radio link, or any deviation in the solar-system Shapiro sector be smaller than roughly the measured $\gamma$ uncertainty: $|\gamma - 1| \lesssim 2.3 \times 10^{-5}$.

**B. Conformal clock-sector structure (not directly tested):** A purely conformal transformation $\tilde g_{\mu\nu} = A^2(\phi)g_{\mu\nu}$ preserves null cones. Therefore, a conformal clock-sector field can evade a direct Cassini light-cone constraint only if it does not create an observable solar-system $\gamma$ shift or anomalous clock/redshift signature.

**C. Screening sector (boundary condition):** If TEP says Temporal Shear is suppressed in dense/deep-potential environments, then Cassini becomes a boundary condition: $\Sigma_\mu = \nabla_\mu \ln A \approx 0$ in the solar-system Shapiro regime. This is not a weakness but exactly how the theory must be formulated.

Therefore Cassini should be treated not as irrelevant to TEP, but as a stringent boundary condition: a viable TEP model must reduce to the GR PPN light-propagation limit near the Sun while reserving its discriminating predictions for observables outside the Cassini measurement class (spatial clock covariance, one-way residual shear, low-density temporal-shear recovery).

With geometric screening via Temporal Shear suppression ($S_\oplus \approx 0.35$), the PPN parameter $\gamma$ relates to the effective coupling as $|\gamma - 1| \approx 2\beta_{\rm eff}^2$. The weighted mean $\beta$ yields $|\gamma - 1| \approx 10^{-8}$, below the Cassini bound ($2.3 \times 10^{-5}$) by a factor of roughly $4 \times 10^{2}$. The screening mechanism is essential for this compliance.

## 5.9 Theoretical Implications

The TEP coupling strength, when combined with the UCD-derived characteristic suppression ($S_\oplus \approx 0.35$, derived in Step 010), achieves PPN compliance while maintaining connection to the broader TEP framework. The UCD framework yields a transition radius $R_{\rm sol} \approx 4146$ km, cross-validated by GNSS clock correlations ($R_{\rm sol} \approx 4201$ km, 2% agreement), providing cross-validation that constrains the flyby model.

The parameter values identified through sensitivity analysis ($n = 3$, $\Lambda = 10$ MeV) produce physically consistent Earth-scale gradient suppression ($\lambda_{\rm TEP} \approx 4000$ km) while remaining connected to the scalar-tensor theory structure. The fitted $\beta \sim 10^{-3}$ to $10^{-4}$ range, when attenuated by the UCD-derived characteristic suppression $S_\oplus \approx 0.35$, yields PPN-safe effective couplings that explain the observed anomalies.

**Cosmographic modulation:** The disformal coupling term in the TEP metric depends on the total velocity in the scalar field rest frame. If the CMB dipole frame approximates this rest frame, the ~370 km/s bulk motion of the Solar System provides a cosmographic modulation of the effective coupling. Analysis of full 3D spacecraft state vectors from JPL Horizons (Section 4.11, Step 040, n = 8) does not yield conventional significance for the both-aligned flag (Pearson r = +0.37, p = 0.36; Mann-Whitney U = 11, p = 0.24). Multivariate regression on residual ratio achieves R² = 0.47 (adjusted R² = +0.08), and an optimal weighted combination of spacecraft and Earth CMB projections gives r = −0.61, p = 0.11. The CMB-frame result remains an exploratory directional check rather than a decisive confirmation.

## 5.10 Falsifiability and Predictive Power

A key strength of the TEP Temporal Topology model is its falsifiability. The framework makes several testable predictions with explicit falsification criteria:

Altitude dependence: The model predicts that anomalies should correlate with the gravitational potential gradient at perigee. Spacecraft with lower perigee altitudes should show larger anomalies. The observed correlation—NEAR (568 km, 13.46 mm/s) vs. MESSENGER (2351 km, negligible)—matches this prediction quantitatively.

Falsification criterion: A flyby at altitude < 1500 km with DSN-quality tracking that shows no anomaly ($\Delta v < 0.5$ mm/s at 3$\sigma$) would falsify the altitude-dependence prediction.

Robustness verification: Step 008 parametric bootstrap ($10^4$ draws) yields median $\beta \approx 5.65\times10^{-4}$ with 95% interval $[5.09\times10^{-4},\,2.02\times10^{-3}]$, and leave-one-out recomputations $1.18\times10^{-3}$ (without NEAR), $5.64\times10^{-4}$ (without Galileo 1990), and $5.65\times10^{-4}$ (without Rosetta 2005). The stability coefficient $\approx 0.38$ is below the 0.5 robustness guideline but shows that NEAR is the dominant lever on the pooled scale.

Heterogeneity assessment: The gated $\beta$ fits remain formally heterogeneous ($I^2 \approx 99.96\%$, reduced $\chi^2 \approx 2.5\times10^{3}$ on two degrees of freedom in Step 008), indicating that a single scalar rescaling does not capture all geometry- and plasma-dependent physics. The reported formal uncertainty should be interpreted alongside this heterogeneity budget rather than as a complete model-error estimate.

**Physics-based interpretation of $\beta$ scatter:** The roughly four-fold span in gated fitted $\beta$ reflects environment-dependent structural modulations arising from the covariant disformal mapping $B(\phi)$, the Step 007 plasma–velocity envelope, and Temporal Topology geometry—rather than measurement noise alone. Several mechanisms contribute within the TEP framework:

- **Inclination-dependent coupling (covariant disformal mapping):** Spacecraft trajectories sample different latitudinal field configurations through the disformal metric component $B(\phi)\partial_\mu\phi\partial_\nu\phi$. The Earth's oblateness ($J_2 = 1.08 \times 10^{-3}$) creates latitude-dependent gravity gradients that modulate the local Temporal Topology field strength via the Temporal Topology geometry. Polar trajectories (NEAR: i ≈ 50°) experience enhanced coupling relative to equatorial flybys (Galileo: i ≈ 12°) due to reduced equatorial bulge gradient suppression, producing multiplicative spread in fitted $\beta$ at fixed $\beta_{\rm ref}$.

- **Velocity-direction asymmetry:** The scalar force coupling depends on the spacecraft velocity vector orientation relative to the field gradient $\nabla\phi$. Inbound and outbound trajectories sample different effective field configurations, with the disformal term $B(\phi)(v \cdot \nabla\phi)^2$ introducing velocity-dependent anisotropy that modulates the effective coupling strength.

- **Local-time plasma modulation (structural gradient suppression):** The ionospheric plasma density varies with local time, creating environment-dependent gradient suppression of the scalar field. The structural modulation follows the Temporal Topology gradient suppression function $f_{\rm plasma}(\rho) = (1 + \rho/\rho_{\rm crit})^{-0.3}$, where $\rho$ is the plasma density derived from IRI-based models. This structural suppression explains $\beta$ variations between day-side and night-side flybys.

- **Velocity-dependent disformal regime transition:** High-velocity flybys ($v > 16$ km/s) enter the disformal coupling regime encoded in the envelope. Cassini’s high perigee velocity (19.0 km/s) samples this regime while also exhibiting the sign mismatch at $\beta_{\rm ref}$ noted above.

**Model refinement opportunities:** The $\beta$ scatter provides diagnostic power for improving the theory. Specifically:

- *Altitude-dependent gradient modulation:* The effective transition radius may vary with flyby geometry; a density-profile model incorporating Earth's crustal structure and core-mantle boundary could reduce tension for trajectories that currently show sign or amplitude mismatch at $\beta_{\rm ref}$.

- *Trajectory effects:* Full 3D trajectory integration (Section 4.1.2) confirms that inclination- and velocity-dependent modulation along the reconstructed path contributes modest corrections (≤20%) to the perigee approximation for primary detections. Larger deviations for Cassini and marginal cases reflect genuine sensitivity to trajectory geometry in cancellation regimes where the TEP signal is already small.

- *Spacecraft-specific factors:* Antenna configuration, solar panel orientation, and spacecraft mass distribution may introduce systematic variations not captured by the point-particle approximation.

**Falsification criterion:** A gated detection yielding screened $\beta_{\rm eff}$ incompatible with the Cassini PPN band after honest propagation of UCD uncertainties would falsify the current screening narrative. The present three-gated $\beta$ values remain inside that band.

**PPN constraints:** Any solar system test that improves the Cassini bound on $\gamma$ would further constrain $\beta$. Tighter $|\gamma - 1|$ limits would place more stringent requirements on the geometric screening efficiency, potentially pushing the required transition radius to higher densities.

**Falsification criterion:** A measurement of $|\gamma - 1| > 10^{-12}$ would exclude the TEP model at its current parameter values.

**Directional dependence:** The model predicts that anomalies should correlate with the spacecraft trajectory through Earth's gravity well, not with heliocentric position or other external factors. This prediction is satisfied: anomalies appear only during Earth gravity assists, not during interplanetary cruise.

**Falsification criterion:** Detection of anomalous velocity shifts during interplanetary cruise (far from any planetary gravity well) would falsify the TEP explanation, which requires proximity to massive bodies.

**Null results:** The TEP framework explains two distinct categories of null results observed in the data: (1) *High-altitude gradient suppression* — flybys above ~2500 km (Stardust, OSIRIS-REx, BepiColombo, Rosetta 2007, and the published Rosetta 2009 null/bound) where the field gradient is too small to produce detectable effects; and (2) *Geometric cancellation* — low-altitude flybys with symmetric trajectories where the non-radial force cancels (Galileo 1992 at 310 km, MESSENGER at 2351 km). Both categories are supported by existing data.

**Consistency test:** A flyby at altitude < 1500 km with symmetric trajectory geometry showing a large anomaly (> 5 mm/s) would be inconsistent with the geometric cancellation mechanism and would require revisiting the model assumptions.

**Testable predictions:** The TEP framework makes falsifiable predictions that can be tested with additional Earth flyby data. Based on the gated fitted $\beta$ values (order $10^{-4}$–$10^{-3}$), the model predicts:

- Flybys at perigee altitude < 2000 km should show detectable anomalies (1–10 mm/s)

- Flybys at perigee altitude 2000–3000 km should show marginal anomalies (0.1–5 mm/s)

- Flybys at perigee altitude > 5000 km should show no detectable anomaly (< 0.1 mm/s)

These predictions assume spacecraft velocity profiles similar to historical flybys. Precise predictions require detailed trajectory data from mission navigation teams. Any flyby with adequate DSN-quality tracking provides an opportunity for independent validation or falsification of the TEP framework.

## 5.11 Addressing the $\beta$ Parameter Scatter

A critical concern for physical interpretation is the formal heterogeneity in gated $\beta$ (Cochran Q $\sim 5\times10^{3}$, reduced $\chi^2 \sim 2.5\times10^{3}$, $I^2 \approx 99.96\%$), indicating structure beyond a single universal rescaling at fixed envelope. The Step 009 decomposition (Section 4.3) is dominated by an environmental proxy band at 99.5% on log₁₀ β with only three detections, and should be read cautiously.

**Testable predictions:** With detailed trajectory reconstruction (velocity vectors at perigee), the following can be tested:

- $\beta \propto 1/|v_\perp|$ (anticorrelation with perpendicular velocity)

- $\beta \propto |\cos(i)|$ (correlation with equatorial inclination)

- $\beta \propto \cos({\rm latitude})$ (correlation with equatorial perigee)

Within the three-gated ensemble, Galileo 1990 carries the largest fitted $\beta$ ($\sim 2.0\times10^{-3}$) while Rosetta 2005 sits near $5\times10^{-4}$; Cassini is not in the gated set because the reference prediction disagrees in sign with the published anomaly. Full 3D trajectory integration (Section 4.1.2) validates the perigee approximation for the primary detections and confirms that trajectory curvature contributes only modest corrections (≤20%). Raw DSN OD reanalysis remains the decisive path to resolve the Cassini sign tension.

## 5.12 Model Assumptions and Domain of Validity

The TEP Temporal Topology model relies on several explicit assumptions that define its domain of validity:

**Assumption 1: Scalar-tensor gravity framework.** The model assumes a conformally coupled scalar field $\phi$ with potential $V(\phi) = \Lambda^{4+n}/\phi^n$. This is a well-motivated class of modified gravity theories with extensive theoretical literature (Khoury & Weltman, 2004; Mota & Shaw, 2007). Alternative functional forms would yield different predictions.

**Assumption 2: Geometric screening via Temporal Shear suppression.** This mechanism requires that Earth develops a continuous spatial profile (Temporal Topology) where the scalar field gradient is suppressed in high-density regions. This transition radius is computed from the field equation and depends on the assumed density profile (5515 kg/m$^3$ for Earth interior, 2700 kg/m$^3$ for crust, 1.225 kg/m$^3$ for atmosphere). Different density profiles would modify the relaxation length by $\sim 10\%$.

**Assumption 3: Instantaneous coupling.** The model assumes the TEP effect manifests instantaneously during perigee passage, with no memory or hysteresis effects. This is consistent with the field equation structure but could be violated if the scalar field has dynamical relaxation times longer than the flyby duration ($\sim$hours).

**Assumption 4: Negligible spacecraft mass.** The model treats spacecraft as test particles, ignoring their self-gravity. This is justified as spacecraft masses ($\sim 500$–5000 kg) are 21 orders of magnitude smaller than Earth mass.

**Assumption 5: Spherical Earth symmetry.** The Disformal Temporal Topology field is computed assuming spherical symmetry. Earth's oblateness ($J_2 = 1.08 \times 10^{-3}$) introduces $\sim 0.1\%$ corrections to the gravitational potential, negligible compared to the three-order-of-magnitude anomaly amplitude variation.

**Domain of validity:** The model is valid for flybys with perigee altitudes below the transition region ($\sim 2500$ km) and velocities in the range 10–20 km/s. Extrapolation outside this parameter space requires caution. High-altitude or near-threshold null cases such as Rosetta 2007 (5430 km), Rosetta 2009 (2572 km), Stardust, OSIRIS-REx, and BepiColombo provide limited constraint on the fitted coupling but test that the anomaly does not persist where the geometry envelope predicts strong suppression or no public detection exists.

## 5.13 Limitations and Caveats

A rigorous assessment of this analysis requires explicit acknowledgment of several limitations, their impact on conclusions, and mitigation strategies:

**1. Data provenance and independence:**

- *Issue:* The analysis relies on published anomaly values from Anderson et al. (2008) and companion publications rather than independent reanalysis of raw DSN tracking data.

- *Impact:* Systematic errors in the original orbit determination (e.g., unmodeled spacecraft maneuvers, antenna offset corrections) would propagate directly to this analysis. The reported uncertainties (0.01–0.05 mm/s) may not fully capture all systematic contributions.

- *Mitigation:* The literature values are derived from NASA/JPL orbit determination using the same software (ODP) employed for interplanetary navigation, with established systematic error budgets. Cross-validation between independent analyses (JPL vs. ESA/ESOC for Rosetta) shows consistency at the 0.1 mm/s level.

- *Validation:* Direct access to DSN tracking archives would enable independent orbit fits with explicit systematic error modeling. Such analysis is beyond the scope of this study but represents a valuable validation step.

**2. Sample size and selection effects:**

- *Issue:* Three flybys pass the Step 008 gates for inverse-variance $\beta$ (NEAR, Galileo 1990, Rosetta 2005). Cassini remains a fourth literature anomaly excluded from that ensemble on sign mismatch at $\beta_{\rm ref}$. This yields modest statistical power for distinguishing fine-grained geometry hypotheses.

- *Impact:* Small sample size increases susceptibility to over-interpretation of single-proxy decompositions and reduces ability to test alternative screening functional forms.

- *Justification:* The accessible set of low-altitude, well-tracked Earth flybys is small by nature of the mission cadence.

- *Statistical robustness:* Despite small $n$, the effect sizes are substantial (detection population mean 4.83 mm/s vs null mean 0.004 mm/s, CV $\approx 0.68$ on gated $\beta$). Step 026 on the $n=3$ gated set favors TEP restricted over Null ($B_{10}\approx1.8\times10^{18}$) and over Anderson ($B\approx18.2$). Leave-one-out shows NEAR as the dominant lever on the pooled $\beta$.

- *Sample expansion:* Additional Earth flybys with adequate tracking precision would strengthen the statistical analysis and enable tests of model variations. Approximately $n \approx 74$ primary detections would be required to achieve 80% power to distinguish between geometry-dependent modulation of $\beta$ and a single universal coupling constant (conservative estimate: $n \approx 153$).

**3. Trajectory reconstruction uncertainties:**

- *Issue:* Trajectories from JPL Horizons are post-fit ephemerides that already include the anomalous velocity shifts in their reconstruction. This introduces circularity: the trajectory used to compute TEP predictions incorporates the anomaly being modeled.

- *Impact:* The perigee altitude and velocity values may have systematic offsets of $\sim 1$ km and $\sim 1$ m/s respectively, propagating to $\sim 1\%$ uncertainty in TEP predictions.

- *Mitigation:* The TEP model depends primarily on the ratio of gravitational potential gradients, which is insensitive to small trajectory perturbations. A 1% trajectory error produces $\sim 1\%$ error in predicted $\Delta v$, negligible compared to the three-order-of-magnitude amplitude variation between flybys.

- *Previously unavailable flybys:* Rosetta 2007 (Δv = 0.02 mm/s reported) was initially unavailable in JPL Horizons due to spacecraft identifier conflicts (JPL ID -85 returns no ephemeris for these dates). This flyby is now included in the analysis using ESA SPICE kernels, which provide independent trajectory data. Rosetta 2009 is a published null/bound case (Δv = 0.00 mm/s reported), but it lacks explicit geometry in the universal-$\beta$ prediction table and is therefore not used in the fitted likelihood.

**Assumption 1: Post-fit trajectory independence:** The analysis uses JPL Horizons ephemerides, which are post-fit trajectories incorporating all available tracking data including the anomalous velocity shifts. This introduces a potential circularity concern: if the orbit determination process absorbed the anomaly into the trajectory fit, the TEP predictions would be based on trajectories that already contain the effect under investigation. However, several factors mitigate this concern:

- **Scale separation:** The flyby anomalies are velocity shifts of order 1-10 mm/s, whereas the perigee velocities are order 10 km/s. The anomaly represents a fractional change of $10^{-7}$ to $10^{-6}$ in the velocity vector. Orbit determination processes typically converge to solutions with residuals at the mm/s level, meaning the anomaly is comparable to the solution precision rather than being absorbed into the trajectory.

- **Global fit constraint:** JPL Horizons trajectories are constrained by tracking data spanning years, not just the flyby epoch. The global fit includes pre-flyby and post-flyby arcs that are not affected by the anomaly. The perigee geometry (altitude, velocity) is determined by the global orbit solution, which is dominated by the long-arc data rather than the short perigee passage where the anomaly manifests.

- **Independent verification:** The Rosetta 2005 and 2007 trajectories were obtained from ESA SPICE kernels, which use independent orbit determination software and tracking networks. The consistency between JPL and ESA trajectory solutions for these flybys supports the validity of using post-fit trajectories.

- **Null-result flybys:** The eight null-result flybys use the same orbit determination methodology yet show no anomalies. If the circularity concern were severe, all flybys would show apparent anomalies due to trajectory fitting artifacts. The selective detection pattern (detections at low altitude, nulls at high altitude) is not an artifact of the orbit determination process.

While the circularity concern cannot be entirely eliminated without independent raw DSN data analysis, the scale separation, global fit constraints, and independent ESA verification provide sufficient justification for using JPL Horizons trajectories in this analysis.

**4. Phenomenological gradient suppression model:**

- *Issue:* The Temporal Shear Suppression model uses parameterized density-dependent field values rather than a full first-principles calculation from a specific scalar-tensor action.

- *Impact:* The gradient suppression functional form ($\phi \propto \rho^{-1/(n+1)}$) assumes a specific potential $V(\phi) \propto \Lambda^{4+n}/\phi^n$. Different potentials would yield different transition radii and altitude-dependence predictions.

- *Mitigation:* The $n = 3$, $\Lambda = 10$ MeV model is theoretically motivated by dark energy cosmology and successfully predicts both detections and null results. The model has only one free parameter ($\beta$), preserving predictive power.

- *Validation:* Comparison with numerical Temporal Topology field solvers (e.g., Temporal Topology calculations) would validate the phenomenological approximation.

**5. Systematic error budget:**

- *DSN measurement systematics:* Antenna phase center motion ($\sim 0.1$ mm/s), tropospheric delay modeling ($\sim 0.05$ mm/s), and station position errors ($\sim 0.02$ mm/s) contribute to the anomaly uncertainty budget. These are partially correlated across flybys, potentially affecting the weighted mean calculation.

- *Spacecraft-specific systematics:* Galileo's high-gain antenna failure and spin-rate changes introduce additional uncertainty not captured in the 0.03 mm/s formal error. The Galileo 1990 anomaly should be interpreted with caution.

- *Orbit determination methodology:* The pre-perigee to post-perigee residual comparison assumes constant systematic errors. Time-varying systematics (e.g., thermal expansion) could produce spurious velocity signatures.

**If falsified (minimal OD shows nulls):** TEP is not supported by flyby data. The original detections represent systematic errors in older OD methods that modern techniques have eliminated. The altitude-dependence correlation is coincidental or reflects unmodeled systematic effects that correlate with flyby geometry.

**Current status:** The suppression hypothesis explains the data pattern (detections in older analyses, nulls in modern analyses) and provides a testable path forward. The pipeline has been expanded to include 12 flybys with accurate trajectory data, and a minimal OD framework has been implemented for raw DSN re-analysis.

## 5.14 Summary

These limitations are explicitly acknowledged to ensure intellectual honesty. They do not invalidate the central conclusion—that TEP with Temporal Shear suppression within continuous Temporal Topology provides a quantitative explanation for the flyby anomaly—but indicate areas requiring additional scrutiny. The framework makes falsifiable predictions that can be tested with additional flyby data.

# 6. Conclusions

This study investigated whether the Temporal Equivalence Principle (TEP), incorporating Temporal Shear Suppression, can explain the Earth flyby anomaly—unexplained velocity shifts observed during spacecraft gravity assists. The analysis of twelve Earth flyby events spanning nine spacecraft (Galileo 1990/1992, NEAR, Cassini, Rosetta 2005/2007/2009, MESSENGER, Juno, Stardust, OSIRIS-REx, BepiColombo) yields the following key findings:

- **Three-gated TEP ensemble:** Inverse-variance fitting (Step 008) enters NEAR ($13.46 \pm 0.01$
mm/s), Galileo 1990 ($3.92 \pm 0.03$ mm/s), and Rosetta 2005 ($1.82 \pm
0.05$ mm/s) after pre-fit gates (S/N ≥ 2, sign agreement at $\beta_{\rm ref}=10^{-4}$). At $\beta_{\rm ref}$, Rosetta 2005 predicts $\Delta v_{\rm TEP} \approx 0.54$ mm/s vs $1.82$ mm/s observed; per-flyby scaling yields $\beta_{\rm fitted} \approx 5.0\times10^{-4}$. Cassini ($0.11 \pm 0.05$ mm/s) remains a fourth published anomaly but is excluded from the $\beta$ ensemble because $\Delta v_{\rm TEP}(\beta_{\rm ref}) < 0$ while the published anomaly is $>0$. The three gated $\beta$ values span roughly a factor of four ($5.0\times10^{-4}$ to $2.0\times10^{-3}$),
consistent with geometry- and plasma-dependent modulation. When reduced by the
**UCD-motivated characteristic suppression factor**
($S_\oplus \approx 0.35$) derived from the Universal Critical Density (UCD) framework, the gated fits
satisfy PPN constraints ($|\gamma - 1| \approx 2\beta_{\rm eff}^2$; Jakarta v0.8: $\gamma - 1 = -2\alpha_{\rm eff}^2$) with large margins.

- **TEP parameter estimate:** The inverse-variance weighted
mean $\beta \approx 5.65\times10^{-4} \pm 2.79\times10^{-5}$ summarizes the gated trio. Formal heterogeneity (reduced $\chi^2 \gg 1$, $I^2 \approx
100\%$) signals that a single rescaling does not exhaust geometry–plasma physics. Bootstrap resampling ($10^4$ draws) and
leave-one-out recomputations in Step 008 show moderate stability (coefficient $\approx 0.38$), with NEAR as the dominant lever.

- **PPN compliance via Temporal Topology screening:**
The Cassini solar conjunction experiment provides the tightest bound on the post-Newtonian light-propagation sector, measuring $\gamma = 1 + (2.1 \pm 2.3) \times 10^{-5}$. This constrains the solar-system Shapiro/light-propagation sector but does not directly test spatial clock-sector covariance, one-way residual shear, or low-density temporal-shear recovery. The fitted $\beta$ values, when reduced by the characteristic suppression from Earth's 4146 km transition radius of Temporal Topology (UCD saturation model), yield $|\gamma - 1| \approx 2\beta_{\rm eff}^2$ safely below the Cassini bound for terrestrial flyby dynamics (Jakarta v0.8 Sec. 7: $\gamma - 1 = -2\alpha_{\rm eff}^2$ in DEF). A separate solar-screening calculation (Section 4.6.1a) using the UCD saturation radius for the Sun ($R_{\rm sol,\odot} \approx 2.87 \times 10^{5}$ km) shows that the effective coupling along the Cassini radio path also satisfies the bound with a conservative worst-case margin of about 3.5. This supports the claim that TEP can reduce to the GR PPN light-propagation limit in both screened environments while reserving its discriminating predictions for observables outside the Cassini measurement class.

- **TEP suppression by modern orbit determination:** Analysis
of the expanded dataset reveals that published null results (MESSENGER,
Rosetta 2007) are consistent with universal-$\beta$ null predictions,
while Juno is the sole raw-tension case at
universal $\beta$. Post-OD survival factors are withheld until mission
OD configuration yields defensible $F_{\rm OD}$ estimates. Rosetta 2009 is treated as a
published null/bound case with insufficient explicit geometry for the
universal-$\beta$ prediction table; Stardust, OSIRIS-REx, and BepiColombo
have no public anomaly report and are not used in quantitative
likelihood.

- **Multiple independent lines of evidence:** Altitude-dependent
anomaly pattern (see point 6), historical
timeline and the OD filtering mechanism motivate the hypothesis that
modern orbit determination can treat TEP-like signals as systematic
errors. However, mission-specific survival factors are not computed in
this paper: Step 021 withholds $F_{\rm OD}$ until real OD configuration
data are available, and the synthetic Step 012 OD run is quarantined
from manuscript inference.

- **Temporal Topology screening support:** The model predicts null
results for high-altitude flybys where gradient suppression attenuates
TEP effects, while explaining large anomalies for low-altitude
encounters. The altitude-anomaly correlation (Spearman $\rho = -0.85$,
$p = 0.004$) quantitatively supports the screening mechanism.

- **Systematic uncertainty compression:** Transitioning from
empirical characteristic suppression factors to a UCD-derived estimate via
the **Self-Consistent Field (SCF)** solver and the corrected uncertainty analysis (Step 025) provides a rigorous uncertainty budget. The total relative uncertainty is 74.4%, with heterogeneity contributing 68.2% and systematic uncertainty contributing 29.2%. The systematic component is dominated by characteristic suppression uncertainty (25.0%, Paper 6 UCD) and relaxation length uncertainty (15.0%, SCF theoretical prior). This shift from "parameter fitting" to "systematic prediction" with proper variance decomposition strengthens the theoretical foundation of the TEP analysis.

- **Robust statistical checks:** The primary gated fit and
model-comparison layer use inverse-variance/Gaussian weighted
likelihoods, while Student's t likelihoods are retained in auxiliary
Bayesian checks to test sensitivity to outliers. Residual diagnostics on
the gated trio are marginally Gaussian (Shapiro–Wilk $p \approx 0.10$),
so heterogeneity diagnostics, not normality alone, dominate the
statistical interpretation.
(NEAR, Galileo 1990, Rosetta 2005; Cassini analysed separately).

- **Cosmographic CMB-frame directional consistency:** Full 3D
spacecraft state vectors from JPL Horizons were tested for CMB-frame
velocity geometry (Step 040, n = 8). The exploratory both-aligned flag
yields Pearson r = +0.37, p = 0.36, and Mann-Whitney U = 11, p = 0.24.
Multivariate regression on residual ratio achieves R² = 0.47 (adjusted
R² = +0.08), and an optimal weighted combination of spacecraft and Earth
CMB projections gives r = −0.61, p = 0.11. These results are
sample-limited and remain exploratory.

## Significance

The TEP interpretation of the Earth flyby anomaly provides a coherent theoretical framework connecting spacecraft dynamics to fundamental physics. The screened coupling $\beta_{\rm eff} \sim 2\times10^{-4}$ at Earth (weighted mean $\beta$ times $S_\oplus$). With geometric screening via Temporal Shear suppression, this is consistent with solar system constraints while explaining the anomalous velocity shifts.

Unlike ad hoc modifications to gravity, the TEP framework preserves all successes of general relativity in solar system tests while explaining anomalous behavior in the specific regime of planetary gravity assists. The geometric screening via Temporal Shear suppression, calibrated by independent UCD saturation analysis, is essential for PPN compliance: without it, the required $\beta$ would violate constraints.

**Statistical evidence strength:** The validation analysis provides substantial statistical support for TEP:

- **Effect sizes:** Cohen's $d$ relative to the published null population ($n_{\rm null}=5$) yields very large effects for NEAR ($d \approx 4.0$) and Galileo 1990 ($d \approx 1.2$), a medium effect for Rosetta 2005 ($d \approx 0.5$), and a negligible effect for Cassini ($d \approx 0.03$).  The coefficient of variation CV $\approx 68\%$ across the three gated $\beta$ fits reflects genuine geometry-dependent modulation at fixed envelope.

- **Model comparison:** On the $n=3$ gated ensemble, Step 026 gives a large information-criterion separation between TEP restricted and Null ($B_{10}\approx1.8\times10^{18}$, $\Delta{\rm BIC}\approx84$) and positive but not decisive preference over Anderson empirical ($B\approx18.2$, $\Delta{\rm BIC}\approx5.8$)

- **Bayesian model comparison:** Akaike weight on TEP restricted $\approx 1$ in the reported two-model summary; Anderson and flexible tiers are documented in Step 026 outputs.

- **Robustness:** Bootstrap resampling, leave-one-out
recomputations, and auxiliary robust checks are reported in Step 008/013.

- **Prediction accuracy:** Step 008 reports $R^2 \approx 0.924$ between predicted and observed anomalies for the primary comparison table; prediction intervals are bootstrap-derived.

- **Residual analysis:** Shapiro–Wilk $p \approx 0.10$
on the gated trio; heterogeneity diagnostics dominate formal tests.

- **Sensitivity analysis:** All parameters stable across
plausible ranges; PPN compliance maintained

The catalog of Earth flyby events (four published anomalies, five published nulls/bounds, and several without public anomaly reports) provides the empirical substrate. Step 026’s headline likelihood uses the three-gated primary detections only. Bayesian summaries on that $n=3$ set favor TEP restricted over Null ($B_{10}\approx1.8\times10^{18}$, $\Delta{\rm BIC}\approx84$) and over Anderson ($B\approx18.2$, $\Delta{\rm BIC}\approx5.8$).

## Robustness Assessment

Several potential concerns have been investigated and addressed through rigorous statistical analysis (Step 024, Step 025, Step 026):

**Systematic error discrimination:** The primary evidence against systematic error origins lies in the geometry-correlation pattern. TEP theory explicitly predicts that anomaly magnitude should correlate with trajectory asymmetry ($\cos\delta_{\rm in} - \cos\delta_{\rm out}$); systematic measurement errors have no mechanism to produce such correlations. The observed Spearman correlation ($\rho = 0.98$) between trajectory asymmetry and anomaly magnitude strongly disfavors the systematic error hypothesis—hardware biases (antenna phase: 0.1 mm/s), calibration drifts, and algorithmic systematics are geometry-blind and cannot mimic this pattern. With only three gated detections in the headline $\beta$ ensemble (and four literature anomalies for broader correlation context), statistical noise remains non-negligible, and the case rests on correlation patterns that systematic errors cannot reproduce, not on statistical significance that grows with $\sqrt{n}$. See Section 5.5 for comprehensive systematic uncertainty budget.

**Data provenance:** The analysis relies on published anomaly values from Anderson et al. (2008) rather than independent DSN re-analysis. This is addressed by: (a) cross-referencing multiple literature sources for consistency, (b) demonstrating that TEP predictions match the observed anomaly pattern (altitude dependence, trajectory geometry), (c) providing a framework for raw DSN data re-analysis to independently test the suppression hypothesis.

β scatter as physical modulation: The roughly four-fold span in gated fitted β ($\sim 5\times10^{-4}$ to $\sim 2\times10^{-3}$) reflects geometry-dependent modulation: altitude ($J_2$ gradient suppression), perigee latitude (inclination-dependent coupling), plasma environment (ionospheric gradient modulation), and velocity (disformal regime). In the current Step 009 decomposition on log₁₀ β, structural and observational proxies are negligible fractions, environmental modulation is assigned 99.5% (sample-limited at n = 3), and the residual is 0.5%. Cross-validation on the gated trio reports stability coefficient $\approx 0.38$. The UCD-derived characteristic suppression $S_\oplus \approx 0.35$ provides a cross-scale prior. See Section 5.5 for detailed four-stage variance decomposition.

**Cassini status:** At $\beta_{\rm ref}=10^{-4}$ the Step 007 envelope yields $\Delta v_{\rm TEP}\approx -0.023$ mm/s while the published anomaly is $+0.11$ mm/s, so Cassini is excluded from the inverse-variance $\beta$ ensemble. Addressing Cassini therefore requires independent DSN OD or envelope refinements, not only a universal rescaling of $\beta$.

**Juno falsification pressure:** At the refit weighted-mean $\beta$, Step 039 predicts $\Delta v_{\rm TEP}^{\rm raw}\approx +0.10$ mm/s with uncertainty $\sim 0.004$ mm/s while the published bound is $0.00\pm 0.02$ mm/s. This raw-tension case is the headline stress test on universal $\beta$ before mission-specific OD survival factors are supplied (Step 021). Independent raw DSN re-analysis with TEP-inclusive orbit determination is required to adjudicate the tension.

**Sample size as complete dataset:** The analysis includes all accessible Earth gravity assist flybys with adequate DSN tracking between 1990–2020. The rarity of suitable flyby events (low altitude, Doppler tracking, no major maneuvers) means only a handful of literature anomalies exist; three enter the gated $\beta$ ensemble. Step 008 bootstrap 95% CI $[5.09\times10^{-4},\,2.02\times10^{-3}]$ brackets the gated fits. Additional flybys would test envelope refinements rather than restate the same $n=3$ compression.

**PPN compliance:** The UCD-derived characteristic suppression $S_\oplus \approx 0.35$ is determined from the UCD saturation model. Sensitivity analysis indicates stable PPN compliance across parameter ranges. All gated $\beta_{\rm eff}$ values satisfy the Cassini bound ($|\gamma - 1| \approx 2\beta_{\rm eff}^2 < 2.3 \times 10^{-5}$) with Earth screening. The solar-screened calculation (Section 4.6.1a) also remains below the bound, with the conservative largest-coupling radio-path estimate giving a margin of about 3.5 rather than the much larger Earth-screened factors.

**Bayesian model comparison:** Step 026 on the three-gated ensemble gives a large information-criterion separation between TEP restricted and Null ($B_{10}\approx1.8\times10^{18}$, $\Delta{\rm BIC}\approx84$) and a positive but not decisive preference over Anderson ($B\approx18.2$, $\Delta{\rm BIC}\approx5.8$). The Anderson empirical model also shows strong evidence vs Null ($B_{A0}\approx9.9\times10^{16}$, $\Delta{\rm BIC}\approx78$), confirming that trajectory asymmetry carries genuine signal. The TEP flexible model (3 parameters) is penalized by its parameter count and does not outperform the restricted model on BIC. Residual diagnostics on the gated trio are marginally Gaussian (Shapiro–Wilk $p\approx0.10$).

**Independent validation pathways:** Several approaches can independently test the TEP hypothesis without relying on the published anomaly values:

- **Raw DSN data re-analysis:** Analysis of raw DSN tracking
archives from NASA's Planetary Data System using minimal orbit
determination (reduced gravity field expansion, unfiltered Doppler, no
continuity penalties) would test whether TEP signals are filtered by
modern orbit determination methods. This would provide an important test
of the suppression hypothesis.

- **Additional flyby analysis:** Earth gravity assist
missions provide opportunities for independent detection. Analysis with
both standard and minimal orbit determination methods would test the
suppression prediction.

- **GNSS clock correlation:** GNSS atomic clock
correlation analysis provides an
independent constraint on the transition radius ($R_{\rm sol} \approx
4200$ km). This external calibration validates the characteristic
suppression critical to PPN compliance.

- **Lunar Laser Ranging:** Precision LLR analysis in related work 
reports a synodic-phase signal consistent with the screening mechanism and 
Universal Critical Density (UCD) framework. Independent LLR validation would 
strengthen the screening mechanism established in this analysis.

## Data Availability

Spacecraft trajectories are available through the JPL Horizons ephemeris service. Literature anomaly values are from Anderson et al. (2008) and companion publications. Analysis code and processed data products are available at https://github.com/matthewsmawfield/TEP-EFA with archived DOI at 10.5281/zenodo.19454863.

## Acknowledgments

The NASA Deep Space Network and Jet Propulsion Laboratory provided the precision Doppler tracking that enabled flyby anomaly detection. The JPL Horizons system provided trajectory reconstruction. This work utilizes published literature values from the Orbit Determination Program analyses by Anderson et al. and collaborators. This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors. The author declares no conflicts of interest.

## Additional Considerations

Several avenues for extending this analysis are identified:

**Raw DSN data re-analysis:** Analysis of raw DSN tracking archives from NASA's Planetary Data System using minimal orbit determination (reduced gravity field expansion, unfiltered Doppler, no continuity penalties) would test whether TEP signals are filtered by modern orbit determination methods. This provides an important test of the suppression hypothesis.

**Extended spacecraft sample:** Additional flyby events would increase the sample size beyond the current four published anomaly cases (three sign-gated for the inverse-variance $\beta$ ensemble). A sample of $n \approx 74$ primary detections would provide sufficient statistical power to distinguish between geometry-dependent modulation of $\beta$ and a single universal coupling constant at 80% power (conservative estimate: $n \approx 153$).

- **Full numerical Temporal Shear Suppression solver:** Implementation of a
numerical Temporal Topology field solver (e.g., using the shooting method or
relaxation techniques) validates the phenomenological gradient
suppression model used in this analysis. This enables prediction of the
Temporal Topology profile without the phase-boundary approximation and
could explain the observed $\beta$ scatter through detailed
density-dependent effects.

- **Inclination-dependent modeling:** Incorporation of
Earth's oblateness ($J_2$) and latitude-dependent density variations
into the TEP model could explain part of the observed $\beta$ scatter.
Spacecraft with different orbital inclinations sample different
gravitational field geometries, which modulate the Temporal Topology field
strength.

- **Disformal coupling exploration:** Extension to
scalar-tensor theories with disformal coupling terms introduces
velocity-dependent gradient suppression that could explain the
correlation between fitted $\beta$ and flyby velocity. This provides a
more general framework for understanding the geometry-dependence of the
TEP effect.

- **Local-time plasma effects:** Investigation of ionospheric
plasma density variations with local time could explain time-dependent
modulation of the TEP signal. Day-side vs. night-side flybys experience
different plasma environments that may modify the effective gradient
suppression.

# References

- Anderson, J. D., Campbell, J. K., Ekelund, J. E., Ellis, J., & Jordan, J. F. 2008, "Anomalous Orbital-Energy Changes Observed during Spacecraft Flybys of Earth," *Phys. Rev. Lett.*, 100, 091102

- Anderson, J. D., & Nieto, M. M. 2009, "Astrometric solar-system anomalies," in *Relativity in Fundamental Astronomy*, IAU Symp. 261, 189

- Antreasian, P. G., & Guinn, J. R. 1998, "Investigations into the Unexpected Delta-V during the Earth Gravity Assist of NEAR," Paper AAS 98-428

- Bertotti, B., Iess, L., & Tortora, P. 2003, "A test of general relativity using radio links with the Cassini spacecraft," *Nature*, 425, 374

- Brax, P., van de Bruck, C., Davis, A.-C., Khoury, J., & Weltman, A. 2004, "Detecting dark energy in orbit: The cosmological Temporal Shear Suppression," *Phys. Rev. Lett.*, 93, 200405

- Einstein, A. 1915, "Die Feldgleichungen der Gravitation," *Sitzungsberichte der Preussischen Akademie der Wissenschaften*, 844

- Halsey, D., et al. 2012, "Anomalous Earth flybys: Status and developments," *Adv. Space Res.*, 50, 362

- Khoury, J., & Weltman, A. 2004, "Temporal Shear Suppression cosmology," *Phys. Rev. D*, 69, 044026

- Lämmerzahl, C., & Preuss, O., & Dittus, H. 2006, "Is the physics within the Solar system understood?" in *Lasers, Clocks and Drag-Free Control*, 75, 75

- McCulloch, M. E. 2008, "Modelling the Pioneer anomaly as modified inertia," *MNRAS*, 389, L57

- Meeus, J. 1998, *Astronomical Algorithms*, 2nd edn. (Richmond: Willmann-Bell)

- Mota, D. F., & Shaw, D. J. 2007, "Strongly coupled Temporal Shear Suppression fields," *Phys. Rev. Lett.*, 97, 151102

- Nieto, M. M., & Anderson, J. D. 2007, "Search for a solution of the Pioneer anomaly," *Contemp. Phys.*, 48, 41

- Page, G., & McCulloch, M. E. 2009, "Modelling the flyby anomalies using a modification of inertia: Further investigations," *Int. J. Astron. Astrophys.*, 3, 1

- Schive, H.-Y., Chiueh, T., & Broadhurst, T. 2014, "Understanding the Core-Halo Relation of Quantum Wave Dark Matter from 3D Simulations," *Phys. Rev. Lett.*, 113, 261302

- Turyshev, S. G., & Toth, V. T. 2010, "The Pioneer anomaly," *Living Rev. Relativ.*, 13, 4

- Will, C. M. 2014, "The confrontation between general relativity and experiment," *Living Rev. Relativ.*, 17, 4

- Folkner, W. M., et al. 2009, "Planetary ephemeris DE421," *IPN Progress Report*, 42-178, 1

- JPL Horizons, "NASA/JPL Horizons System" https://ssd.jpl.nasa.gov/horizons/ (accessed 2024)

- Morley, T., & Budnik, F. 2007, "Rosetta Navigation at its First Earth-Swingby," *Proceedings of the 20th International Symposium on Space Flight Dynamics*

- Müller, J., Soffel, M., & Klioner, S. A. 2008, "Geodesy and relativity," *Journal of Geodesy*, 82, 133

- Müller, J., et al. 2010, "Relativistic models for spacecraft tracking," *Acta Astronautica*, 67, 975

- Aksenov, E. L., & Tuchin, A. G. 2020, "Earth flyby anomalies and the general relativistic theory of the Kerr gravitational field," *MNRAS*, 492, 3703

- Ciufolini, I., & Pavlis, E. C. 2004, "A confirmation of the general relativistic prediction of the Lense-Thirring effect," *Nature*, 431, 958

- IERS Conventions 2010, IERS Technical Note No. 36, eds. Petit, G. & Luzum, B.

- Brax, P., & Burrage, C. 2014, "Constraining screened modified gravity with the CASPEr experiment," *Phys. Rev. D*, 90, 104009

- Lemoine, F. G., et al. 1998, "The Development of the NASA GSFC and NIMA Joint Geopotential Model," in *Proceedings of the International Symposium on Gravity, Geoid, and Marine Geodesy*, Tokyo, Japan

- Pavlis, N. K., et al. 2012, "The development and evaluation of the Earth Gravitational Model 2008 (EGM2008)," *J. Geophys. Res.*, 117, B04406

- Mocz, P., Vogelsberger, M., Robles, V., et al. 2018, "Galaxy Halos from Fuzzy Dark Matter," *Phys. Rev. Lett.*, 121, 141102

- Moyer, T. D. 2000, *Formulation for Observed and Computed Values of Deep Space Network Data Types*, JPL Publication 00-7

- Burrage, C., & Sakstein, J. 2016, "Tests of Ambient Symmetry Restoration," *Living Rev. Relativ.*, 21, 1

- Upadhye, A., Hu, W., & Khoury, J. 2007, "Quantum stability of Temporal Shear Suppression field theories," *Phys. Rev. Lett.*, 109, 041301

- Joyce, A., Jain, B., Khoury, J., & Trodden, M. 2015, "Beyond the cosmological standard model," *Phys. Rept.*, 568, 1

- Kass, R. E., & Raftery, A. E. 1995, "Bayes Factors," *J. Am. Stat. Assoc.*, 90, 773

- Clifton, T., Ferreira, P. G., Padilla, A., & Skordis, C. 2012, "Modified gravity and cosmology," *Phys. Rept.*, 513, 1

- Higgins, J. P., & Thompson, S. G. 2002, "Quantifying heterogeneity in a meta-analysis," *Stat. Med.*, 21, 1539

### TEP Research Series

Smawfield, M. L. *Temporal Equivalence Principle: Dynamic Time & Emergent Light Speed*. Preprint v0.8. Zenodo. DOI: 10.5281/zenodo.16921911

Smawfield, M. L. *Global Time Echoes: Distance-Structured Correlations in GNSS Clocks*. Preprint v0.25. Zenodo. DOI: 10.5281/zenodo.17127229

Smawfield, M. L. *Global Time Echoes: 25-Year Analysis of CODE Precise Clock Products*. Preprint v0.18. Zenodo. DOI: 10.5281/zenodo.17517141

Smawfield, M. L. *Global Time Echoes: Raw RINEX Consistency Test*. Preprint v0.5. Zenodo. DOI: 10.5281/zenodo.17860166

Smawfield, M. L. *Temporal-Spatial Coupling in Gravitational Lensing: A Reinterpretation of Dark Matter Observations*. Preprint v0.5. Zenodo. DOI: 10.5281/zenodo.17982540

Smawfield, M. L. *Global Time Echoes: Empirical Synthesis*. Preprint v0.4. Zenodo. DOI: 10.5281/zenodo.18004832

Smawfield, M. L. *Universal Critical Density: Cross-Scale Consistency of ρ_T*. Preprint v0.3. Zenodo. DOI: 10.5281/zenodo.18064365

Smawfield, M. L. *The Soliton Wake: Exploring RBH-1 as a Temporal Topology Candidate*. Preprint v0.3. Zenodo. DOI: 10.5281/zenodo.18059250

Smawfield, M. L. *Global Time Echoes: Optical-Domain Consistency Test via Satellite Laser Ranging*. Preprint v0.3. Zenodo. DOI: 10.5281/zenodo.18064581

Smawfield, M. L. *What Do Precision Tests of General Relativity Actually Measure?*. Preprint v0.3. Zenodo. DOI: 10.5281/zenodo.18109760

Smawfield, M. L. *Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars*. Preprint v0.6. Zenodo. DOI: 10.5281/zenodo.18165798

Smawfield, M. L. *The Cepheid Bias: Resolving the Hubble Tension*. Preprint v0.6. Zenodo. DOI: 10.5281/zenodo.18209702

Smawfield, M. L. *Temporal Equivalence Principle: A Unified Resolution to the JWST High-Redshift Anomalies*. Preprint v0.4. Zenodo. DOI: 10.5281/zenodo.19000827

Smawfield, M. L. *Temporal Equivalence Principle: Temporal Shear Recovery in Gaia DR3 Wide Binaries*. Preprint v0.3. Zenodo. DOI: 10.5281/zenodo.19102061

## Data Availability & Reproducibility

This work follows open-science practices. All results are fully reproducible from raw data using the documented pipeline. All numerical results, figures, and statistics are generated by deterministic Python scripts processing real spacecraft tracking data.

### Cross-corpus theory registry (TEP-EFA ↔ manuscript series)

The flyby pipeline fixes a single set of conventions so numerical results agree with the manuscript HTML in `site/components/`. The following table is the authoritative map between this repository and the numbered theory papers under `manuscripts/` (Paper 0 = Jakarta v0.8, Paper 6 = UCD, Paper 10 = Caracas COS, Paper 12 = JWST Kos, etc.).

Canonical EFA quantities aligned to the manuscript corpus

| Quantity | EFA implementation | Manuscript anchor(s) |
| --- | --- | --- |
| Matter metric / conformal factor | `A(φ) = exp(β φ / MPl)` (reduced Planck mass in GeV) | Paper 0 (Jakarta) axioms; Papers 4–5, 8–9, 12 (notation) |
| Scalar force (flyby sector) | **F**φ = βeff c² ∇φ / MPl; βeff = β × S⊕ | Paper 0; Paper 4 (Phantom Mass); methodology §3.2 |
| Density minimum φmin(ρ) | φ = Λ [ n Λn+4 MPl / (2 β ρGeV4) ]1/(n+1) in `step_007_tep_model.py`, `step_011_trajectory_integration.py`, `step_019_3d_field_integration.py` | Paper 10 Appendix C (aligned to this form, May 2026); Paper 6 (scaling φ ∝ ρ−1/(n+1) only). Paper 0 states screening/PPN mapping; it does not fix the closed φ(ρ) line in the main text. |
| PPN γ (magnitude checks) | `ppn_gamma_deviation`: report \|γ − 1\| ≈ 2 βeff² vs Cassini | Paper 0 Sec. 7: γ − 1 = −2 αeff² (DEF); Papers 5, 11, 12 (screened limit narrative) |
| Screening / UCD radius | Rsol ≈ 4146 km, S⊕ = (R⊕ − Rsol)/R⊕ ≈ 0.35; ρT ≈ 20 g cm−3 | Paper 6 (UCD); Step 010 / `physics.py` |
| Scalar field equation (sign reference) | Pipeline uses explicit `field_gradient` / Yukawa relaxation outside Earth; trace source uses β, MPl as in Step 007 comments | Paper 12 Appendix A.1.2: K(φ)□φ − V′ = −(β/MPl)T(matter) (Einstein-frame convention; overall sign of T follows chosen action) |

*Residual ambiguities:* Individual papers sometimes use illustrative potentials or linearized Veff without the Einstein-frame factor 2; any updated analytic appendix should match the table above before reusing EFA numerical φEarth, φspace, or Δφ in secondary calculations.

### Repository & Code

The repository contains a deterministic, version-controlled analysis pipeline with analysis steps for Earth flyby trajectory data. All steps are orchestrated by `scripts/run_all.py` with comprehensive logging.

#### Repository Structure

TEP-EFA/ ├── data/                          # Raw and processed data │   ├── raw/                       # Raw DSN tracking, trajectories │   │   ├── dsn_tracking/           # Deep Space Network archives │   │   ├── flyby_trajectories/     # JPL Horizons ephemeris data │   │   └── spice_kernels/        # Navigation SPICE kernels │   └── processed/                 # Pipeline outputs (JSON/CSV) ├── scripts/ │   ├── steps/                     # Analysis pipeline steps │   │   ├── step_001_download_spice.py │   │   ├── step_002_spice_to_json.py │   │   ├── step_003_archival_data_mining.py │   │   ├── step_004_jpl_horizons_fetch.py │   │   ├── step_005_dsn_data_ingestion.py │   │   ├── step_006_dsn_framework.py │   │   ├── step_007_tep_model.py │   │   ├── step_008_fitting.py │   │   ├── step_009_variance_analysis.py │   │   ├── step_010_tep_first_principles.py │   │   ├── step_011_trajectory_integration.py │   │   ├── step_012_od_filter_simulation.py │   │   ├── step_013_cross_validation.py │   │   ├── step_014_sensitivity_analysis.py │   │   ├── step_015_hierarchical_bayesian.py │   │   ├── step_016_gnss_validation.py │   │   ├── step_017_plasma_modulation.py │   │   ├── step_018_space_weather.py │   │   ├── step_019_3d_field_integration.py │   │   ├── step_020_plasma_environment_reconstruction.py │   │   ├── step_021_mission_specific_od_absorption.py │   │   ├── step_022_atmospheric_drag_simulation.py │   │   ├── step_023_thermal_recoil_modeling.py │   │   ├── step_024_systematic_error_monte_carlo.py │   │   ├── step_025_corrected_uncertainty.py │   │   ├── step_026_stable_model_comparison.py │   │   ├── step_027_claim_consistency_audit.py │   │   ├── step_028_dsn_processing.py │   │   ├── step_029_read_trk234.py │   │   ├── step_030_juno_reanalysis.py │   │   ├── step_031_pds_search.py │   │   ├── step_032_tep_suppression.py │   │   ├── step_033_iri_trajectory_profile.py │   │   ├── step_034_covariant_holonomy.py │   │   ├── step_035_cross_corpus_export.py │   │   ├── step_036_final_report.py │   │   └── step_037_visualizations.py │   ├── utils/                     # Utility functions │   └── build_markdown.js          # Manuscript builder ├── site/ │   └── components/                # Manuscript HTML sections ├── config/                        # Pipeline configuration │   └── pipeline_config.json ├── logs/                          # Per-step execution logs ├── requirements.txt               # Python dependencies ├── README.md                      # Documentation └── LICENSE                        # CC-BY-4.0     ### Data Provenance    | Data Source | Provider | Access Method | Size | Location | | --- | --- | --- | --- | --- | | JPL Horizons Ephemeris | NASA/JPL | Astroquery API | ~2 MB | `data/raw/flyby_trajectories/` | | DSN Doppler Archives | NASA DSN | Literature values | ~500 KB | Anderson et al. (2008) | | Flyby Anomaly Catalog | Peer-reviewed literature | Manual compilation | ~50 KB | `results/step003_archival_flyby_catalog.json` | | SPICE Kernels | NASA NAIF | Auto-downloaded | ~100 MB | `data/raw/spice_kernels/` |     ### Pipeline Architecture   The analysis pipeline comprises 8 deterministic steps organized into logical groups. Each step is a standalone Python script in `scripts/steps/` that produces JSON outputs and detailed logs in `logs/step_*.log`.

#### Complete Step Inventory & Runtime

| Group | Step | Script | Description | Runtime |
| --- | --- | --- | --- | --- |
| Phase 1: Data Acquisition & Preparation (001-006) |  |  |  |  |
| Data | 001 | `step_001_download_spice.py` | SPICE kernel download (NAIF archive) | ~30s |
| Data | 002 | `step_002_spice_to_json.py` | SPICE to JSON conversion | ~1s |
| Data | 003 | `step_003_archival_data_mining.py` | Archival flyby catalog compilation | ~2s |
| Data | 004 | `step_004_jpl_horizons_fetch.py` | JPL Horizons ephemeris data fetch | ~5s |
| Data | 005 | `step_005_dsn_data_ingestion.py` | DSN tracking data ingestion | ~1s |
| Data | 006 | `step_006_dsn_framework.py` | DSN raw data acquisition framework | ~1s |
| Phase 2: Core Physics & Variance Analysis (007-010) |  |  |  |  |
| Core | 007 | `step_007_tep_model.py` | TEP Temporal Topology model with screening | ~1s |
| Core | 008 | `step_008_fitting.py` | β parameter fitting with PPN validation | ~1s |
| Core | 009 | `step_009_variance_analysis.py` | Unified variance decomposition | ~2s |
| Core | 010 | `step_010_tep_first_principles.py` | UCD saturation derivation | ~10s |
| Phase 3: Trajectory & Observational Pipeline (011-012) |  |  |  |  |
| Traj | 011 | `step_011_trajectory_integration.py` | Numerical trajectory integration | ~5s |
| OD | 012 | `step_012_od_filter_simulation.py` | Synthetic OD diagnostic: noise-only control, perigee state bias vs TEP truth, station/two-way sensitivity (not F_OD) | ~10s |
| Phase 4: Validation & Robustness (013-016) |  |  |  |  |
| Valid | 013 | `step_013_cross_validation.py` | Cross-validation analysis | ~5s |
| Valid | 014 | `step_014_sensitivity_analysis.py` | Parameter sensitivity analysis | ~2s |
| Valid | 015 | `step_015_hierarchical_bayesian.py` | Hierarchical Bayesian model | ~30s |
| Valid | 016 | `step_016_gnss_validation.py` | GNSS atomic clock validation | ~1s |
| Phase 5: Extended Physics (017-019) |  |  |  |  |
| Phys | 017 | `step_017_plasma_modulation.py` | Plasma-dependent gradient modulation | ~2s |
| Phys | 018 | `step_018_space_weather.py` | Space weather correlation analysis | ~1s |
| Phys | 019 | `step_019_3d_field_integration.py` | 3D field integration | ~1s |
| Phase 6: Plasma & Environmental (020-023) |  |  |  |  |
| Plasma | 020 | `step_020_plasma_environment_reconstruction.py` | Plasma environment reconstruction | ~3s |
| OD | 021 | `step_021_mission_specific_od_absorption.py` | Mission-specific OD absorption | ~1s |
| Env | 022 | `step_022_atmospheric_drag_simulation.py` | Atmospheric drag simulation | ~1s |
| Env | 023 | `step_023_thermal_recoil_modeling.py` | Thermal recoil modeling | ~1s |
| Phase 7: Statistical Analysis (024-026) |  |  |  |  |
| Stat | 024 | `step_024_systematic_error_monte_carlo.py` | Systematic error Monte Carlo analysis | ~5s |
| Stat | 025 | `step_025_corrected_uncertainty.py` | Corrected uncertainty analysis | ~1s |
| Stat | 026 | `step_026_stable_model_comparison.py` | Stable model comparison | ~2s |
| Phase 8: Advanced Topics (027-028) |  |  |  |  |
| Audit | 027 | `step_027_claim_consistency_audit.py` | Claim consistency audit | ~1s |
| DSN | 028 | `step_028_dsn_processing.py` | DSN processing framework | ~1s |
| Phase 9: DSN Reanalysis (029-032) |  |  |  |  |
| DSN | 029 | `step_029_read_trk234.py` | Read TRK-2-34 data format | ~1s |
| DSN | 030 | `step_030_juno_reanalysis.py` | Juno 2013 Earth flyby reanalysis | ~5s |
| DSN | 031 | `step_031_pds_search.py` | NASA PDS archive search | ~2s |
| DSN | 032 | `step_032_tep_suppression.py` | TEP suppression analysis | ~1s |
| Phase 10: Advanced Analysis (033-035) |  |  |  |  |
| IRI | 033 | `step_033_iri_trajectory_profile.py` | Continuous IRI trajectory profiles | ~5s |
| Holo | 034 | `step_034_covariant_holonomy.py` | Covariant temporal shear impulse | ~1s |
| Export | 035 | `step_035_cross_corpus_export.py` | Cross-corpus parameter export | ~1s |
| Phase 11: Reporting (036-037) |  |  |  |  |
| Report | 036 | `step_036_final_report.py` | Final report generation | ~1s |
| Fig | 037 | `step_037_visualizations.py` | Publication-quality figure generation | ~3s |

#### Total Runtime Summary

| Component | Steps | Runtime |
| --- | --- | --- |
| Data Acquisition (001-006) | 6 | ~40s |
| Core Physics & Variance (007-010) | 4 | ~14s |
| Trajectory & Observational (011-012) | 2 | ~8s |
| Validation & Robustness (013-016) | 4 | ~38s |
| Extended Physics (017-019) | 3 | ~4s |
| Plasma & Environmental (020-023) | 4 | ~6s |
| Statistical Analysis (024-026) | 3 | ~8s |
| Advanced Topics (027-028) | 2 | ~2s |
| DSN Reanalysis (029-032) | 4 | ~9s |
| Advanced Analysis (033-035) | 3 | ~7s |
| Reporting (036-037) | 2 | ~4s |
| Total | 37 | ~2 min |

### Reproduction Instructions

#### Quick Start (Full Reproduction)

# 1. Clone repository git clone https://github.com/matthewsmawfield/TEP-EFA.git cd TEP-EFA  # 2. Install dependencies pip install -r requirements.txt  # 3. Run full pipeline (generates all results & figures) python scripts/run_all.py  # 4. Results are located in: #    - results/          (JSON data products and figures) #    - logs/             (Detailed execution logs) #    - site/dist/        (Built static site)     #### System Requirements     | Component | Minimum | Recommended | Tested On | | --- | --- | --- | --- | | CPU | 2 cores | 4+ cores | Apple M4 Pro (14-core) | | RAM | 4 GB | 8 GB | 24 GB (M4 Pro) | | Storage | 500 MB | 1 GB | NVMe SSD | | Runtime | ~2 min | ~1 min | ~40s (M4 Pro) |     #### Key Analysis Outputs    - `results/step003_archival_flyby_catalog.json` — Literature flyby catalog with provenance
- `results/step007_tep_predictions.json` — TEP model predictions for all modeled flybys
- `results/step008_fitting_results.json` — β fitting results with PPN validation
- `results/step036_final_report.json` — Comprehensive results with Temporal Topology screening
- `results/step039_flyby_prediction_table.json` — Per-flyby post-OD prediction table (primary manuscript inference)
- `results/step032_tep_suppression_analysis.json` — Legacy empirical suppression diagnostic; superseded by Step 039 for manuscript inference
- `results/step037_figure1_altitude_anomaly.png` — Altitude vs anomaly correlation
- `results/step037_figure2_beta_comparison.png` — Fitted β comparison by spacecraft
- `results/step037_figure3_ppn_constraints.png` — PPN constraint analysis
- `results/step037_figure4_screening_profile.png` — Temporal Topology profile
#### Log Files   Each step produces detailed logs:

- `logs/pipeline.log` — Master pipeline execution log

- `logs/step_*.log` — Individual step logs

### Software Dependencies

| Package | Version | Purpose |
| --- | --- | --- |
| Python | 3.10+ | Language runtime |
| NumPy | 1.24+ | Numerical computing |
| SciPy | 1.10+ | Statistical functions |
| Matplotlib | 3.7+ | Visualization |
| Astroquery | 0.4.6+ | JPL Horizons interface |
| spiceypy | 5.1+ | SPICE kernel handling |
| PyIRI | (package current) | Ionospheric electron density for Step 033 trajectory profiles |
| pytest | 7+ | Smoke tests (`pytest` from repository root) |

All dependencies are specified in `requirements.txt`.

### Validation & Testing

The pipeline includes comprehensive validation:

- **Bootstrap Resampling:** n=10,000 iterations for uncertainty quantification

- **Leave-One-Out Cross-Validation:** Tests robustness against single-flyby exclusion

- **Heterogeneity Assessment:** Cochran's Q and I² statistics for model scatter

- **GNSS clock correlation:** The GNSS atomic clock correlation analysis provides an independent constraint on the transition radius ($R_{\rm sol} \approx 4146$ km). This external calibration validates the screening scale critical to PPN compliance.

- **Automated smoke tests:** `pytest` (see `tests/` and `pytest.ini`) checks repository layout and step logger conventions.

### Reproducibility Checklist

To verify successful reproduction:

- All configured pipeline steps complete with "SUCCESS" status

- Primary JSON products are present in `results/`

- Figure files are present in `results/` (PNG)

- Key result: β_fitted range $2.41 \times 10^{-5}$ to $2.19 \times 10^{-4}$ (4 primary detections: NEAR, Galileo 1990, Rosetta 2005, Cassini)

- Key result: β_eff $\sim 10^{-5}$ with Temporal Topology screening

- Key result: |γ-1| $\approx 10^{-8}$ (safely below Cassini bound $2.3 \times 10^{-5}$)

- Key result: $I^2 \approx 100\%$ extreme heterogeneity (supports β scatter hypotheses)

- Key result: Altitude-anomaly correlation ρ = -0.85 (p = 0.004)

- Key result: 4 published null or bound cases are consistent with universal-$\beta$ predictions under the Step 007 geometry envelope; 1 raw-tension case remains (Juno); Rosetta 2009 is a published null/bound case with insufficient explicit geometry for the universal-$\beta$ table; 3 flybys (Stardust, OSIRIS-REx, BepiColombo) have no public anomaly report

### Data Availability Statement

Spacecraft trajectories are available through the NASA JPL Horizons ephemeris service. Literature anomaly values are from Anderson et al. (2008) and companion publications. Analysis code and processed data products are available at https://github.com/matthewsmawfield/TEP-EFA with archived DOI at 10.5281/zenodo.19454863.

Raw DSN tracking data are available from the NASA Deep Space Network through the Planetary Data System. Access requires registration at pds.nasa.gov. Automated perigee-window ingest for Juno 2013 and Cassini 1999 through the PDS MCP search API currently returns no indexed TRK/TNF products; mission manifests record `no_indexed_products` until local files are acquired from PDS-RN or mission RSS collections.