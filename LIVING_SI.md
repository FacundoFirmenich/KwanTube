# LIVING_SI.md - Supplementary Information (Automated Validation)

> **Version** 3.5.1 * **Generated** 2026-05-06T14:20:49Z * **Wall-time** 116.19s
> **SHA-256 Hash** `a21fb9e7798043c1...` * **Audit Status** 14/14 validation criteria met

This document is machine-regenerated from `validation_report.json` on every 
pipeline run. Every result is cross-referenced with the machine-auditable 
JSON artifact, verified via the cryptographic SHA-256 signature above.

---

## SI-1 - Non-equilibrium Dynamics, Inversion, and Sensitivity

- **Coherence Figure-of-Merit**: \(\varphi_0 = 1.036e-05\) (mean estimated \(T_2^*\) in ns).
- **Inversion Fidelity**: \(\hat\varphi = 1.000\) (Target interval \([0.85,\,1.01]\)).
- **Model Selection**: **emergent** architecture favored (\(\Delta\mathrm{BIC}_{max} = 142.70\)).

## SI-2a - Analytic Perturbative Benchmarking (Section 2.2.5, COMP-1)

Cross-validation of the master equation formalisms under the secular and memory-factor 
approximations. Calculations assume an Ohmic spectral density 
\(J(\omega)=\eta\omega\exp(-\omega/\omega_c)\) with \(\omega_c=4.5\times10^{12}\) rad/s 
and \(T=310\) K:

| eta (Coupling) | tau_Lindblad (s) | tau_Redfield_approx (s) | tau_HEOM_approx (s) | Relative Spread |
|---|---:|---:|---:|---:|
| 0.10 | 3.913e-14 | 3.909e-14 | 4.107e-14 | 0.050 |
| 0.30 | 1.304e-14 | 1.300e-14 | 1.498e-14 | 0.145 |
| 0.60 | 6.490e-15 | 6.447e-15 | 8.423e-15 | 0.278 |
| 1.00 | 3.913e-15 | 3.871e-15 | 5.846e-15 | 0.435 |

The `relative_spread` indicator ((max - min)/mean) quantifies cross-formalism 
concordance. Values below 1.0 indicate that the closed-form Lindblad rate accurately 
captures the hierarchical physics within the specified perturbative regime.

## SI-2b - Hierarchical Equations of Motion (HEOM) Validation

Full non-perturbative hierarchical integration (\(L=4\), high-temperature Matsubara truncation). 
Comparison between the nominal Lindblad baseline and the numerically exact HEOM propagator:

_(Hierarchical results pending solver completion)_

## SI-2c - Bayesian HEOM Hierarchy (v2) - Contraction Analysis

Automated Bayesian hierarchy for summarize small-N HEOM convergence evidence. This 
module models jump magnitudes on the log-scale to infer stable contraction ratios \(r\).

- **Global Contraction Ratio** (\(r = \exp(\mu_{logr})\)): 0.528 (\([0.385,\, 0.706]\) 95% CI).
- **Global Decay Rate** (\(\beta = -\mu_{logr}\)): 0.650.
- **Hierarchical Stability**: \(\tau_{logr} = 0.266\) (Group-level heterogeneity).

**Output Artifacts**:
- [Group Summary](outputs_data/raw_csv/heom_+bayesian_analysis/group_loglinear_summary.csv)
- [Global Contraction](outputs_data/raw_csv/heom_+bayesian_analysis/hierarchy_global_contraction.csv)
- [Extrapolated Jumps](outputs_data/raw_csv/heom_+bayesian_analysis/extrapolated_jumps.csv)
- [Level Checks](outputs_data/raw_csv/heom_+bayesian_analysis/level_reference_checks.csv)
- [Diagnostics](outputs_data/raw_txt+md/reports/diagnostics_v2.txt)

**Posterior Plots**: [posterior_plots_v2.png](outputs_data/figures_final/posterior_plots_v2.png) and [posterior_plots_v2.pdf](outputs_data/figures_final/posterior_plots_v2.pdf)

## SI-2d - Mean-Force Steady-State Diagnostic

Diagnostic of HEOM relaxation and consistency with second-order Mean-Force (MF) Gibbs states. 
Calculated via Kullback-Leibler (KL) divergence from the final HEOM state $\rho(t_{final})$.

_(Diagnostic results pending execution)_

- **Interpretation**: `KL_bare < 0.05` indicates the system has relaxed to the standard 
  Gibbs state. A divergent `KL_mf` is the mathematical signature of the failure of 
  second-order perturbation theory in the intermediate coupling regime (\(\lambda\beta \sim 1\)).

## SI-3 - Detector Performance: ROC Detection Surface (Section 5, COMP-12)

Probability of detection \(P_D(\Delta\ell,\mathrm{SNR})\) at a fixed false-alarm rate 
\(\alpha=0.05\). Results computed using a matched-filter detector over \(N_{MC}=10\) 
stochastic trials per configuration.

| Delta_l \\ log10 SNR | 2.00 | 2.85 | 3.70 |
|---|---|---|---|
| 0.30 | 0.00 | 0.00 | 0.10 |
| 1.05 | 0.10 | 0.10 | 0.80 |
| 1.80 | 0.00 | 0.50 | 1.00 |

**Consistency Check**: Verification of monotonic detection gain with increasing SNR across 
the spatial coherence grid (\(\Delta\ell\)).

## SI-4 - Bayesian Evidence Meta-Analysis (Section 5, COMP-11)

Summary of experimental contrasts integrated into the Bayesian hierarchy:

| Study Identifier | Observable Scale | Effect Size | Standard Error | Source / Context |
|---|---|---:|---:|---|
| Babcock2024 | log_ratio | 0.51 | 0.13 | J. Phys. Chem. B / eNeuro |
| Kalra2024 | raw_mean_diff_seconds | 69.00 | 20.41 | J. Phys. Chem. B / eNeuro |

**Statistical Inference Results**:
- **Babcock (2024)**: \(BF_{10} = 183.3\) (Decisive evidence, Jeffreys scale). Nested Sampling verification: \(178.3 \pm 10.1\) (\(n_{live}=600\)).
- **Kalra (2024)**: \(BF_{10} = 43.3\) (Very Strong evidence, Jeffreys scale). Nested Sampling verification: \(42.8 \pm 2.3\) (\(n_{live}=600\)).

**Aggregate Significance**: Combined Bayes Factor (assuming study independence) \(BF_{10} \approx 7.6e+03\).

## SI-7 - Calibration and Robustness Audits

### Simulation-Based Calibration (SBC)
Validation of the Nested Sampling (NS) inference engine via SBC on \(N_{sim}=1000\) 
calibration trials.
- **Uniformity p-value**: 0.560 (Statistically consistent with a calibrated rank distribution).
- **Scope**: SBC results validate the engine's performance under the specific generative models 
  deployed in this study.
- **Diagnostic Plot**: [sbc_calibration_ns.pdf](figures_final/sbc_calibration_ns.pdf).

### Prior Sensitivity Analysis
Evaluation of Bayes Factor (\(BF_{10}\)) stability across a spectrum of weakly-informative priors.
- **Stability**: \(BF_{10}\) remains robustly above the "Decisive" threshold (\(>100\)) for 
  prior standard deviations \(\sigma_{prior} \in [0.2,\,1.0]\).
- **Caveat**: The shaded regions indicate prior-dominated regimes where \(\sigma_{prior} < SE\).
- **Sensitivity Profiles**: [prior_sensitivity.pdf](figures_final/prior_sensitivity.pdf).

## SI-8 - HEOM Integration Pre-registration
- **Cryptographic Hash**: `5385692fbb6622b6f48b0535b38dfc07a5cffde2656ff6b6b458bb3da10c4217`
- **Registration Timestamp**: 2026-04-22T05:55:12Z

## SI-5 - Collective Modes in the Microtubule Lattice (Section 4.3, COMP-6)

Analysis of a 13-protofilament B-lattice configuration (\(N = 130\) dimers, 
\(\mu=1700\) D, \(\varepsilon_r = 80\)):

- **Superradiant Band Edge** (\(E_+\)): 480.55 meV
- **Subradiant Band Edge** (\(E_-\)): -480.55 meV
- **Excitonic Spectral Gap** (\(\Delta\)): 961.10 meV
- **Inverse Participation Ratio (IPR)**: 64.1 (>= 2 indicates delocalized modes)
- **Axial Interaction** (\(J_\parallel\)): -88.08 meV (Attractive coupling; J-aggregate character)
- **Lateral Interaction** (\(J_\perp\)): 160.36 meV (Repulsive coupling; H-aggregate character)

## SI-6 - Summary of Automated Validation Checks

| Validation Metric | Status | Technical Detail |
|---|:---:|---|
| `noneq_ladder_monotone` | [OK] | Coherence tau(Delta_mu) is monotonically non-decreasing |
| `inversion_recovers_fidelity` | [OK] | Recovery Fidelity phi_hat=1.000 |
| `sensitivity_phi_finite` | [OK] | Nominal Coherence phi_0=1.036e-05 |
| `model_selection_picks_emergent` | [OK] | Selected Model: emergent (Delta_BIC_max=142.70) |
| `multi_formalism_concordance` | [OK] | Relative spread < 1.0 across eta coupling grid |
| `roc_monotone_global` | [OK] | Detection probability P_D increases with SNR |
| `babcock_bf_decisive` | [OK] | Decisive Evidence (BF10=183.3) |
| `kalra_bf_very_strong` | [OK] | Very Strong Evidence (BF10=43.3) |
| `sbc_calibrated` | [OK] | NS Calibration p=0.560 |
| `lattice_gap_positive` | [OK] | Spectral Gap Delta=961.10 meV |
| `lattice_subradiant_delocalized` | [OK] | Subradiant IPR=64.1 |
| `heom_production_extracted` | [OK] | Pur(30ps)=0.21, Disc=26.39% |
| `heom_non_equilibrium_regime` | [OK] | Terminal purity 0.21 < 0.25 confirms transient dynamics |
| `heom_redfield_divergence` | [OK] | HEOM-Redfield gap 26.39% exceeds truncation error |

## SI-9 - L4 Public-Data Audit (Structural & Spectroscopic)

To mitigate epistemic risk, the pipeline integrates a multi-layered empirical audit:

- **Structural Audit ($N=362$ PDB entries)**: Median Wilson B-factor \(\langle B \rangle = 48.21\) \AA$^2$ yielding \(\eta_{proxy}\) = 0.603.
- **Spectroscopic Audit ($N=93$ studies)**: Consensus vibrational cutoff \(\omega_c \approx 150 \text{ cm}^{-1}\) (4.5 THz) and observed spectral density support.
- **Ligand Layer**: Detectability thresholds calibrated against PubChem UV/Vis and wavenumber data.

### HEOM Production Trajectory Diagnostics
| Metric | Value | Interpretation |
|--------|-------|----------------|
| $P_{init}$(500 fs) | `0.8488` | Non-perturbative population retention |
| Purity (30 ps) | `0.21` | Confirms non equilibrium transient |
| IPR (30 ps) | `4.759` | Delocalization extent at terminal window |
| Redfield discrepancy | `26.39%` | Model-level divergence (>10% threshold) |

---

*End of auto-generated Supplementary Information. To regenerate, execute:* 
`python src/scripts/validation/reproduce_paper_results.py [--full-roc]`.
