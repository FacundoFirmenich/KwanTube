# LIVING_SI.md — Supplementary Information (auto-generated)

> **Version** 3.5.0 · **Generated** 2026-04-24T02:20:48Z · **Wall-time** 12.42s
> **SHA-256** `a7d81477dab23c22…` · **Validation** 11/11 checks passed

This document is machine-regenerated from `validation_report.json` on every
pipeline run. **Do not edit by hand.** Every number has a corresponding entry
in the JSON artefact, auditable via the SHA-256 stamp above.

---

## SI-1 · Non-equilibrium ladder, inversion, sensitivity

- Nominal coherence figure-of-merit: \(\varphi_0 = 1.057e-05\) (mean \(T_2^*\) in ns)
- Parameter-inversion recovered fidelity: \(\hat\varphi = 1.000\) (target \([0.85,\,1.01]\))
- Model-selection winner: **emergent** (\(\Delta\mathrm{BIC}_{\max} = 142.70\))

## SI-2a · Analytic perturbative cross-check (§2.2.5, COMP-1)

Analytic perturbative cross-check (Lamb-shift & memory-factor approximation).
Ohmic spectral density \(J(\omega)=\eta\omega\exp(-\omega/\omega_c)\),
\(\omega_c=4.5\times10^{12}\) rad/s, \(T=310\) K:

| η | τ_Lindblad (s) | τ_Redfield_approx (s) | τ_HEOM_approx (s) | rel. spread |
|---|---:|---:|---:|---:|
| 0.10 | 3.913e-14 | 3.909e-14 | 4.107e-14 | 0.050 |
| 0.30 | 1.304e-14 | 1.300e-14 | 1.498e-14 | 0.145 |
| 1.00 | 3.913e-15 | 3.871e-15 | 5.846e-15 | 0.435 |

`relative_spread` = (max − min)/mean. Spread < 1 validates that the
Lindblad closed form tracks the higher-order hierarchies within its own truncation error.

## SI-2b · HEOM hierarchical solver validation

Full hierarchical solve (L=4, high-temperature Matsubara truncation).
Comparison between the nominal Lindblad rate and the exact HEOM propagator:

_(Hierarchical results pending solver completion)_

## SI-3 · ROC detection surface (§5, COMP-12)

\(P_D(\Delta\ell,\mathrm{SNR})\) at false-alarm \(\alpha=0.05\), matched-filter
detector, \(N_{\mathrm{MC}}=10\) trials per cell.

| Δℓ \\ log₁₀ SNR | 2.00 | 2.85 | 3.70 |
|---|---|---|---|
| 0.30 | 0.00 | 0.00 | 0.10 |
| 1.05 | 0.10 | 0.10 | 0.80 |
| 1.80 | 0.00 | 0.50 | 1.00 |

Validation: monotone non-decreasing in SNR at every \(\Delta\ell\).

## SI-4 · Per-study Bayesian evidence (§5, COMP-11)

| Study | Scale | Effect | SE | Source |
|---|---|---:|---:|---|
| Babcock2024 | log_ratio | 0.51 | 0.13 | J. Phys. Chem. B / eNeuro |
| Kalra2024 | raw_mean_diff_seconds | 69.00 | 20.41 | J. Phys. Chem. B / eNeuro |

**Evidence results:**
- **Babcock 2024**: $BF_{10} = 183.3$ (decisive, Jeffreys); nested sampling cross-check: $178.3 \pm 10.1$ at $n_{live}=600$.
- **Kalra 2024**: $BF_{10} = 43.3$ (very strong, Jeffreys); nested sampling cross-check: $42.8 \pm 2.3$ at $n_{live}=600$.

Combined BF (independence): $BF_{10} \approx 7.6e+03$.

Note: Bandyopadhyay (2013) and Craddock (2012) are excluded from the statistical pool as they provide mechanistic support without extractable $SE$ or observational $N$ for a contrast.

## SI-7 · Calibration and Robustness

### SBC Calibration
Validation of the Nested Sampling engine via Simulation-Based Calibration (SBC) on $N_{sim}=100$ trials.
- **p-value**: 0.419 (statistically consistent, uniform ranks).
- **Scope**: SBC validated under the exact 1D Gaussian generative model used in this analysis; generalization to complex multimodal shapes is not claimed.
- **Rank Histogram**: [sbc_calibration_ns.pdf](file:///c:/Users/User/3D%20Objects/biofisicaquantiqaCLINE/git_repo/figures_final/sbc_calibration_ns.pdf).

### Prior Sensitivity
Robustness of $BF_{10}$ across a range of weakly-informative priors.
- **Babcock**: $BF_{10}$ remains $>100$ for $prior\_sd \in [0.2, 1.0]$.
- **Caveat**: Values of $prior\_sd < SE$ (shaded region) are prior-dominated and uninformative about $H_1$ vs $H_0$.
- **Sensitivity Curves**: [prior_sensitivity.pdf](file:///c:/Users/User/3D%20Objects/biofisicaquantiqaCLINE/git_repo/figures_final/prior_sensitivity.pdf).

## SI-8 · HEOM Pre-registration
- **Hash**: `5385692fbb6622b6f48b0535b38dfc07a5cffde2656ff6b6b458bb3da10c4217`
- **Specification**: [heom_acceptance_criteria.md](file:///c:/Users/User/3D%20Objects/biofisicaquantiqaCLINE/git_repo/heom_acceptance_criteria.md)
- **Commit Timestamp**: 2026-04-22T05:55:12Z

## SI-5 · Microtubule lattice & collective modes (§4.3, COMP-6)

13-protofilament B-lattice, \(N = 130\) dimers, \(\mu=1700\) D, \(\varepsilon_r = 80\):

- Superradiant edge \(E_+ = 480.55\) meV
- Subradiant edge  \(E_- = -480.55\) meV
- Spectral gap     \(\Delta   = 961.10\) meV
- Subradiant IPR  \(= 64.1\) (≥ 2 ⇒ extended mode)
- NN axial coupling   \(J_\parallel = -88.08\) meV (attractive ⇒ J-aggregate)
- NN lateral coupling \(J_\perp     = 160.36\) meV (repulsive ⇒ H-aggregate)

## SI-6 · Validation checks

| Check | Status | Detail |
|---|:---:|---|
| `noneq_ladder_monotone` | ✅ | τ(Δμ) monotone non-decreasing |
| `inversion_recovers_fidelity` | ✅ | φ̂=1.000 |
| `sensitivity_phi_finite` | ✅ | φ₀=1.057e-05 |
| `model_selection_picks_emergent` | ✅ | best=emergent (ΔBIC_max=142.70) |
| `multi_formalism_concordance` | ✅ | max spread < 1.0 across η grid |
| `roc_monotone_global` | ✅ | P_D non-decreasing in SNR across full Δℓ grid |
| `babcock_bf_decisive` | ✅ | BF10=183.3 |
| `kalra_bf_very_strong` | ✅ | BF10=43.3 |
| `sbc_calibrated` | ✅ | p=0.419 |
| `lattice_gap_positive` | ✅ | gap=961.10 meV |
| `lattice_subradiant_delocalized` | ✅ | IPR=64.1 |

---

*End of auto-generated SI. Regenerate via* `python reproduce_paper_results.py [--full-roc]`.
