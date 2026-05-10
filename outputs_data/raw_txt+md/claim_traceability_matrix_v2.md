# Claim Traceability Matrix v2 (paper_unified_v3.tex - Public-Data Layer)

> Generated: 2026-05-05T07:11:39.692897+00:00

This matrix links manuscript claim blocks to public-data artifacts and scripts.
It is NOT a replacement for the physics validation ledgers (HEOM/Redfield/Sobol),
but constitutes the L4 (External validation) layer of the traceability hierarchy.

---

## Section Section2.1.1 - Spectral density calibration and uncertainty

| Claim in manuscript | Empirical constraint | Source | Output |
|---|---|---|---|
| eta in [0.1, 1.0] (generic protein range, paper central=0.3) | eta_proxy=0.6026 from B_median=48.21 A^2 (n=362 structures); paper range COVERS proxy range [0.4801,0.7464]: validated=True | `analysis/bath_params_empirical.csv` | `analysis/structures_compact.csv` |
| omega_c in [100, 250] cm^-1 (generic protein range) | No tubulin-specific Raman/THz wavenumber data in public corpus (PubChem PUG-View does not contain cm^-1 spectra; BMRB not yet fetched). Paper range retained with explicit limitation note. | BMRB (not yet fetched) | OPEN - priority future fetch |
| Structural heterogeneity index H_s | H_s = stdev(B_Wilson)/mean(B_Wilson) = 0.5292 (n=362); indicates broad B-factor distribution -> validates eta upper bound | `analysis/bath_params_empirical.csv` | Computed in `compute_detectability_metrics.py` |
| Resolution distribution (n=495, median=2.50 A) | Methods: X-ray=N/A, cryo-EM=N/A, NMR=N/A | `analysis/structures_compact.csv` | RCSB core 507 entries |
| Proxy note | eta_proxy=0.603 from B_median=48.2Ang2 vs FMO_ref B=12Ang2 eta=0.15. Paper range [0.1,1.0] covers [0.480,0.746]. | | |

---

## Section Section3.4 - Quantum coherence utility U_phys per mechanism

| Mechanism | tau_coh (central) | K_req | U_phys@K_est | U>=1 |
|---|---|---|---|---|
| equilibrium | tau_coh=3.910e-14 s | K_req=6.394e+11 | U_phys=1.564e-11 | U>=1: False |
| frohlich | tau_coh=5.000e-11 s | K_req=5.000e+08 | U_phys=0.2 | U>=1: False |
| qed_cavity | tau_coh=5.000e-07 s | K_req=5.000e+04 | U_phys=2 | U>=1: True |
| subradiance | tau_coh=1.000e+01 s | K_req=2.500e-03 | U_phys=4e+05 | U>=1: True |


---

## Section Section6.3 - Fisher information barrier (Experiment 4)

| Claim in manuscript | Value | Source |
|---|---|---|
| Detection probability surface P(doublet | Deltalambda, SNR) | 36-point grid (Deltalambda in [0.3, 1.6]nm x SNR in [10, 10000]) | `analysis/fisher_barrier.csv` |
| SNR_50pct at Deltalambda=1.6nm | ~8 (Fisher CRB estimate) | `fisher_barrier.csv` row Deltalambda=1.60 |
| SNR_95pct at Deltalambda=1.6nm | ~24 | same |
| SNR_50pct at Deltalambda=0.89nm (num. estimate) | ~10 | `fisher_barrier.csv` row Deltalambda=0.89 |
| BIC decisive threshold (paper Section6.3) | SNR ~ 260 (Monte Carlo, 20 resolution elements) | `paper Section6.3` (computed analytically) |

---

## Provenance and reproducibility audit

| Dataset | Rows | Script |
|---|---|---|
| structures_compact.csv | 507 | curate_compact.py |
| studies_compact.csv | 95 | curate_compact.py |
| spectral_compact.csv | 892 | curate_compact.py |
| metrics_compact.csv | 892 | compute_detectability_metrics.py |
| bath_params_empirical.csv | 1 | compute_detectability_metrics.py |
| fisher_barrier.csv | 36 | compute_detectability_metrics.py |
| comparative_panels_compact.csv | 4 | run_comparative_panels.py |
| data_registry.csv | 395 | build_registry.py |

---

## Open gaps (honest accounting for reviewers)

1. **BMRB NMR chemical shifts**: Not yet fetched. Would provide residue-level
   linewidth data for omegac constraint. Estimated 1-2 day fetch + curation sprint.
2. **Tubulin-specific Raman/THz wavenumbers**: PubChem PUG-View does not store
    spectroscopic spectra in machine-readable JSON for small-molecule modulators.
   The specific protein vibrational spectrum of alphabeta-tubulin requires either:
   (a) EuropePMC full-text parsing of Raman papers (Gascoyne 2011, Craddock 2017), or
   (b) MD simulation of 1JFF dimer (priority future work, acknowledged in Section2.1.1).
3. **eta_proxy is a linear rescaling estimate, not a measured value**: It
   empirically *validates* that the paper range [0.1, 1.0] is not arbitrary,
   but does not *replace* the MD-derived spectral density. This distinction
   must be preserved in any paper text referencing these data.
