# HEOM Bayesian hierarchy — modeled observables and ready-to-run CSV

## Scope
This file registers the **observable groups that can already be modeled immediately**
with `bayesian_heom_hierarchy.py` using only results that are explicitly available
from the current chat session.

It deliberately excludes single-point diagnostics that do not yet have at least
two distinct `NC` levels, because the current Bayesian script requires at least
two `NC` values per group.

## Modeling principle

The hierarchy script fits, for each group \(g\),

\[
y_{g,n} \sim \mathcal{N}\left(\theta_g + \alpha_g e^{-\beta_g (n-n_{0,g})},\ \sigma_g\right),
\]

with posterior summaries for:

- asymptotic limit \(\theta_g\)
- contraction speed \(\beta_g\)
- convergence ratio \(r_g = e^{-\beta_g}\)
- residual numerical scale \(\sigma_g\)

## Important semantic rule for “jump” observables

For all `jump_*` groups, the `nc` column is the **upper NC level of the pair**:

- `nc = 4` means the jump from `NC=3 -> 4`
- `nc = 5` means the jump from `NC=4 -> 5`
- etc.

This is intentional. These jump sequences are modeled as quantities that should
contract toward **zero** as `NC` increases.

Therefore, for all `jump_*` groups, the `reference` column is set to `0.0`.

## Groups included in the ready CSV

### Tier A — strong fragment convergence groups
These have 5 NC levels each and are the most informative groups for immediate
hierarchical fitting.

1. `fragment6dpu_jump_dPop_max@500fs`
   - values: 5.51e-02, 2.90e-02, 1.29e-02, 5.45e-03, 2.09e-03
   - reference: 0.0

2. `fragment6dpu_jump_dCoh_max@500fs`
   - values: 3.57e-03, 2.87e-03, 1.43e-03, 9.10e-04, 5.06e-04
   - reference: 0.0

3. `fragment6dpu_jump_dFrob_max@500fs`
   - values: 9.97e-02, 4.91e-02, 2.24e-02, 9.45e-03, 3.76e-03
   - reference: 0.0

### Tier B — fragment level observable with Richardson anchor
4. `fragment6dpu_site0_population_level@500fs`
   - NC=7: 0.486401
   - NC=8: 0.486348
   - external asymptotic reference from Richardson: 0.486315

This group is weaker than Tier A because it only has two NC levels, but it is
still admissible for the current Bayesian script and connects directly to the
Richardson analysis.

### Tier C — full-system 1JFF jump observable
5. `full1jff_jump_dPop_site5@500fs`
   - NC=5: 9.15e-03  (jump 4->5)
   - NC=6: 5.02e-03  (jump 5->6)
   - reference: 0.0

This is also only a two-point group, so inference will be weaker, but it is
useful as a first hierarchical anchor for the full 1JFF system.

## Diagnostics known but NOT included in the ready CSV
These are real results, but they are currently **not eligible** for the existing
hierarchy script because they lack enough NC levels per group.

### Single-pair ultrashort 1JFF probe at 40 fs
- `dPop max = 2.88e-08`
- `dCoh max = 1.10e-05`
- `dFrob max = 1.92e-05`

Reason for exclusion:
- only one pair (`NC=5 -> 6`) is available.

### HEOM vs Redfield 1JFF site-5 level at ~500 fs
- Redfield: 0.6248
- HEOM NC=6: 0.8516
- absolute discrepancy: 2.27e-01
- relative discrepancy: 26.63%

Reason for exclusion:
- only one HEOM NC level is explicitly available for this level observable.

### Living SI summary-level statistics
- Sobol/Saltelli `S1(eta) = 1.0216 [0.935, 1.113]`
- SBC p-value = 0.80
- Bayes factor range = 12.5–266.7
- KL divergence = 0.3271 nats
- full-system convergence ratio `r = 0.549`
- projected `eps7 = 0.612%`
- projected 30 ps cost ≈ 51.2 h

Reason for exclusion:
- these are summary diagnostics, not NC-indexed observable series.

## Ready CSV file
The ready-to-run dataset is:

- `heom_bayes_ready_observables.csv`

## Recommended execution
Run:

```bash
python bayesian_heom_hierarchy.py heom_bayes_ready_observables.csv --output-dir bayes_heom_out
```

## Interpretation guidance

### What should be most informative?
The three Tier A fragment jump groups should drive the strongest posterior signal.

### What should be treated cautiously?
The two-point groups:
- `fragment6dpu_site0_population_level@500fs`
- `full1jff_jump_dPop_site5@500fs`

are useful but should not be overclaimed.

### What would most improve the hierarchy next?
The single best improvement would be to add:
- one more NC level for the full 1JFF jump observable,
- and one more HEOM NC level for a full-system level observable at 500 fs.

That would move the full-system side from “anchored” to “seriously informative”.

## Provenance of the numbers
All values in the CSV come from results explicitly stated in this chat:
- `observable_level_convergence_refined.py`
- `richardson_fragment.py`
- `jff_calibration_check.py`
- `compare_heom_vs_redfield_1jff.py`

No values were invented or interpolated.