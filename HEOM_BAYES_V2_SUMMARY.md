# HEOM Bayesian Hierarchy v2 — Contraction Summary (current ledger)

This summary is generated from the canonical v2 outputs in:

- `heom_bayes_out_v2/hierarchy_global_contraction.csv`
- `heom_bayes_out_v2/hierarchy_group_shrinkage.csv`
- `heom_bayes_out_v2/diagnostics_v2.txt`

## Global contraction

- `r = exp(mu_logr)` posterior mean: **0.526**
- 95% interval: **[0.385, 0.706]**

## Group-level shrunk contraction ratios

| Group | Median r | 95% interval | Note |
|---|---:|---:|---|
| fragment6dpu dCoh | 0.645 | [0.587, 0.708] | slower contraction |
| fragment6dpu dFrob | 0.464 | [0.448, 0.481] | strong contraction |
| fragment6dpu dPop | 0.475 | [0.456, 0.494] | strong contraction |
| full1jff site5 jump | 0.524 | [0.350, 0.794] | 2 points only; ratio calibration |

## Interpretation guardrail

The v2 hierarchy stabilizes contraction inference for small-$N$ groups by modeling positive jump observables on the log scale and pooling `log(r_g)` through random effects.

It does **not** create new HEOM evidence and does **not** replace new full-system HEOM runs.
