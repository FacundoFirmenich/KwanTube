# HEOM Acceptance Criteria — Pre-Registered Specification

**Build target:** v3.7.0 (post-HEOM integration into SI-2)
**Status:** PRE-REGISTERED — written before HEOM solver output is inspected.
**Purpose:** Fix acceptance thresholds for HEOM vs. perturbative baselines
*a priori*, so the verdict on cross-formalism consistency is mechanical
and not adjusted to the observed numbers.

---

## 1. Rationale

The hierarchical equations of motion (HEOM; Tanimura & Kubo 1989;
Ishizaki & Fleming 2009) are a numerically exact solver for open-quantum-
system dynamics with a Drude-Lorentz (overdamped Brownian) bath, in the
limit of sufficient hierarchy depth and Matsubara truncation. Secular
Redfield and Lindblad are perturbative approximations that are known to
fail progressively as the system-bath coupling strength \(\eta\) increases
and as the bath correlation time \(\tau_c\) approaches the system timescale
\(\hbar/\Delta E\).

Literature benchmarks establish the following expectations, which this
document codifies as thresholds:

- Ishizaki & Fleming 2009 (*J. Chem. Phys.* 130, 234111): Redfield
  agrees with HEOM within ~10–15% on transfer timescales in the weak-
  coupling / Markovian regime (\(\eta \lesssim 0.1\), \(\tau_c \ll \hbar/\Delta E\)),
  and diverges by factors of 2–5 at \(\eta \gtrsim 1\) (strong coupling,
  non-Markovian).
- de la Lande et al. 2023 and Kreisbeck & Kramer 2012: at intermediate
  coupling (\(\eta \sim 0.3\)), perturbative methods typically capture
  qualitative trends (ordering of timescales, existence of oscillations)
  but quantitative errors of 20–40% on \(\tau_{\text{tr}}\) and
  \(\tau_{\text{coh}}\) are normal.

The criteria below encode these literature-established tolerances. A
failure of Redfield to match HEOM at \(\eta=1.0\) is **not** a failure
of the build — it is an **expected physical result** that Redfield is
outside its regime of validity. What would be a failure is HEOM itself
violating physical constraints (positivity, trace preservation,
thermalization).

---

## 2. Definitions

Let:

- \(\tau_{\text{tr}}^{X}(\eta)\) = transfer time (1/e decay of donor
  population to its long-time value) computed by method \(X \in
  \{\text{HEOM}, \text{Redfield}, \text{Lindblad}\}\) at coupling \(\eta\).
- \(\tau_{\text{coh}}^{X}(\eta)\) = coherence lifetime (1/e decay of
  \(|\rho_{12}(t)|\) from its transient peak) by method \(X\) at coupling
  \(\eta\).
- Relative deviation:
  $$
  \Delta_X(\eta) \;\equiv\; \frac{\tau^{X}(\eta) - \tau^{\text{HEOM}}(\eta)}{\tau^{\text{HEOM}}(\eta)}.
  $$
- Relative spread across all three methods:
  $$
  S(\eta) \;\equiv\; \frac{\max_X \tau^{X}(\eta) - \min_X \tau^{X}(\eta)}{\text{median}_X \tau^{X}(\eta)}.
  $$

HEOM is the reference (\(\Delta_{\text{HEOM}} \equiv 0\) by construction).

---

## 3. Acceptance thresholds — Regime I: weak coupling (\(\eta = 0.1\))

Redfield is expected to be quantitatively accurate here.

| Quantity | Threshold | Rationale |
|---|---|---|
| \(\lvert \Delta_{\text{Redfield}}(\tau_{\text{tr}}) \rvert\) | \(< 0.15\) | Ishizaki-Fleming benchmark |
| \(\lvert \Delta_{\text{Redfield}}(\tau_{\text{coh}}) \rvert\) | \(< 0.20\) | Coherence more sensitive to secular approximation |
| \(\lvert \Delta_{\text{Lindblad}}(\tau_{\text{tr}}) \rvert\) | \(< 0.25\) | Lindblad loses Lamb shift |
| \(S(\eta=0.1)\) | \(< 0.30\) | Global consistency |

**Verdict rule:** all four must pass → `regime_weak_passed = True`.

---

## 4. Acceptance thresholds — Regime II: intermediate coupling (\(\eta = 0.3\))

Redfield is expected to deviate moderately. Qualitative agreement
(same ordering of \(\tau_{\text{tr}}^{\text{1JFF}} > \tau_{\text{tr}}^{\text{6DPU}}\))
is required; quantitative agreement is relaxed.

| Quantity | Threshold | Rationale |
|---|---|---|
| \(\lvert \Delta_{\text{Redfield}}(\tau_{\text{tr}}) \rvert\) | \(< 0.40\) | Intermediate regime tolerance |
| \(S(\eta=0.3)\) | \(< 0.60\) | Spread allowed to grow |
| Sign of \(\tau_{\text{tr}}^{\text{1JFF}} - \tau_{\text{tr}}^{\text{6DPU}}\) under HEOM | must be **positive** | Same ordering as Redfield result reported in build v3.6.0 (894 fs vs. 83 fs) |

**Verdict rule:** all three must pass → `regime_intermediate_passed = True`.

The ordering check is the critical one: if HEOM inverts the ordering,
the entire "6DPU as enhanced regime" narrative collapses and must be
reframed or retracted.

---

## 5. Acceptance thresholds — Regime III: strong coupling (\(\eta = 1.0\))

Redfield is **outside its regime of validity**. No quantitative agreement
with HEOM is required or expected. What is required is that HEOM itself
behaves as a physical density-matrix evolution.

| Quantity | Threshold | Rationale |
|---|---|---|
| \(\min_t \lambda_{\min}(\rho^{\text{HEOM}}(t))\) | \(> -10^{-6}\) | Positivity (numerical) |
| \(\lvert \text{Tr}(\rho^{\text{HEOM}}(t)) - 1 \rvert_{\max}\) | \(< 10^{-4}\) | Trace preservation |
| \(\lVert \rho^{\text{HEOM}}(t \to \infty) - \rho^{\text{Gibbs}} \rVert_{\text{F}}\) | \(< 0.05\) | Thermalization to Boltzmann |
| Hierarchy convergence: \(\lVert \rho^{\text{HEOM}}_{N_{\max}} - \rho^{\text{HEOM}}_{N_{\max}-1} \rVert_{\text{F}}\) at \(t = t_{\text{final}}\) | \(< 10^{-3}\) | Depth truncation adequate |
| Matsubara convergence: same norm, varying \(K\) | \(< 10^{-3}\) | Bath expansion converged |

**Verdict rule:** all five must pass → `regime_strong_passed = True`.

If `regime_strong_passed = False`, the HEOM solver has a numerical
problem and the \(\eta=1\) row of SI-2 must be flagged as "solver did
not converge" rather than reported as a physical result.

---

## 6. Global verdict logic

```
heom_integration_accepted = (
    regime_weak_passed
    AND regime_intermediate_passed
    AND regime_strong_passed
)
```

- If `heom_integration_accepted = True`: SI-2b is promoted from
  supplementary to main-text cross-validation. The narrative
  "Redfield is validated against numerically exact HEOM in the regime
  where it applies" is permitted.
- If only Regime I and/or II pass: SI-2b is kept supplementary. The
  narrative must read "HEOM confirms Redfield in weak-to-intermediate
  coupling; strong-coupling regime requires HEOM and is reported as
  such." No retraction needed.
- If Regime I fails: **stop**. Either Redfield or HEOM has a bug.
  Do not publish v3.7.0. Debug before proceeding.

---

## 7. What would count as a *negative* result worth publishing

If HEOM reveals that:

1. At \(\eta=0.3\), \(\tau_{\text{tr}}^{\text{6DPU,HEOM}} / \tau_{\text{tr}}^{\text{1JFF,HEOM}} > 0.5\)
   (the one-order-of-magnitude gap collapses), or
2. At \(\eta=1.0\), \(\tau_{\text{coh}}^{\text{6DPU,HEOM}} < 500\) fs
   (the claimed "prolonged quantum stability" disappears in the
   non-perturbative regime),

then the central physical claim of the manuscript is weakened. These
are **not** build failures — they are scientific findings. The
manuscript must be rewritten to report them honestly rather than
suppressed. A `narrative_revision_required` flag shall be raised in
`LIVING_SI.md` and SI-2b reported as the primary result, with Redfield
demoted to a baseline comparison.

---

## 8. Pre-registration hash

This document will be hashed (SHA-256) and the hash written into
`LIVING_SI.md` **before** the HEOM solver output is read. Any
post-hoc modification of the thresholds above invalidates the
pre-registration and must be disclosed as such in the manuscript.

```
pre_registration_sha256: 5385692fbb6622b6f48b0535b38dfc07a5cffde2656ff6b6b458bb3da10c4217
commit_timestamp:        2026-04-22T05:55:12Z
heom_output_timestamp:   <must be strictly later than commit_timestamp>
```

---

## 9. Non-adjustable parameters

The following cannot be retuned after seeing HEOM output:

- Coupling grid: \(\eta \in \{0.1, 0.3, 1.0\}\)
- Temperature: \(T = 310\) K (physiological)
- Bath cutoff: \(\omega_c = 100\) cm\(^{-1}\) (Drude)
- Hierarchy depth target: \(N_{\max}\) such that convergence criterion
  in §5 is met (chosen by solver, not by author)
- Matsubara terms target: \(K\) such that convergence criterion in §5
  is met

Any deviation requires a new pre-registration with its own hash.

---

*End of pre-registered specification.*
