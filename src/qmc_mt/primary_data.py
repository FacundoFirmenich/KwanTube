"""
Primary-literature evidence registry for microtubule quantum-coherence
claims. Each entry documents EXACTLY what the primary paper reports and
what is extractable, with page/table provenance.

NO pooling of heterogeneous scales. NO imputation of missing SE/N.
Studies flagged `usable_for_logeffect=False` are NOT passed to any
meta-analytic estimator.
"""
from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Optional
import math


@dataclass(frozen=True)
class StudyRecord:
    key: str
    citation: str
    doi_or_url: str
    claim_domain: str               # "optical" | "behavioral" | "electronic" | "docking"
    scale: str                      # "log_ratio" | "raw_mean_diff_seconds" | "none"
    effect: Optional[float]         # point estimate in the stated scale
    se: Optional[float]
    n: Optional[int]
    usable_for_logeffect: bool
    provenance: str                 # table/figure/page pointer
    notes: str = ""


BABCOCK_2024 = StudyRecord(
    key="Babcock2024",
    citation=("Babcock, Kurian et al. (2024). Ultraviolet Superradiance "
              "from Mega-Networks of Tryptophan in Biological Architectures. "
              "J. Phys. Chem. B 128, 4035-4046."),
    doi_or_url="10.1021/acs.jpcb.3c07936",
    claim_domain="optical",
    scale="log_ratio",
    effect=math.log(17.6 / 10.6),                          # 0.50700...
    se=math.sqrt((2.1/17.6)**2 + (0.6/10.6)**2),           # 0.13207...
    n=None,                                                 # N of replicates not tabulated
    usable_for_logeffect=True,
    provenance=("Table 1, p. 4040 (PDF p. 5): QY-Trp @ 280 nm, "
                "MT = 17.6 +/- 2.1, TuD = 10.6 +/- 0.6. Delta method applied."),
    notes=("Effect = ln(QY_MT / QY_TuD). SE via first-order delta method "
           "on independent Gaussian errors of numerator and denominator. "
           "Likely positive instrumental correlation would reduce SE and "
           "increase BF; current value is therefore a conservative lower bound."),
)


KALRA_2024 = StudyRecord(
    key="Kalra2024",
    citation=("Khan, Kalra, Hameroff et al. (2024). Microtubule-Stabilizer "
              "Epothilone B Delays Anesthetic-Induced Unconsciousness in Rats. "
              "eNeuro 11(8), ENEURO.0291-24.2024."),
    doi_or_url="10.1523/ENEURO.0291-24.2024",
    claim_domain="behavioral",
    scale="raw_mean_diff_seconds",
    effect=69.0,
    se=max(69.0 - 37.4, 109.0 - 69.0) / 1.96,               # 20.408...
    n=8,
    usable_for_logeffect=False,                             # no pre/post means
    provenance=("p. 7 main text (mean diff = 69 s, p = 0.0016, N = 8); "
                "Table 2, p. 23/29 (95% CI [37.4, 109] s); "
                "Figure 1 legend, p. 22 (SE across sessions)."),
    notes=("95% CI [37.4, 109] is asymmetric around mean=69 (likely bootstrap). "
           "SE derived from the larger half-width as a conservative Gaussian "
           "approximation. Using 1.96 assumes normality; the asymmetric CI "
           "suggests a skewed distribution where Gaussian SE is conservative."),
)


BANDYOPADHYAY_2013 = StudyRecord(
    key="Sahu2013",
    citation=("Sahu, Ghosh, Hirata, Fujita, Bandyopadhyay (2013). Atomic "
              "water channel controlling remarkable properties of a single "
              "brain microtubule. Biosensors and Bioelectronics 47, 141-148."),
    doi_or_url="10.1016/j.bios.2013.02.050",
    claim_domain="electronic",
    scale="none",
    effect=None, se=None, n=None,
    usable_for_logeffect=False,
    provenance=("Section2.6 and Fig. 5, pp. 145-146. Fig. 5f lists dominant "
                "peaks {12, 20, 22, 30, 101, 113, 185, 204} MHz as "
                "'statistically most-occurred', selected from 'thousands "
                "of noise and actual peaks'."),
    notes=("No mean +/- SE and no observational N tabulated for a "
           "resonant/non-resonant conductance contrast. Mechanistic "
           "evidence only. Does NOT enter any pooled estimator."),
)


CRADDOCK_2012 = StudyRecord(
    key="Craddock2012",
    citation=("Craddock, St. George, Freedman, Barakat, Damaraju, Hameroff, "
              "Tuszynski (2012). Computational Predictions of Volatile "
              "Anesthetic Interactions with the Microtubule Cytoskeleton. "
              "PLoS ONE 7(6): e37251."),
    doi_or_url="10.1371/journal.pone.0037251",
    claim_domain="docking",
    scale="none",
    effect=None, se=None, n=None,
    usable_for_logeffect=False,
    provenance=("Table 2, p. 4 (per-site halothane binding energies, e.g. "
                "2.54, 2.70, 2.91, 3.39, 3.31 kcal/mol with persistence %); "
                "p. 6 (9/47 sites persisted >70% of 5-ns MD; 5/47 all of it; "
                "estimated Kd 6-16 mM)."),
    notes=("Computational docking + MD, not an observational series. No "
           "mean +/- SE and no N for a binding-shift contrast. Mechanistic "
           "evidence only. Does NOT enter any pooled estimator."),
)


REGISTRY: tuple[StudyRecord, ...] = (
    BABCOCK_2024, KALRA_2024, BANDYOPADHYAY_2013, CRADDOCK_2012,
)


def usable_logeffect_studies() -> list[StudyRecord]:
    return [s for s in REGISTRY if s.usable_for_logeffect]


def provenance_table() -> list[dict]:
    return [asdict(s) for s in REGISTRY]


if __name__ == "__main__":
    import json
    print(json.dumps({
        "registry": provenance_table(),
        "usable_for_logeffect": [s.key for s in usable_logeffect_studies()],
        "n_usable": len(usable_logeffect_studies()),
    }, indent=2, default=str))
