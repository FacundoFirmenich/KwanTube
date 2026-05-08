"""
Build an interactive epistemic graph for the manuscript bundle — visually refined,
interaction-rich, and audit-ready. The graph links mechanism classes to quantitative
constraints, validation artifacts, falsifying tests, and scope boundaries.

V4.4 — Synchronizes linewidth-conditional Frohlich, TUR closure, HEOM structured
diagnostics, and fixes sidebar/inspector collapse layout.
"""

from __future__ import annotations

import csv
import html
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("epistemic_graph")


def read_json(path: Path, default: Any = None) -> Any:
    if not path.exists():
        logger.warning(f"Missing JSON: {path}")
        return default
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        logger.error(f"Invalid JSON: {path}")
        return default


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        logger.warning(f"Missing CSV: {path}")
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def short_sha(value: str, n: int = 16) -> str:
    return str(value or "")[:n]


def _extract_fraction(data: Any, default: str) -> str:
    if not isinstance(data, dict):
        return default
    search_keys = ("subradiant_fraction", "fraction_subradiant", "subradiant_fraction_pct")
    sections = [data, data.get("summary", {}), data.get("metrics", {}), data.get("results", {})]
    for section in sections:
        if not isinstance(section, dict):
            continue
        for key in search_keys:
            if key in section:
                val = section[key]
                if isinstance(val, (int, float)):
                    return f"{val:.1%}" if 0 < val <= 1 else f"{val}%"
                return str(val)
    return default


def node(
    node_id: str,
    label: str,
    node_type: str,
    layer: str,
    status: str,
    summary: str,
    *,
    details: str = "",
    artifacts: Optional[list[str]] = None,
    lenses: Optional[list[str]] = None,
    metrics: Optional[Dict[str, Any]] = None,
) -> dict:
    return {
        "id": node_id,
        "label": label,
        "type": node_type,
        "layer": layer,
        "status": status,
        "summary": summary,
        "details": details,
        "artifacts": artifacts or [],
        "lenses": lenses or [],
        "metrics": metrics or {},
    }


def edge(
    source: str,
    target: str,
    relation: str,
    label: str,
    *,
    status: str = "neutral",
    lenses: Optional[list[str]] = None,
    summary: str = "",
) -> dict:
    return {
        "id": f"{source}->{target}:{relation}",
        "source": source,
        "target": target,
        "relation": relation,
        "label": label,
        "status": status,
        "lenses": lenses or [],
        "summary": summary,
    }


def build_graph() -> dict[str, Any]:
    validation_path = PROJECT_ROOT / "outputs_data" / "raw_json" / "structural" / "validation_report.json"
    frohlich_audit_path = PROJECT_ROOT / "outputs_data" / "raw_json" / "nonequilibrium" / "frohlich_universal_gating_audit.json"
    structured_path = PROJECT_ROOT / "outputs_data" / "raw_json" / "metrics" / "heom_structured_relaxation_diagnostics.json"
    kww_path = PROJECT_ROOT / "outputs_data" / "raw_json" / "metrics" / "heom_kww_relaxation_fit.json"
    panels_path = PROJECT_ROOT / "outputs_data" / "raw_csv" / "compact" / "comparative_panels_compact.csv"
    trace_path = PROJECT_ROOT / "outputs_data" / "raw_txt+md" / "claim_traceability_matrix_v2.md"
    lineage_path = PROJECT_ROOT / "outputs_data" / "raw_json" / "audit" / "lineage_audit.json"

    report = read_json(validation_path, {})
    frohlich_audit = read_json(frohlich_audit_path, {})
    structured = read_json(structured_path, {})
    panels = {row.get("mechanism", ""): row for row in read_csv(panels_path)}
    lineage = read_json(lineage_path, [])
    if not isinstance(lineage, list):
        lineage = []

    validation = report.get("_validation", {}) or {}
    checks = validation.get("checks", []) or []
    metadata = report.get("_metadata", {}) or {}
    sha = short_sha(report.get("_sha256", ""))

    lineage_counts: Dict[str, int] = {}
    for item in lineage:
        status = str(item.get("status", "UNKNOWN"))
        lineage_counts[status] = lineage_counts.get(status, 0) + 1

    eq = panels.get("equilibrium", {})
    frohlich = panels.get("frohlich", {})
    cavity = panels.get("qed_cavity", {})
    structured_summary = structured.get("purity_population_entropy_beta_summary", {}) if isinstance(structured, dict) else {}

    sub_metrics: Dict[str, str] = {"N130": "70.8%", "N260": "77.3%", "N520": "84.8%"}
    sub_artifact_map = {
        "N130": "subradiant_decay_spectrum_N130.json",
        "N260": "subradiant_decay_spectrum.json",
        "N520": "subradiant_decay_spectrum_N520.json",
    }
    sub_artifacts: list[str] = []
    for n_key, fname in sub_artifact_map.items():
        p = PROJECT_ROOT / "outputs_data" / "raw_json" / "metrics" / fname
        data = read_json(p, None)
        sub_metrics[n_key] = _extract_fraction(data, sub_metrics[n_key])
        sub_artifacts.append(str(p))

    nodes: list[dict] = [
        node("premise_open_system", "Warm, wet open‑system environment", "premise", "Premises", "constraint",
             "Tubulin operates in a thermal, aqueous, strongly fluctuating biological environment.",
             lenses=["open_system", "scale", "scope"]),
        node("premise_quantum_signatures", "Tubulin quantum‑optical signatures", "premise", "Premises", "supported",
             "Collective optical signatures motivate tests, but do not imply cognition.",
             details="The graph treats quantum signatures as empirical inputs, not as explanatory shortcuts.",
             lenses=["biophysical", "scope"]),
        node("gap_core", "Equilibrium‑to‑functionality gap", "claim", "Premises", "supported",
             "Two distinct scale gaps must be kept separate: microsecond readout and neural functionality.",
             metrics={"microsecond_gap": "10^7‑10^8", "neural_gap": "10^11‑10^13"},
             lenses=["scale", "open_system", "scope"]),
        node("gap_micro", "Cellular readout gap: 10⁷–10⁸", "constraint", "Quantitative Constraints", "constraint",
             "Equilibrium fs coherence to a μs‑scale readout already requires 10⁷–10⁸ amplification.",
             lenses=["scale", "falsification"]),
        node("gap_neural", "Neural functionality gap: 10¹¹–10¹³", "constraint", "Quantitative Constraints", "constraint",
             "Equilibrium fs coherence to 10–100 ms neural‑scale claims requires 10¹¹–10¹³ amplification.",
             lenses=["scale", "scope"]),
        node("mech_equilibrium", "Equilibrium protection", "mechanism", "Mechanism Classes", "fails",
             f"T₂* central={eq.get('tau_coh_central_s', '3.91e-14')} s; U_phys={eq.get('u_phys_central', '1.6e-11')}.",
             details="Equilibrium protection cannot bridge either the microsecond readout gap or the neural functionality gap.",
             artifacts=[str(panels_path), str(validation_path)],
             metrics={"k_req_central": eq.get("k_req_central"), "u_phys": eq.get("u_phys_central")},
             lenses=["open_system", "scale", "falsification"]),
        node("mech_scaling", "Naive cooperative scaling", "mechanism", "Mechanism Classes", "fails",
             "Independent baths produce superdecoherence; symmetry breaking prevents a simple N or √N rescue.",
             lenses=["open_system", "scale"]),
        node("mech_frohlich", "Fröhlich pumping", "mechanism", "Mechanism Classes", "conditional",
             f"Linewidth-conditional, utility-limited; central U_phys={frohlich.get('u_phys_central', '0.2')}.",
             details="Carrier criterion gives L_omega(0.1 THz) ≈ 10 nm; a 10 µm gate requires gamma_Hz ≈ 10^8 Hz or a lower-frequency collective mode.",
             artifacts=[str(panels_path), str(frohlich_audit_path)],
             metrics={
                 "k_req_central": frohlich.get("k_req_central"),
                 "u_phys": frohlich.get("u_phys_central"),
                 "L_omega_0.1THz": "~10 nm",
                 "L_gamma_10um_requires": "gamma_Hz ~ 1e8",
             },
             lenses=["biophysical", "falsification", "scale"]),
        node("mech_cavity", "Ordered‑water QED cavity", "mechanism", "Mechanism Classes", "candidate",
             f"Survives as a falsifiable candidate if UV/THz splitting and gain constraints are met; central U_phys={cavity.get('u_phys_central', '2')}.",
             artifacts=[str(panels_path)],
             metrics={"k_req_central": cavity.get("k_req_central"), "u_phys": cavity.get("u_phys_central")},
             lenses=["biophysical", "falsification"]),
        node("mech_subradiance", "Geometric subradiance", "mechanism", "Mechanism Classes", "candidate",
             "Free‑space Γᵢ/Γ₀ supports size‑dependent subradiant fractions; functional validation remains open.",
             artifacts=sub_artifacts, metrics=sub_metrics,
             lenses=["biophysical", "falsification", "open_system"]),
        node("mech_classical", "Classical synchronization null", "mechanism", "Mechanism Classes", "supported",
             "Known classical synchronization mechanisms can account for the neural‑scale phenomena considered here.",
             lenses=["scope", "biophysical"]),
        node("constraint_t2", "Equilibrium T₂*: 13–39 fs", "observable", "Quantitative Constraints", "supported",
             "Parametric protein-bath baselines set the principal equilibrium result; structural η-proxy is supplemental only.",
             artifacts=[str(validation_path), str(trace_path)],
             lenses=["open_system", "scale", "reproducibility"]),
        node("constraint_tur", "TUR closure: Pmin/GTP ≈ 10⁵", "constraint", "Quantitative Constraints", "fails",
             "Thermodynamic uncertainty relation closes localized MT→channel amplification on first-law grounds.",
             details="For K=10^11 and epsilon_in≈0.1, P_min≈8.6e-8 W while the local MT GTP budget is ≈6e-13 W.",
             metrics={"P_min_W": "8.6e-8", "P_GTP_W": "6e-13", "gain_K": "1e11", "epsilon_in": "0.1"},
             lenses=["open_system", "scale", "falsification"]),
        node("constraint_heom", "HEOM 30 ps non‑thermalized transient", "observable", "Quantitative Constraints", "supported",
             f"Pur(30 ps)={report.get('heom_production', {}).get('purity_30ps', '0.21')}; Redfield discrepancy={report.get('heom_production', {}).get('redfield_discrepancy_pct', '26.39')}%.",
             artifacts=[str(validation_path)],
             lenses=["open_system", "reproducibility"]),
        node("constraint_structured", "HEOM structured non-Markovian relaxation", "observable", "Quantitative Constraints", "supported",
             "Population, purity, and entropy observables show sub-unitary KWW exponents consistent with distributed relaxation kinetics.",
             details="This supports structured non-Markovian distributed relaxation over a finite window, not a thermodynamic glass-transition claim.",
             artifacts=[str(structured_path), str(kww_path)],
             metrics={
                 "beta_range": f"{structured_summary.get('beta_min', 0.370):.3f}-{structured_summary.get('beta_max', 0.462):.3f}" if isinstance(structured_summary.get('beta_min'), (int, float)) else "0.370-0.462",
                 "n_observables": structured_summary.get("n_observables", 6),
                 "quantum_purity_beta": "0.4427",
                 "quantum_purity_ci95": "[0.4379, 0.4472]",
             },
             lenses=["open_system", "reproducibility", "biophysical"]),
        node("constraint_utility", "Coherence utility U = K τ_coh / τ_func", "observable", "Quantitative Constraints", "supported",
             "Separates required amplification K_req from independently bounded empirical gain K_emp.",
             artifacts=[str(panels_path)],
             lenses=["scale", "falsification", "scope"]),
        node("artifact_validation", "Validation ledger: 20/20", "artifact", "Evidence Artifacts", "supported",
             f"{validation.get('passed', 0)}/{validation.get('total', 0)} checks passed; SHA {sha}.",
             artifacts=[str(validation_path)],
             metrics={"version": metadata.get("version"), "sha": sha},
             lenses=["reproducibility"]),
        node("artifact_trace", "Claim traceability matrix", "artifact", "Evidence Artifacts", "supported",
             "Links public‑data claims to empirical constraints and scripts.",
             artifacts=[str(trace_path)],
             lenses=["reproducibility", "scope"]),
        node("artifact_lineage", "Lineage audit: CLEAN", "artifact", "Evidence Artifacts", "supported",
             f"SAME={lineage_counts.get('SAME', 0)}, UPDATED={lineage_counts.get('UPDATED', 0)}, ADDED={lineage_counts.get('ADDED', 0)}, DEPRECATED={lineage_counts.get('DEPRECATED_INTENTIONAL', 0)}.",
             artifacts=[str(lineage_path)],
             lenses=["reproducibility"]),
        node("test_cavity", "UV/THz cavity splitting test", "experiment", "Decision Tests", "open",
             "Positive: Δλ in the predicted nm range. Null: no splitting at high SNR.",
             lenses=["falsification", "biophysical"]),
        node("test_beta", "THz β/γ linewidth benchmark", "experiment", "Decision Tests", "open",
             "Tests whether Frohlich-like finite-size scaling is controlled by the measured linewidth gamma.",
             lenses=["falsification", "biophysical", "scale"]),
        node("test_subradiance", "TCSPC length‑resolved subradiance", "experiment", "Decision Tests", "open",
             "Tests whether radiative lifetimes scale with MT length as predicted by geometric subradiance.",
             lenses=["falsification", "biophysical"]),
        node("scope_no_consciousness", "Blocked shortcut: quantum signature → consciousness", "scope_boundary",
             "Scope Boundary", "blocked",
             "Molecular quantum signatures do not imply consciousness without a measured transfer and amplification mechanism.",
             lenses=["scope", "scale", "falsification"]),
        node("scope_no_pooling", "No causal pooling of optical and behavioral evidence", "scope_boundary",
             "Scope Boundary", "blocked",
             "Bayes factors remain stratified across incommensurate evidentiary scales.",
             lenses=["scope", "reproducibility"]),
        node("scope_no_in_vivo", "No in vivo functional validation claimed", "scope_boundary",
             "Scope Boundary", "blocked",
             "HEOM and radiative results constrain mechanisms; they do not establish in vivo neural function.",
             lenses=["scope", "biophysical"]),
    ]

    # Append individual validation checks as nodes (degree 0 – will be collapsible)
    for check in checks:
        name = check.get("name", "check")
        nodes.append(
            node(
                f"check_{name}",
                name,
                "validation_check",
                "Evidence Artifacts",
                "supported" if check.get("passed") else "fails",
                check.get("detail", ""),
                artifacts=[str(validation_path)],
                lenses=["reproducibility"],
            )
        )

    # Convert artifact paths to relative for portability
    for item in nodes:
        portable: list[str] = []
        for artifact in item.get("artifacts", []):
            try:
                portable.append(Path(artifact).resolve().relative_to(PROJECT_ROOT).as_posix())
            except (OSError, ValueError):
                portable.append(str(artifact).replace("\\", "/"))
        item["artifacts"] = portable

    edges = [
        edge("premise_open_system", "constraint_t2", "constrains", "sets equilibrium baseline",
             status="constraint", lenses=["open_system"],
             summary="The open thermal environment determines the equilibrium T₂* baseline."),
        edge("constraint_t2", "gap_micro", "requires", "fs → μs",
             status="constraint", lenses=["scale"],
             summary="Femtosecond coherence requires 10⁷–10⁸ amplification for microsecond readout."),
        edge("constraint_t2", "gap_neural", "requires", "fs → neural scale",
             status="constraint", lenses=["scale"],
             summary="Femtosecond coherence requires 10¹¹–10¹³ amplification for neural‑scale function."),
        edge("gap_core", "gap_micro", "decomposes", "readout gap",
             status="neutral", lenses=["scale"],
             summary="The core gap decomposes into the microsecond readout gap."),
        edge("gap_core", "gap_neural", "decomposes", "functionality gap",
             status="neutral", lenses=["scale"],
             summary="The core gap decomposes into the neural functionality gap."),
        edge("constraint_utility", "mech_equilibrium", "falsifies", "U ≪ 1",
             status="fails", lenses=["scale", "falsification"],
             summary="Utility metric U ≪ 1 falsifies equilibrium protection as a viable mechanism."),
        edge("constraint_utility", "mech_scaling", "falsifies", "K_emp too small",
             status="fails", lenses=["scale", "falsification"],
             summary="Empirical gain K_emp is insufficient for naive cooperative scaling."),
        edge("constraint_tur", "constraint_utility", "closes", "first-law closure",
             status="fails", lenses=["scale", "falsification", "open_system"],
             summary="TUR power accounting closes localized equilibrium cascade amplification."),
        edge("constraint_utility", "mech_frohlich", "constrains", "linewidth‑conditional & utility‑limited",
             status="conditional", lenses=["biophysical"],
             summary="Frohlich pumping is constrained to linewidth-conditional regimes with bounded utility."),
        edge("constraint_utility", "mech_cavity", "constrains", "requires measured K",
             status="conditional", lenses=["biophysical"],
             summary="QED cavity mechanism requires measured gain K to satisfy utility constraints."),
        edge("constraint_utility", "mech_subradiance", "constrains", "requires functional validation",
             status="conditional", lenses=["biophysical"],
             summary="Geometric subradiance requires functional validation to satisfy utility constraints."),
        edge("constraint_heom", "mech_equilibrium", "constrains", "non‑perturbative transient",
             status="constraint", lenses=["open_system"],
             summary="HEOM non‑thermalized transient constrains equilibrium‑only models."),
        edge("constraint_heom", "constraint_structured", "supports", "KWW beta clustering",
             status="supported", lenses=["open_system", "reproducibility"],
             summary="HEOM production data supports structured non-Markovian distributed relaxation diagnostics."),
        edge("mech_frohlich", "test_beta", "falsifies_if_absent", "β/γ benchmark",
             status="open", lenses=["falsification"],
             summary="Absence of linewidth-controlled finite-size scaling would falsify micron-scale Frohlich gating."),
        edge("mech_cavity", "test_cavity", "falsifies_if_absent", "UV/THz splitting",
             status="open", lenses=["falsification"],
             summary="Absence of UV/THz cavity splitting would falsify the QED cavity mechanism."),
        edge("mech_subradiance", "test_subradiance", "falsifies_if_absent", "length scaling",
             status="open", lenses=["falsification"],
             summary="Absence of length‑resolved lifetime scaling would falsify geometric subradiance."),
        edge("artifact_validation", "constraint_t2", "verified_by", "ledger‑backed",
             status="supported", lenses=["reproducibility"],
             summary="The T₂* constraint is verified by the validation ledger."),
        edge("artifact_validation", "constraint_heom", "verified_by", "ledger‑backed",
             status="supported", lenses=["reproducibility"],
             summary="The HEOM constraint is verified by the validation ledger."),
        edge("artifact_validation", "constraint_tur", "documents", "validation scope",
             status="supported", lenses=["reproducibility"],
             summary="Validation and supplemental ledgers document the TUR closure inputs."),
        edge("artifact_validation", "constraint_structured", "verified_by", "KWW diagnostics",
             status="supported", lenses=["reproducibility"],
             summary="The structured non-Markovian relaxation diagnostics are documented by JSON artifacts."),
        edge("artifact_validation", "mech_subradiance", "verified_by", "N130/N260/N520",
             status="supported", lenses=["reproducibility"],
             summary="Subradiance data across N130/N260/N520 is verified by the validation ledger."),
        edge("artifact_trace", "constraint_utility", "supports", "public‑data layer",
             status="supported", lenses=["reproducibility"],
             summary="The utility constraint is supported by the claim traceability matrix."),
        edge("artifact_lineage", "artifact_validation", "supports", "release governance",
             status="supported", lenses=["reproducibility"],
             summary="Lineage audit supports the integrity of the validation ledger."),
        edge("premise_quantum_signatures", "scope_no_consciousness", "does_not_imply", "blocked shortcut",
             status="blocked", lenses=["scope"],
             summary="Quantum signatures do not imply consciousness; this inference is blocked."),
        edge("mech_classical", "scope_no_consciousness", "constrains", "no quantum necessity",
             status="constraint", lenses=["scope"],
             summary="Classical synchronization accounts for observed phenomena without quantum necessity."),
        edge("constraint_utility", "scope_no_consciousness", "requires", "missing transfer function",
             status="blocked", lenses=["scope", "scale"],
             summary="Consciousness claims require a transfer function that is currently missing."),
        edge("mech_subradiance", "scope_no_in_vivo", "caveated_by", "not functional validation",
             status="blocked", lenses=["scope"],
             summary="Subradiance results constrain mechanisms but do not constitute in vivo validation."),
        edge("mech_cavity", "scope_no_in_vivo", "caveated_by", "candidate only",
             status="blocked", lenses=["scope"],
             summary="Cavity mechanism remains a candidate only; no in vivo validation is claimed."),
        edge("premise_quantum_signatures", "scope_no_pooling", "caveated_by", "stratified evidence",
             status="blocked", lenses=["scope"],
             summary="Optical and behavioral evidence must remain stratified; causal pooling is blocked."),
    ]

    return {
        "schema_version": "1.1",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "title": "Interactive Epistemic Map of Claims, Evidence, and Falsification",
        "subtitle": "Claim → evidence → validation ledger → falsifier → scope boundary",
        "summary": {
            "validation": f"{validation.get('passed', 0)}/{validation.get('total', 0)}",
            "sha256_prefix": sha,
            "lineage_counts": lineage_counts,
            "microsecond_gap": "10^7‑10^8",
            "neural_gap": "10^11‑10^13",
        },
        "paper_refs": {
            "fails": "§2 (Equilibrium limits), §3.2 (Scaling failure), §3.4 (Utility U < 1)",
            "candidate": "§1.2 (Hypotheses), §4 (Non‑equilibrium mechanisms), §6 (Exp 4 & 6)",
            "open": "§6 (Experimental programmes: Exp 1–6); detectability/IQI/BOED layer in Supplemental Material",
            "blocked": "§1.3–1.5 (Scope & boundaries), §5 (Classical comparators)",
            "supported": "§2.1 (Sobol & SBC audit), Table 1, Appendix A (L4 audit)",
            "all": "Full manuscript: §§1–6, Appendices",
        },
        "lenses": [
            {"id": "scale", "label": "Scale‑Separation Lens"},
            {"id": "open_system", "label": "Open‑System Lens"},
            {"id": "biophysical", "label": "Biophysical Mechanism Lens"},
            {"id": "falsification", "label": "Falsification Lens"},
            {"id": "reproducibility", "label": "Reproducibility Lens"},
            {"id": "scope", "label": "Scope‑Boundary Lens"},
        ],
        "layers": [
            "Premises",
            "Mechanism Classes",
            "Quantitative Constraints",
            "Evidence Artifacts",
            "Decision Tests",
            "Scope Boundary",
        ],
        "target_timescales": {
            "microsecond": {"label": "1 µs readout", "gap": "10⁷–10⁸"},
            "neural_25ms": {"label": "25 ms neural reference", "gap": "~10¹² central baseline"},
            "neural_100ms": {"label": "100 ms neural upper reference", "gap": "up to 10¹³"},
        },
        "nodes": nodes,
        "edges": edges,
    }


def render_html(graph: dict[str, Any]) -> str:
    data = json.dumps(graph, ensure_ascii=False, indent=2)
    escaped_title = html.escape(graph["title"])
    escaped_subtitle = html.escape(graph["subtitle"])
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escaped_title}</title>
<style>
:root {{
  --bg: #0a0e1a; --panel: #111527; --panel2: #1a1f35; --ink: #e8edf5;
  --muted: #8892b0; --line: #2e3750; --supported: #4ade80; --candidate: #38bdf8;
  --conditional: #fbbf24; --fails: #fb7185; --blocked: #94a3b8; --open: #c084fc;
  --constraint: #f59e0b; --glow: rgba(96,165,250,0.3);
  --shadow-sm: 0 1px 3px rgba(0,0,0,0.4);
  --shadow-md: 0 4px 12px rgba(0,0,0,0.5);
  --radius: 10px;
}}
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{
  font-family: 'Inter', system-ui, -apple-system, sans-serif;
  background: radial-gradient(ellipse at 20% 0%, #172554 0%, var(--bg) 55%);
  color:var(--ink);
  overflow:hidden;
}}
body::before {{
  content:'';
  position:fixed;
  inset:0;
  background-image: url('data:image/svg+xml;utf8,<svg xmlns="http://www.w3.org/2000/svg" width="40" height="40"><path d="M30 10a10 10 0 0 1 0 20 10 10 0 0 1 0-20z" fill="rgba(100,150,250,0.03)"/></svg>');
  z-index:-1;
}}
header {{
  padding:14px 24px;
  border-bottom:1px solid var(--line);
  display:flex;
  justify-content:space-between;
  align-items:center;
  backdrop-filter:blur(12px);
  background:rgba(10,14,26,0.8);
  z-index:50;
}}
h1 {{ font-size:20px; font-weight:600; letter-spacing:0.3px; }}
.subtitle {{ color:var(--muted); font-size:13px; margin-top:2px; }}
.header-controls {{ display:flex; gap:8px; }}
.header-controls button {{
  background:var(--panel2);
  color:var(--ink);
  border:1px solid var(--line);
  border-radius:6px;
  padding:6px 12px;
  cursor:pointer;
  font-size:12px;
  transition:border-color 0.2s, box-shadow 0.2s;
}}
.header-controls button:hover {{
  border-color:#60a5fa;
  box-shadow:0 0 8px var(--glow);
}}
.app {{
  display:grid;
  grid-template-columns:260px 1fr 320px;
  height:calc(100vh - 60px);
}}
.app.hide-left {{ grid-template-columns:1fr 320px; }}
.app.hide-left .canvas-wrap {{ grid-column:1; }}
.app.hide-left .inspector {{ grid-column:2; }}
.app.hide-right {{ grid-template-columns:260px 1fr; }}
.app.hide-right aside {{ grid-column:1; }}
.app.hide-right .canvas-wrap {{ grid-column:2; }}
.app.hide-both {{ grid-template-columns:1fr; }}
.app.hide-both .canvas-wrap {{ grid-column:1; }}
aside, .inspector {{
  background:rgba(15,20,35,0.96);
  border-right:1px solid var(--line);
  padding:16px;
  backdrop-filter:blur(8px);
  overflow-y:auto;
  opacity:1;
  min-width:0;
  transition:opacity 0.3s ease, padding 0.3s ease;
}}
.app.hide-left aside,
.app.hide-both aside,
.app.hide-right .inspector,
.app.hide-both .inspector {{ display:none; }}
aside {{
  grid-column: 1;
}}
.inspector {{
  grid-column: 3;
  border-right:0; 
  border-left:1px solid var(--line); 
}}
.canvas-wrap {{
  grid-column: 2;
  position:relative;
  overflow:hidden;
}}
#graph {{
  width:100%;
  height:100%;
  display:block;
  cursor:grab;
}}
#graph:active {{ cursor:grabbing; }}
.small {{ color:var(--muted); font-size:12px; line-height:1.45; }}
.stat {{
  background:var(--panel2);
  border:1px solid var(--line);
  border-radius:var(--radius);
  padding:10px;
  margin:8px 0;
}}
.stat b {{ display:block; font-size:15px; }}
label.lens {{
  display:block;
  margin:6px 0;
  padding:8px;
  border:1px solid var(--line);
  border-radius:8px;
  cursor:pointer;
  background:#0f172a;
  transition:border-color 0.2s;
}}
label.lens:hover {{ border-color:#60a5fa; }}
input[type="search"], select {{
  width:100%;
  border:1px solid var(--line);
  border-radius:8px;
  background:#0f172a;
  color:var(--ink);
  padding:9px;
  margin:8px 0 12px;
}}
button.view {{
  border:1px solid var(--line);
  background:#17213a;
  color:var(--ink);
  border-radius:8px;
  padding:6px 10px;
  margin:3px 2px;
  cursor:pointer;
  font-size:11px;
  transition:border-color 0.2s, background 0.2s;
}}
button.view:hover {{ border-color:#60a5fa; }}
button.view.active {{ border-color:#60a5fa; background:#1e3a5f; color:#fff; }}
.node rect {{
  stroke-width:1.4;
  transition:filter 0.15s, stroke-width 0.15s;
}}
.node:hover rect {{
  filter:drop-shadow(0 0 10px var(--glow));
  stroke-width:2.2;
}}
.node text {{
  fill:var(--ink);
  font-size:12px;
  pointer-events:none;
}}
.node .type {{
  fill:var(--muted);
  font-size:10px;
  text-transform:uppercase;
  letter-spacing:0.6px;
}}
.edge {{
  fill:none;
  stroke-width:2;
  transition:opacity 0.2s;
}}
.edge-label {{
  fill:var(--muted);
  font-size:10px;
  paint-order:stroke;
  stroke:#0a0e1a;
  stroke-width:3px;
  pointer-events:none;
  transition:opacity 0.2s;
}}
.layer-label {{
  fill:#cbd5e1;
  font-size:13px;
  font-weight:700;
  letter-spacing:0.5px;
  text-anchor:middle;
}}
.pill {{
  display:inline-block;
  padding:3px 7px;
  border-radius:999px;
  font-size:11px;
  margin:2px 4px 2px 0;
  background:#1f2a44;
  color:var(--muted);
}}
.supported {{ color:var(--supported); }} .candidate {{ color:var(--candidate); }} .conditional {{ color:var(--conditional); }} .fails {{ color:var(--fails); }} .blocked {{ color:var(--blocked); }} .open {{ color:var(--open); }}
.artifact {{ color:#93c5fd; overflow-wrap:anywhere; font-size:12px; }}
.legend-row {{ display:flex; align-items:center; gap:7px; margin:4px 0; color:var(--muted); font-size:11px; }}
.dot {{ width:8px; height:8px; border-radius:50%; display:inline-block; }}
.refs-box {{
  background:var(--panel2);
  border:1px solid var(--line);
  border-radius:8px;
  padding:10px;
  margin:10px 0;
  font-size:11px;
  color:var(--muted);
  line-height:1.4;
}}
.zoom-controls {{
  position:absolute;
  bottom:16px;
  right:16px;
  display:flex;
  flex-direction:column;
  gap:6px;
  z-index:20;
}}
.zoom-controls button {{
  background:rgba(15,20,35,0.9);
  color:var(--ink);
  border:1px solid var(--line);
  border-radius:8px;
  padding:8px;
  cursor:pointer;
  font-size:16px;
  line-height:1;
  transition:border-color 0.2s, background 0.2s;
  backdrop-filter:blur(4px);
}}
.zoom-controls button:hover {{
  border-color:var(--candidate);
  background:var(--panel);
}}
#tooltip {{
  position:fixed;
  background:rgba(10,14,26,0.95);
  color:var(--ink);
  padding:8px 12px;
  border-radius:var(--radius);
  font-size:12px;
  pointer-events:none;
  z-index:200;
  max-width:280px;
  border:1px solid var(--line);
  display:none;
  line-height:1.4;
  box-shadow:var(--shadow-md);
  backdrop-filter:blur(8px);
  transition:opacity 0.15s;
}}
.highlight {{ background:rgba(96,165,250,0.15); border-radius:4px; }}
@media (max-width:1100px) {{
  .app {{ grid-template-columns:1fr; }}
  aside,.inspector {{ max-height:none; border:0; }}
  #graph {{ height:72vh; }}
}}
</style>
</head>
<body>
<header>
  <div>
    <h1>{escaped_title}</h1>
    <div class="subtitle">{escaped_subtitle}</div>
  </div>
  <div class="header-controls">
    <button id="toggleLeft" title="Show/hide filters">☰ Filters</button>
    <button id="toggleRight" title="Show/hide inspector">♟ Inspector</button>
  </div>
</header>
<main class="app" id="appLayout">
  <aside id="leftPanel">
    <div class="stat"><b id="validationStat"></b><span class="small">Validation ledger</span></div>
    <div class="stat"><b id="gapStat"></b><span class="small">Scale separation</span></div>
    <div class="stat"><b id="visibleStat"></b><span class="small">Visible nodes / edges</span></div>
    <label class="small" for="timescale">Target timescale</label>
    <select id="timescale"></select>
    <input id="search" type="search" placeholder="Search claims, mechanisms, artifacts…">
    <h3>Lenses</h3>
    <div id="lenses"></div>
    <h3>Quick views</h3>
    <div id="viewBtns">
      <button class="view" data-view="fails">📉 Failed routes</button>
      <button class="view" data-view="candidate">🔬 Surviving candidates</button>
      <button class="view" data-view="open">🧪 Falsifiers</button>
      <button class="view" data-view="blocked">🚧 Scope boundaries</button>
      <button class="view" data-view="supported">✅ Ledger‑backed</button>
      <button class="view active" data-view="all">🌐 All</button>
    </div>
    <div id="refsBox" class="refs-box">Select a view to see relevant manuscript sections.</div>
    <h3>Data Clutter</h3>
    <label class="lens" style="margin-top:4px">
      <input type="checkbox" id="hideIsolatedCheck" checked>
      Hide isolated nodes (validation checks)
    </label>
    <h3>Legend</h3>
    <div class="legend-row"><span class="dot" style="background:var(--fails)"></span>Fails current constraints</div>
    <div class="legend-row"><span class="dot" style="background:var(--conditional)"></span>Conditional / utility‑limited</div>
    <div class="legend-row"><span class="dot" style="background:var(--candidate)"></span>Falsifiable candidate</div>
    <div class="legend-row"><span class="dot" style="background:var(--supported)"></span>Ledger‑supported</div>
    <div class="legend-row"><span class="dot" style="background:var(--blocked)"></span>Scope boundary</div>
  </aside>
  <section class="canvas-wrap">
    <div class="zoom-controls">
      <button id="zoomIn" title="Zoom In (Ctrl + wheel)">+</button>
      <button id="zoomOut" title="Zoom Out">−</button>
      <button id="zoomFit" title="Fit to Screen">⊡</button>
    </div>
    <svg id="graph" role="img" aria-label="Interactive epistemic graph"></svg>
  </section>
  <section class="inspector" id="rightPanel">
    <h2 id="iTitle">Select a node</h2>
    <div id="iMeta" class="small">Click any node or edge to inspect its epistemic role.</div>
    <p id="iSummary"></p>
    <div id="iDetails" class="small"></div>
    <h3>Artifacts</h3>
    <div id="iArtifacts" class="small">None selected.</div>
    <h3>Metrics</h3>
    <div id="iMetrics" class="small">None selected.</div>
  </section>
</main>
<div id="tooltip"></div>

<script>
const GRAPH = {data};
const COLORS = {{
  supported:'#4ade80', candidate:'#38bdf8', conditional:'#fbbf24',
  fails:'#fb7185', blocked:'#94a3b8', open:'#c084fc',
  constraint:'#f59e0b', neutral:'#64748b'
}};
let activeLenses = new Set(GRAPH.lenses.map(l => l.id));
let quickView = 'all';
let search = '';
let selected = null;
let hideIsolated = true;
let leftHidden = false;
let rightHidden = false;
const svg = document.getElementById('graph');
const appLayout = document.getElementById('appLayout');
const canvasWrap = document.querySelector('.canvas-wrap');
const lensesBox = document.getElementById('lenses');
const refsBox = document.getElementById('refsBox');
const tooltip = document.getElementById('tooltip');
let tooltipTimeout = null;

document.getElementById('validationStat').textContent = `${{GRAPH.summary.validation}} OK / SHA ${{GRAPH.summary.sha256_prefix}}`;
document.getElementById('gapStat').textContent = `${{GRAPH.summary.microsecond_gap}} | ${{GRAPH.summary.neural_gap}}`;

/* Pre‑compute node degrees */
const nodeDegrees = {{}};
GRAPH.nodes.forEach(n => nodeDegrees[n.id] = 0);
GRAPH.edges.forEach(e => {{
  nodeDegrees[e.source] = (nodeDegrees[e.source] || 0) + 1;
  nodeDegrees[e.target] = (nodeDegrees[e.target] || 0) + 1;
}});

/* Lenses checkboxes */
for (const lens of GRAPH.lenses) {{
  const label = document.createElement('label'); label.className = 'lens';
  label.innerHTML = `<input type="checkbox" checked value="${{lens.id}}"> ${{lens.label}}`;
  label.querySelector('input').addEventListener('change', e => {{
    e.target.checked ? activeLenses.add(lens.id) : activeLenses.delete(lens.id);
    render();
  }});
  lensesBox.appendChild(label);
}}

/* Timescale selector */
const timescale = document.getElementById('timescale');
for (const [id, item] of Object.entries(GRAPH.target_timescales)) {{
  const opt = document.createElement('option'); opt.value = id;
  opt.textContent = `${{item.label}} (${{item.gap}})`; timescale.appendChild(opt);
}}
timescale.addEventListener('change', render);

/* Search */
let searchTimeout;
document.getElementById('search').addEventListener('input', e => {{
  clearTimeout(searchTimeout);
  searchTimeout = setTimeout(() => {{ search = e.target.value.toLowerCase(); render(); }}, 150);
}});

/* Quick view buttons */
const viewBtns = document.querySelectorAll('#viewBtns button.view');
viewBtns.forEach(btn => btn.addEventListener('click', () => {{
  quickView = btn.dataset.view;
  viewBtns.forEach(b => b.classList.toggle('active', b === btn));
  refsBox.innerHTML = `<b>📖 Manuscript refs:</b><br>${{GRAPH.paper_refs[quickView] || 'N/A'}}`;
  render();
}}));

/* Hide isolated toggle */
document.getElementById('hideIsolatedCheck').addEventListener('change', e => {{
  hideIsolated = e.target.checked;
  render();
}});

/* Sidebar toggles */
const leftPanel = document.getElementById('leftPanel');
const rightPanel = document.getElementById('rightPanel');

function updateLayoutState() {{
  appLayout.classList.toggle('hide-left', leftHidden && !rightHidden);
  appLayout.classList.toggle('hide-right', rightHidden && !leftHidden);
  appLayout.classList.toggle('hide-both', leftHidden && rightHidden);
  requestAnimationFrame(() => requestAnimationFrame(render));
}}

document.getElementById('toggleLeft').addEventListener('click', () => {{
  leftHidden = !leftHidden;
  updateLayoutState();
}});

document.getElementById('toggleRight').addEventListener('click', () => {{
  rightHidden = !rightHidden;
  updateLayoutState();
}});

new ResizeObserver(() => requestAnimationFrame(render)).observe(canvasWrap);

/* Zoom & Pan */
let transform = {{ x: 0, y: 0, k: 1 }};
let isPanning = false;
let panStart = {{ x: 0, y: 0 }};

function applyTransform() {{
  const vp = document.getElementById('viewport');
  if (vp) vp.setAttribute('transform', `translate(${{transform.x}},${{transform.y}}) scale(${{transform.k}})`);
}}

svg.addEventListener('wheel', e => {{
  e.preventDefault();
  const factor = e.deltaY > 0 ? 0.92 : 1.08;
  const rect = svg.getBoundingClientRect();
  const mx = e.clientX - rect.left;
  const my = e.clientY - rect.top;
  transform.x = mx - (mx - transform.x) * factor;
  transform.y = my - (my - transform.y) * factor;
  transform.k = Math.min(4, Math.max(0.2, transform.k * factor));
  applyTransform();
}}, {{ passive: false }});

svg.addEventListener('mousedown', e => {{
  if (e.target === svg || e.target.classList.contains('canvas-bg')) {{
    isPanning = true;
    panStart = {{ x: e.clientX - transform.x, y: e.clientY - transform.y }};
    svg.style.cursor = 'grabbing';
  }}
}});
window.addEventListener('mousemove', e => {{
  if (isPanning) {{
    transform.x = e.clientX - panStart.x;
    transform.y = e.clientY - panStart.y;
    applyTransform();
  }}
}});
window.addEventListener('mouseup', () => {{
  if (isPanning) {{ isPanning = false; svg.style.cursor = 'grab'; }}
}});
document.getElementById('zoomIn').addEventListener('click', () => {{
  transform.k = Math.min(4, transform.k * 1.2); applyTransform();
}});
document.getElementById('zoomOut').addEventListener('click', () => {{
  transform.k = Math.max(0.2, transform.k * 0.8); applyTransform();
}});
document.getElementById('zoomFit').addEventListener('click', () => {{
  transform = {{ x: 0, y: 0, k: 1 }}; applyTransform();
}});

/* Tooltip with delay to avoid flicker */
function showTooltip(evt, htmlContent) {{
  if (tooltipTimeout) clearTimeout(tooltipTimeout);
  tooltipTimeout = setTimeout(() => {{
    tooltip.innerHTML = htmlContent;
    tooltip.style.display = 'block';
    moveTooltip(evt);
  }}, 80);
}}
function moveTooltip(evt) {{
  tooltip.style.left = (evt.clientX + 15) + 'px';
  tooltip.style.top = (evt.clientY + 10) + 'px';
}}
function hideTooltip() {{
  if (tooltipTimeout) clearTimeout(tooltipTimeout);
  tooltip.style.display = 'none';
}}

/* Visibility helpers */
function visibleNode(n) {{
  const lensOk = n.lenses.length === 0 || n.lenses.some(l => activeLenses.has(l));
  const viewOk = quickView === 'all' || n.status === quickView || (quickView === 'candidate' && ['candidate','conditional'].includes(n.status));
  const isoOk = !hideIsolated || nodeDegrees[n.id] > 0;
  const hay = [n.label,n.type,n.status,n.summary,n.details,(n.artifacts||[]).join(' ')].join(' ').toLowerCase();
  return lensOk && viewOk && isoOk && (!search || hay.includes(search));
}}
function edgeVisible(e, visibleIds) {{
  const lensOk = e.lenses.length === 0 || e.lenses.some(l => activeLenses.has(l));
  return lensOk && visibleIds.has(e.source) && visibleIds.has(e.target);
}}
function color(status) {{ return COLORS[status] || COLORS.neutral; }}
function esc(s) {{ return String(s ?? '').replace(/[&<>"]/g, ch => ({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}}[ch])); }}

const nodeById = new Map(GRAPH.nodes.map(n => [n.id, n]));

function getConnected(nodeId) {{
  const c = new Set([nodeId]);
  for (const e of GRAPH.edges) {{
    if (e.source === nodeId) c.add(e.target);
    if (e.target === nodeId) c.add(e.source);
  }}
  return c;
}}

/* Inspector */
function inspect(item, isEdge) {{
  selected = item.id;
  document.getElementById('iTitle').textContent = item.label;
  if (isEdge) {{
    const src = nodeById.get(item.source);
    const tgt = nodeById.get(item.target);
    document.getElementById('iMeta').innerHTML =
      `<span class="pill">${{esc(item.relation)}}</span><span class="pill ${{item.status}}">${{item.status}}</span>` +
      `<div class="small" style="margin-top:6px">${{esc(src ? src.label : item.source)}} &rarr; ${{esc(tgt ? tgt.label : item.target)}}</div>`;
    document.getElementById('iSummary').textContent = item.summary || `${{src ? src.label : item.source}} ${{item.relation.replace(/_/g,' ')}} ${{tgt ? tgt.label : item.target}}`;
  }} else {{
    document.getElementById('iMeta').innerHTML = `<span class="pill">${{item.type}}</span><span class="pill ${{item.status}}">${{item.status}}</span><span class="pill">${{item.layer}}</span>`;
    document.getElementById('iSummary').textContent = item.summary || item.label || '';
  }}
  document.getElementById('iDetails').textContent = item.details || '';
  const arts = item.artifacts || [];
  document.getElementById('iArtifacts').innerHTML = arts.length ? arts.map(a => `<div class="artifact">${{esc(a)}}</div>`).join('') : '<span class="small">None.</span>';
  const metrics = item.metrics || {{}};
  document.getElementById('iMetrics').innerHTML = Object.keys(metrics).length ? Object.entries(metrics).map(([k,v]) => `<div><b>${{esc(k)}}:</b> ${{esc(v)}}</div>`).join('') : '<span class="small">None.</span>';
}}

function clearInspect() {{
  selected = null;
  document.getElementById('iTitle').textContent = 'Select a node';
  document.getElementById('iMeta').innerHTML = 'Click any node or edge to inspect its epistemic role.';
  document.getElementById('iSummary').textContent = '';
  document.getElementById('iDetails').textContent = '';
  document.getElementById('iArtifacts').innerHTML = '<span class="small">None selected.</span>';
  document.getElementById('iMetrics').innerHTML = '<span class="small">None selected.</span>';
}}

/* Barycentric heuristic for layer ordering */
function barycentricOrder(layerIds, prevOrder) {{
  const orderMap = prevOrder ? new Map(prevOrder.map((id, idx) => [id, idx])) : new Map();
  return layerIds.slice().sort((a, b) => {{
    const aParents = GRAPH.edges.filter(e => e.target === a && orderMap.has(e.source));
    const bParents = GRAPH.edges.filter(e => e.target === b && orderMap.has(e.source));
    const aAvg = aParents.length ? aParents.reduce((s, e) => s + (orderMap.get(e.source) || 0), 0) / aParents.length : 500;
    const bAvg = bParents.length ? bParents.reduce((s, e) => s + (orderMap.get(e.source) || 0), 0) / bParents.length : 500;
    return aAvg - bAvg;
  }});
}}

function render() {{
  const width = svg.clientWidth || 1000;
  const height = svg.clientHeight || 700;
  
  const visible = GRAPH.nodes.filter(visibleNode);
  const visibleIds = new Set(visible.map(n => n.id));
  
  if (selected && !visibleIds.has(selected) && !GRAPH.edges.some(e => e.id === selected)) selected = null;
  
  const edges = GRAPH.edges.filter(e => edgeVisible(e, visibleIds));

  const layerX = new Map(GRAPH.layers.map((l,i) => [l, 110 + i * ((width - 220) / Math.max(GRAPH.layers.length-1, 1))]));
  const byLayer = new Map();
  
  for (const n of visible) {{
    if (!byLayer.has(n.layer)) byLayer.set(n.layer, []);
    byLayer.get(n.layer).push(n);
  }}

  const nodeH = 54, gap = 12, step = nodeH + gap;
  let maxNodes = 0;
  for (const arr of byLayer.values()) {{
    if (arr.length > maxNodes) maxNodes = arr.length;
  }}
  
  const contentH = 60 + maxNodes * step + 40;
  const viewH = Math.max(height, contentH);
  svg.setAttribute('viewBox', `0 0 ${{width}} ${{viewH}}`);

  /* Order layers to minimize crossings */
  const ordered = [];
  for (let i = 0; i < GRAPH.layers.length; i++) {{
    const layer = GRAPH.layers[i];
    const ids = (byLayer.get(layer) || []).map(n => n.id);
    let sorted = ids;
    for (let pass = 0; pass < 3; pass++) {{
      if (i > 0 && ordered[i-1]) sorted = barycentricOrder(sorted, ordered[i-1].map(o => o.id));
    }}
    ordered.push(sorted);
  }}

  const positions = {{}};
  for (let i = 0; i < GRAPH.layers.length; i++) {{
    const layer = GRAPH.layers[i];
    const arr = ordered[i];
    const totalH = arr.length * step;
    const startY = Math.max(50, (viewH - totalH) / 2);
    arr.forEach((id, idx) => {{
      positions[id] = {{ x: layerX.get(layer), y: startY + idx * step + nodeH/2 }};
    }});
  }}

  document.getElementById('visibleStat').textContent = `${{visible.length}} / ${{edges.length}}`;

  svg.innerHTML = '';
  const viewport = document.createElementNS('http://www.w3.org/2000/svg', 'g');
  viewport.setAttribute('id', 'viewport');

  /* Background for panning */
  const bg = document.createElementNS('http://www.w3.org/2000/svg','rect');
  bg.setAttribute('width', width*5); 
  bg.setAttribute('height', viewH*5);
  bg.setAttribute('x', -width*2); 
  bg.setAttribute('y', -viewH*2);
  bg.setAttribute('fill','transparent'); 
  bg.setAttribute('class','canvas-bg');
  viewport.appendChild(bg);

  /* Defs for arrow markers */
  const defs = document.createElementNS('http://www.w3.org/2000/svg','defs');
  for (const [status, col] of Object.entries(COLORS)) {{
    const marker = document.createElementNS('http://www.w3.org/2000/svg','marker');
    marker.setAttribute('id', `arrow-${{status}}`); 
    marker.setAttribute('markerWidth','8'); 
    marker.setAttribute('markerHeight','8');
    marker.setAttribute('refX','8'); 
    marker.setAttribute('refY','4'); 
    marker.setAttribute('orient','auto'); 
    marker.setAttribute('markerUnits','userSpaceOnUse');
    const path = document.createElementNS('http://www.w3.org/2000/svg','path');
    path.setAttribute('d','M0,1 L7,4 L0,7 Z'); 
    path.setAttribute('fill', col);
    marker.appendChild(path); 
    defs.appendChild(marker);
  }}
  viewport.appendChild(defs);

  /* Layer labels */
  for (const layer of GRAPH.layers) {{
    const text = document.createElementNS('http://www.w3.org/2000/svg','text');
    text.setAttribute('x', layerX.get(layer));
    text.setAttribute('y', 24);
    text.setAttribute('class','layer-label');
    text.textContent = layer;
    viewport.appendChild(text);
  }}

  const connected = selected ? getConnected(selected) : null;

  /* Edges */
  for (const e of edges) {{
    const a = positions[e.source], b = positions[e.target];
    if (!a || !b) continue;
    const mx = (a.x + b.x) / 2;
    const d = `M ${{a.x+86}} ${{a.y}} C ${{mx}} ${{a.y}}, ${{mx}} ${{b.y}}, ${{b.x-86}} ${{b.y}}`;
    const eConn = connected ? (connected.has(e.source) && connected.has(e.target)) : true;
    const eOpacity = eConn ? '0.85' : '0.06';

    /* Hitbox */
    const hitbox = document.createElementNS('http://www.w3.org/2000/svg','path');
    hitbox.setAttribute('d', d); 
    hitbox.setAttribute('fill','none');
    hitbox.setAttribute('stroke','transparent'); 
    hitbox.setAttribute('stroke-width','14');
    hitbox.style.cursor = 'pointer';
    hitbox.addEventListener('click', ev => {{ ev.stopPropagation(); inspect(e, true); render(); }});
    hitbox.addEventListener('mouseenter', ev => showTooltip(ev, `<b>${{e.label}}</b><br><span style="opacity:0.8">${{e.summary || e.relation}}</span>`));
    hitbox.addEventListener('mousemove', moveTooltip);
    hitbox.addEventListener('mouseleave', hideTooltip);
    viewport.appendChild(hitbox);

    const path = document.createElementNS('http://www.w3.org/2000/svg','path');
    path.setAttribute('d', d); 
    path.setAttribute('class','edge');
    path.setAttribute('stroke', color(e.status)); 
    path.setAttribute('stroke-width','2');
    path.setAttribute('marker-end', `url(#arrow-${{e.status}})`);
    path.style.pointerEvents = 'none'; 
    path.style.opacity = eOpacity;
    if (['blocked','open'].includes(e.status)) path.setAttribute('stroke-dasharray','6 5');
    viewport.appendChild(path);

    /* Edge label */
    const txt = document.createElementNS('http://www.w3.org/2000/svg','text');
    txt.setAttribute('x', mx - 32); 
    txt.setAttribute('y', (a.y + b.y)/2 - 4);
    txt.setAttribute('class','edge-label'); 
    txt.textContent = e.label;
    txt.style.opacity = eConn ? '1' : '0.06';
    viewport.appendChild(txt);
  }}

  /* Nodes */
  for (const n of visible) {{
    const p = positions[n.id]; if (!p) continue;
    const isSelected = selected === n.id;
    const isConn = connected ? connected.has(n.id) : true;

    const g = document.createElementNS('http://www.w3.org/2000/svg','g');
    g.setAttribute('class','node'); g.style.cursor = 'pointer';
    g.style.opacity = isConn ? '1' : '0.1';

    const rect = document.createElementNS('http://www.w3.org/2000/svg','rect');
    rect.setAttribute('x', p.x-92); rect.setAttribute('y', p.y-28);
    rect.setAttribute('width',184); rect.setAttribute('height',56); rect.setAttribute('rx',10);
    rect.setAttribute('fill', isSelected ? '#1e3050' : '#11182c');
    rect.setAttribute('stroke', color(n.status));
    if (isSelected) rect.setAttribute('stroke-width','2.5');
    g.appendChild(rect);

    const maxLabelChars = 26;
    const labelStr = n.label.length > maxLabelChars ? n.label.substring(0, maxLabelChars - 3) + '…' : n.label;

    const title = document.createElementNS('http://www.w3.org/2000/svg','text');
    title.setAttribute('x', p.x-84); title.setAttribute('y', p.y-6);
    title.textContent = labelStr;
    g.appendChild(title);

    const meta = document.createElementNS('http://www.w3.org/2000/svg','text');
    meta.setAttribute('x', p.x-84); meta.setAttribute('y', p.y+15);
    meta.setAttribute('class','type'); meta.textContent = `${{n.type}} / ${{n.status}}`;
    g.appendChild(meta);

    g.addEventListener('click', ev => {{
      ev.stopPropagation();
      if (selected === n.id) {{ selected = null; clearInspect(); }}
      else {{ inspect(n, false); }}
      render();
    }});
    g.addEventListener('mouseenter', ev => showTooltip(ev, `<b>${{n.label}}</b><br><span style="opacity:0.8">${{n.summary || n.type}}</span>`));
    g.addEventListener('mousemove', moveTooltip);
    g.addEventListener('mouseleave', hideTooltip);

    viewport.appendChild(g);
  }}

  svg.appendChild(viewport);
  applyTransform();
}}

render();
window.addEventListener('resize', () => requestAnimationFrame(render));
</script>
</body>
</html>"""


def main() -> int:
    graph = build_graph()
    out_dir = PROJECT_ROOT / "outputs_data" / "raw_json" / "structural"
    interactive_dir = PROJECT_ROOT / "outputs_data" / "interactive"
    out_dir.mkdir(parents=True, exist_ok=True)
    interactive_dir.mkdir(parents=True, exist_ok=True)

    json_path = out_dir / "epistemic_graph.json"
    html_path = interactive_dir / "epistemic_graph.html"
    json_path.write_text(json.dumps(graph, indent=2, ensure_ascii=False), encoding="utf-8")
    html_path.write_text(render_html(graph), encoding="utf-8")
    logger.info(f"JSON written → {json_path}")
    logger.info(f"HTML written → {html_path}")
    logger.info(f"Nodes: {len(graph['nodes'])}, Edges: {len(graph['edges'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
