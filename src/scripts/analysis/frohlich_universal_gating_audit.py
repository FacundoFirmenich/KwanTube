"""Dimensional audit for finite-size Frohlich gating.

This script separates two physically distinct crossover criteria:

1. Carrier-wavelength criterion: L_omega = pi * v_g / omega_F.
2. Linewidth-continuum criterion: L_gamma = pi * v_g / gamma.

The first asks whether the driven carrier wavelength fits in the polymer.
The second asks whether finite-size mode spacing is smaller than the dissipative
linewidth, so that a continuum bath approximation is justified.

The calculation is deliberately lightweight and model-agnostic: it does not
claim that Frohlich condensation occurs, only where a length-gated continuum
argument can be dimensionally self-consistent.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class PolymerCase:
    name: str
    spacing_m: float
    group_velocity_m_s: float
    note: str


POLYMERS = [
    PolymerCase(
        name="microtubule",
        spacing_m=8.0e-9,
        group_velocity_m_s=2.0e3,
        note="Tubulin dimer repeat; representative protein acoustic velocity.",
    ),
    PolymerCase(
        name="F-actin",
        spacing_m=2.75e-9,
        group_velocity_m_s=2.0e3,
        note="Actin axial monomer rise; velocity kept equal for class comparison.",
    ),
    PolymerCase(
        name="collagen",
        spacing_m=0.286e-9,
        group_velocity_m_s=3.0e3,
        note="Peptide axial rise scale; slightly higher axial sound velocity.",
    ),
    PolymerCase(
        name="generic_dipolar_chain",
        spacing_m=1.0e-9,
        group_velocity_m_s=2.0e3,
        note="Dimensionless reference chain.",
    ),
]


FREQUENCY_HZ = [1.0e6, 1.0e8, 1.0e9, 1.0e10, 1.0e11, 1.0e12]
LINEWIDTH_HZ = [1.0e6, 1.0e7, 1.0e8, 1.0e9, 1.0e10, 1.0e11]
TARGET_LENGTHS_M = [10e-9, 100e-9, 1e-6, 10e-6, 100e-6]


def rad_s(freq_hz: float) -> float:
    return 2.0 * math.pi * freq_hz


def crossover_length(group_velocity_m_s: float, angular_rate_s: float) -> float:
    return math.pi * group_velocity_m_s / angular_rate_s


def beta_eff_interpolation(length_m: float, crossover_m: float) -> float:
    """Smooth proxy for beta: beta->1 below crossover, beta->0 above it."""
    x = max(length_m / crossover_m, 1e-30)
    return 1.0 / (1.0 + x * x)


def required_linewidth_hz_for_target(group_velocity_m_s: float, target_length_m: float) -> float:
    # L_gamma = pi v / gamma_rad = v / (2 gamma_hz)
    return group_velocity_m_s / (2.0 * target_length_m)


def audit_case(case: PolymerCase) -> dict:
    carrier = []
    for freq in FREQUENCY_HZ:
        l_omega = crossover_length(case.group_velocity_m_s, rad_s(freq))
        carrier.append(
            {
                "frequency_hz": freq,
                "L_omega_m": l_omega,
                "L_omega_um": l_omega * 1e6,
                "monomers_at_L_omega": l_omega / case.spacing_m,
            }
        )

    linewidth = []
    for gamma_hz in LINEWIDTH_HZ:
        l_gamma = crossover_length(case.group_velocity_m_s, rad_s(gamma_hz))
        linewidth.append(
            {
                "linewidth_hz": gamma_hz,
                "L_gamma_m": l_gamma,
                "L_gamma_um": l_gamma * 1e6,
                "monomers_at_L_gamma": l_gamma / case.spacing_m,
            }
        )

    beta_table = []
    reference_gamma_hz = 1.0e8
    l_gamma_ref = crossover_length(case.group_velocity_m_s, rad_s(reference_gamma_hz))
    for length in TARGET_LENGTHS_M:
        beta_table.append(
            {
                "length_m": length,
                "length_um": length * 1e6,
                "beta_eff_gamma_100MHz": beta_eff_interpolation(length, l_gamma_ref),
            }
        )

    return {
        "case": asdict(case),
        "carrier_wavelength_criterion": carrier,
        "linewidth_continuum_criterion": linewidth,
        "beta_eff_proxy": beta_table,
        "linewidth_required_for_10um_hz": required_linewidth_hz_for_target(
            case.group_velocity_m_s, 10e-6
        ),
    }


def main() -> None:
    project_root = Path(__file__).resolve().parents[3]
    out_dir = project_root / "outputs_data" / "raw_json" / "nonequilibrium"
    out_dir.mkdir(parents=True, exist_ok=True)

    results = {
        "audit": "frohlich_universal_gating_dimensional_audit",
        "criteria": {
            "carrier_wavelength": "L_omega = pi v_g / omega_F = v_g / (2 f_F)",
            "linewidth_continuum": "L_gamma = pi v_g / gamma = v_g / (2 gamma_Hz)",
            "interpretation": "Micron-scale gating at THz carrier frequencies is not wavelength-controlled; it requires a narrow effective linewidth or a lower-frequency collective mode.",
        },
        "cases": [audit_case(case) for case in POLYMERS],
    }

    out_path = out_dir / "frohlich_universal_gating_audit.json"
    out_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    print(f"[frohlich-audit] wrote {out_path}")
    for item in results["cases"]:
        case = item["case"]
        thz = next(row for row in item["carrier_wavelength_criterion"] if row["frequency_hz"] == 1.0e11)
        print(
            f"{case['name']}: L_omega(0.1 THz)={thz['L_omega_um']:.4g} um; "
            f"gamma required for L_gamma=10 um={item['linewidth_required_for_10um_hz']:.3g} Hz"
        )


if __name__ == "__main__":
    main()
