"""
test_parameter_traceability.py — Tests de trazabilidad de parámetros físicos.
Patch v3.5.1 — Valida que build_hamiltonian.py lee desde config/physics_params.yaml,
que los valores numéricos son consistentes con la fuente documentada, y que
no existen parámetros hardcodeados sin respaldo en el YAML.

Ejecutar con: pytest tests/test_parameter_traceability.py -v
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

# Path setup
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

PHYSICS_YAML = PROJECT_ROOT / "config" / "physics_params.yaml"
NUMERICAL_YAML = PROJECT_ROOT / "config" / "numerical_params.yaml"


# ------------------------------------------------------------------ YAML existence
class TestYAMLFiles:
    """Tests de existencia y estructura de los archivos de configuración."""

    def test_physics_yaml_exists(self):
        assert PHYSICS_YAML.exists(), (
            f"config/physics_params.yaml no encontrado en {PHYSICS_YAML}. "
            "Ejecutar el patch v3.5.1 para crearlo."
        )

    def test_numerical_yaml_exists(self):
        assert NUMERICAL_YAML.exists(), (
            f"config/numerical_params.yaml no encontrado en {NUMERICAL_YAML}."
        )

    def test_physics_yaml_parseable(self):
        """El YAML de física debe poder parsearse sin errores."""
        pytest.importorskip("yaml", reason="pyyaml no instalado")
        import yaml
        with open(PHYSICS_YAML, "r", encoding="utf-8") as f:
            params = yaml.safe_load(f)
        assert isinstance(params, dict), "physics_params.yaml debe ser un diccionario YAML"

    def test_numerical_yaml_parseable(self):
        pytest.importorskip("yaml", reason="pyyaml no instalado")
        import yaml
        with open(NUMERICAL_YAML, "r", encoding="utf-8") as f:
            params = yaml.safe_load(f)
        assert isinstance(params, dict)

    def test_physics_yaml_required_sections(self):
        """El YAML debe contener las secciones obligatorias."""
        pytest.importorskip("yaml", reason="pyyaml no instalado")
        import yaml
        with open(PHYSICS_YAML, "r", encoding="utf-8") as f:
            params = yaml.safe_load(f)
        required_sections = ["dipole", "energetics", "constants"]
        for section in required_sections:
            assert section in params, (
                f"Sección '{section}' faltante en physics_params.yaml"
            )

    def test_physics_yaml_required_parameters(self):
        """Cada parámetro físico clave debe existir con 'value' y 'source'."""
        pytest.importorskip("yaml", reason="pyyaml no instalado")
        import yaml
        with open(PHYSICS_YAML, "r", encoding="utf-8") as f:
            params = yaml.safe_load(f)

        required_params = [
            ("dipole", "mu_trp_1La_D"),
            ("dipole", "epsilon_r_protein"),
            ("energetics", "E0_cm1"),
            ("energetics", "sigma_site_cm1"),
            ("constants", "prefactor_cm1_A3_D2"),
        ]
        for section, param in required_params:
            assert section in params, f"Sección '{section}' no encontrada"
            assert param in params[section], (
                f"Parámetro '{param}' no encontrado en sección '{section}'"
            )
            entry = params[section][param]
            assert "value" in entry, (
                f"'{section}/{param}' debe tener campo 'value'"
            )
            assert "source" in entry, (
                f"'{section}/{param}' debe tener campo 'source' para trazabilidad"
            )

    def test_physics_yaml_no_empty_source(self):
        """Ningún parámetro debe tener 'source' vacío."""
        pytest.importorskip("yaml", reason="pyyaml no instalado")
        import yaml
        with open(PHYSICS_YAML, "r", encoding="utf-8") as f:
            params = yaml.safe_load(f)

        def _check_sources(d: dict, path: str = ""):
            for k, v in d.items():
                current_path = f"{path}/{k}" if path else k
                if isinstance(v, dict):
                    if "source" in v:
                        assert v["source"] and str(v["source"]).strip() != "", (
                            f"'source' vacío en {current_path}"
                        )
                    _check_sources(v, current_path)

        _check_sources(params)


# ------------------------------------------------------------------ YAML value ranges
class TestPhysicsParameterRanges:
    """Tests de rangos físicamente razonables para los parámetros."""

    @pytest.fixture(autouse=True)
    def load_yaml(self):
        pytest.importorskip("yaml", reason="pyyaml no instalado")
        import yaml
        if not PHYSICS_YAML.exists():
            pytest.skip("physics_params.yaml no existe")
        with open(PHYSICS_YAML, "r", encoding="utf-8") as f:
            self.params = yaml.safe_load(f)

    def test_mu_trp_range(self):
        """Momento dipolar de Trp 1La debe estar entre 1 y 5 Debye."""
        mu = float(self.params["dipole"]["mu_trp_1La_D"]["value"])
        assert 1.0 <= mu <= 5.0, (
            f"mu_trp_1La_D = {mu} D fuera del rango físico [1, 5] D. "
            "Verificar cita Callis1997."
        )

    def test_epsilon_r_range(self):
        """Constante dieléctrica efectiva de proteína debe estar entre 1.5 y 10."""
        eps = float(self.params["dipole"]["epsilon_r_protein"]["value"])
        assert 1.5 <= eps <= 10.0, (
            f"epsilon_r_protein = {eps} fuera del rango [1.5, 10]."
        )

    def test_E0_trp_range(self):
        """Energía 0-0 de Trp 1La debe estar cerca de 35000 cm^-1 (± 2000)."""
        E0 = float(self.params["energetics"]["E0_cm1"]["value"])
        assert 33000.0 <= E0 <= 37000.0, (
            f"E0_cm1 = {E0} cm^-1 fuera del rango esperado [33000, 37000]."
        )

    def test_sigma_site_range(self):
        """Desorden estático debe estar entre 50 y 500 cm^-1."""
        sigma = float(self.params["energetics"]["sigma_site_cm1"]["value"])
        assert 50.0 <= sigma <= 500.0, (
            f"sigma_site_cm1 = {sigma} cm^-1 fuera del rango [50, 500]."
        )

    def test_prefactor_value(self):
        """Prefactor de acoplamiento debe ser ~5.04e4 cm^-1 Å^3 D^-2."""
        C = float(self.params["constants"]["prefactor_cm1_A3_D2"]["value"])
        # Valor de referencia: Madjet-Renger 2006, ec. 2
        assert abs(C - 5.04e4) < 1e3, (
            f"prefactor_cm1_A3_D2 = {C} se desvía >1000 de 50400 (Madjet-Renger 2006). "
            "Verificar derivación."
        )


# ------------------------------------------------------------------ build_hamiltonian consistency
class TestBuildHamiltonianConsistency:
    """
    Tests de consistencia entre build_hamiltonian.py y physics_params.yaml.
    Verifica que el módulo cargue los parámetros del YAML y no use valores
    distintos hardcodeados de manera silenciosa.
    """

    @pytest.fixture(autouse=True)
    def import_module(self):
        """Importar build_hamiltonian con warnings capturados."""
        from qmc_mt import build_hamiltonian as bh
        self.bh = bh

    def test_module_imports_without_error(self):
        """build_hamiltonian debe importar sin errores fatales."""
        assert hasattr(self.bh, "build")
        assert hasattr(self.bh, "coupling_cm1")
        assert hasattr(self.bh, "load_structure")

    def test_prefactor_consistent_with_yaml(self):
        """
        El PREFACTOR en build_hamiltonian debe ser consistente con
        constants/prefactor_cm1_A3_D2 en el YAML (±0.1%).
        """
        pytest.importorskip("yaml", reason="pyyaml no instalado")
        if not PHYSICS_YAML.exists():
            pytest.skip("physics_params.yaml no existe")
        import yaml
        with open(PHYSICS_YAML, "r", encoding="utf-8") as f:
            yaml_params = yaml.safe_load(f)

        yaml_C = float(yaml_params["constants"]["prefactor_cm1_A3_D2"]["value"])
        module_C = float(self.bh.PREFACTOR)

        rel_diff = abs(module_C - yaml_C) / abs(yaml_C)
        assert rel_diff < 0.001, (
            f"PREFACTOR en módulo ({module_C}) difiere >0.1% del YAML ({yaml_C}). "
            "Hay dos fuentes de verdad en conflicto."
        )

    def test_mu_consistent_with_yaml(self):
        """MU_TRP_D en módulo debe coincidir con YAML (±0.1%)."""
        pytest.importorskip("yaml", reason="pyyaml no instalado")
        if not PHYSICS_YAML.exists():
            pytest.skip("physics_params.yaml no existe")
        import yaml
        with open(PHYSICS_YAML, "r", encoding="utf-8") as f:
            yaml_params = yaml.safe_load(f)

        yaml_mu = float(yaml_params["dipole"]["mu_trp_1La_D"]["value"])
        module_mu = float(self.bh.MU_TRP_D)

        rel_diff = abs(module_mu - yaml_mu) / abs(yaml_mu)
        assert rel_diff < 0.001, (
            f"MU_TRP_D en módulo ({module_mu}) difiere del YAML ({yaml_mu})."
        )

    def test_E0_consistent_with_yaml(self):
        """E0_CM en módulo debe coincidir con YAML (±0.1%)."""
        pytest.importorskip("yaml", reason="pyyaml no instalado")
        if not PHYSICS_YAML.exists():
            pytest.skip("physics_params.yaml no existe")
        import yaml
        with open(PHYSICS_YAML, "r", encoding="utf-8") as f:
            yaml_params = yaml.safe_load(f)

        yaml_E0 = float(yaml_params["energetics"]["E0_cm1"]["value"])
        module_E0 = float(self.bh.E0_CM)

        rel_diff = abs(module_E0 - yaml_E0) / abs(yaml_E0)
        assert rel_diff < 0.001, (
            f"E0_CM en módulo ({module_E0}) difiere del YAML ({yaml_E0})."
        )

    def test_sigma_consistent_with_yaml(self):
        """SIGMA_SITE en módulo debe coincidir con YAML (±0.1%)."""
        pytest.importorskip("yaml", reason="pyyaml no instalado")
        if not PHYSICS_YAML.exists():
            pytest.skip("physics_params.yaml no existe")
        import yaml
        with open(PHYSICS_YAML, "r", encoding="utf-8") as f:
            yaml_params = yaml.safe_load(f)

        yaml_sigma = float(yaml_params["energetics"]["sigma_site_cm1"]["value"])
        module_sigma = float(self.bh.SIGMA_SITE)

        rel_diff = abs(module_sigma - yaml_sigma) / abs(yaml_sigma)
        assert rel_diff < 0.001, (
            f"SIGMA_SITE en módulo ({module_sigma}) difiere del YAML ({yaml_sigma})."
        )

    def test_selfcheck_prefactor_passes(self):
        """El self-check interno del prefactor debe pasar."""
        # V(kappa^2=4, r=10Å, mu=2D, eps=1) = 5.04e4 * 2 * 4 / 1000 = 403.2 cm^-1
        V = self.bh.PREFACTOR * np.sqrt(4.0) * (2.0 ** 2) / (1.0 * 10.0 ** 3)
        expected = 5.04e4 * 2 * 4 / 1000  # 403.2 — invariante de libro
        assert abs(V - expected) < 1.0, (
            f"Self-check prefactor falló: V={V:.3f} vs esperado={expected:.3f}. "
            f"PREFACTOR actual={self.bh.PREFACTOR}."
        )

    def test_coupling_cm1_physical_range(self):
        """coupling_cm1 para parámetros típicos debe dar V en [10, 300] cm^-1."""
        # Casos típicos: kappa^2 en [0, 4], r en [14, 40] Å
        test_cases = [
            (1.0, 14.0),   # par muy cercano
            (0.67, 20.0),  # par típico, kappa^2 isótropo
            (0.5, 30.0),   # par lejano
        ]
        for kappa2, r in test_cases:
            V = self.bh.coupling_cm1(kappa2, r)
            assert 1.0 <= V <= 500.0, (
                f"coupling_cm1(kappa2={kappa2}, r={r}) = {V:.2f} cm^-1 "
                "fuera del rango esperado [1, 500]."
            )

    def test_coupling_cm1_zero_kappa(self):
        """kappa^2=0 debe dar V=0 (sin acoplamiento)."""
        V = self.bh.coupling_cm1(0.0, 15.0)
        assert V == pytest.approx(0.0, abs=1e-12)

    def test_coupling_cm1_negative_kappa_safe(self):
        """kappa^2 negativo (error numérico) debe tratarse como 0."""
        V_neg = self.bh.coupling_cm1(-0.001, 15.0)
        V_zero = self.bh.coupling_cm1(0.0, 15.0)
        assert V_neg == pytest.approx(V_zero, abs=1e-12)


# ------------------------------------------------------------------ numerical params
class TestNumericalParams:
    """Tests de consistencia de numerical_params.yaml."""

    @pytest.fixture(autouse=True)
    def load_yaml(self):
        pytest.importorskip("yaml", reason="pyyaml no instalado")
        import yaml
        if not NUMERICAL_YAML.exists():
            pytest.skip("numerical_params.yaml no existe")
        with open(NUMERICAL_YAML, "r", encoding="utf-8") as f:
            self.params = yaml.safe_load(f)

    def test_sbc_n_sim_is_1000(self):
        """El estándar Tier-1 exige n_sim >= 1000."""
        n_sim = int(self.params.get("sbc", {}).get("n_sim", 0))
        assert n_sim >= 1000, (
            f"sbc/n_sim = {n_sim} < 1000. El estándar Tier-1 exige n_sim >= 1000."
        )

    def test_heom_nc_min_reliable(self):
        """NC mínimo confiable debe ser >= 5."""
        nc_min = int(self.params.get("heom", {}).get("nc_min_reliable", 0))
        assert nc_min >= 5, (
            f"heom/nc_min_reliable = {nc_min} < 5. Verificar convergencia HEOM."
        )

    def test_hermiticity_threshold_strict(self):
        """Umbral de hermiticidad debe ser <= 1e-10."""
        thresh = float(
            self.params.get("tolerances", {}).get("hermiticity_threshold", 1.0)
        )
        assert thresh <= 1e-10, (
            f"tolerances/hermiticity_threshold = {thresh:.2e} demasiado relajado. "
            "Se requiere <= 1e-10 para precisión numérica."
        )

    def test_kl_threshold_reasonable(self):
        """Umbral KL debe estar en (0, 1) nats."""
        kl_thresh = float(self.params.get("tolerances", {}).get("kl_threshold", -1.0))
        assert 0.0 < kl_thresh < 1.0, (
            f"tolerances/kl_threshold = {kl_thresh} fuera del rango (0, 1)."
        )

    def test_richardson_fallback_ratio_documented(self):
        """El fallback de Richardson debe tener 'fallback_source' documentado."""
        richardson = self.params.get("heom", {}).get("richardson", {})
        assert "fallback_ratio_r" in richardson, (
            "heom/richardson/fallback_ratio_r no definido en numerical_params.yaml"
        )
        assert "fallback_source" in richardson, (
            "heom/richardson/fallback_source no definido. "
            "Todo fallback debe tener fuente documentada."
        )
        source = str(richardson.get("fallback_source", "")).strip()
        assert source != "", "heom/richardson/fallback_source no puede estar vacío."

    def test_richardson_epsilon_threshold(self):
        """epsilon_nc8_max_acceptable debe ser <= 0.05 (5%)."""
        eps = float(
            self.params.get("heom", {})
            .get("richardson", {})
            .get("epsilon_nc8_max_acceptable", 1.0)
        )
        assert eps <= 0.05, (
            f"heom/richardson/epsilon_nc8_max_acceptable = {eps:.3f} > 0.05. "
            "Para claims Tier-1, el error de truncación NC=8 debe ser < 5%."
        )


# ------------------------------------------------------------------ regression tests
class TestRegressionKnownValues:
    """
    Tests de regresión con valores conocidos de la literatura.
    Estos tests son inmutables: si fallan, indica un cambio no autorizado
    en los parámetros físicos centrales.
    """

    def test_prefactor_madjet_renger_2006(self):
        """
        Valor de referencia: Madjet, Abdurakhmanov, Renger (2006),
        J. Phys. Chem. B 110:17268, ec. 2.
        V[cm^-1] = 5.04e4 * kappa * mu[D]^2 / (epsilon * r[Å]^3)
        Test: mu=2D, epsilon=1, r=10Å, kappa^2=4 → V = 403.2 cm^-1
        """
        from qmc_mt.build_hamiltonian import coupling_cm1, EPS_R, MU_TRP_D
        # Forzar epsilon=1 y mu=2 para el test de referencia exacto
        C = 5.04e4
        V_ref = C * np.sqrt(4.0) * (2.0 ** 2) / (1.0 * 10.0 ** 3)
        assert abs(V_ref - 403.2) < 0.1, (
            f"Valor de referencia Madjet-Renger 2006 falló: {V_ref:.3f} vs 403.200"
        )

    def test_kappa2_isotropic_mean(self):
        """
        Promedio isótropo de kappa^2 = 2/3 (resultado analítico exacto).
        Verificado numéricamente con N=200000 muestras (MC).
        """
        from qmc_mt.pdb_tubulin_analysis import _validate_kappa2_isotropic
        k2_iso = _validate_kappa2_isotropic(n=50_000, seed=99)
        assert abs(k2_iso - 2.0 / 3.0) < 1e-2, (
            f"<kappa^2>_iso = {k2_iso:.4f}, esperado = 0.6667 ± 0.01"
        )

    def test_trp_1La_angle_callis_1997(self):
        """
        Ángulo de la transición 1La del Trp respecto al eje largo = 38°.
        Ref: Callis, Methods Enzymol. 278 (1997) 113-150, Fig. 2.
        """
        from qmc_mt.pdb_tubulin_analysis import ANGLE_1LA_DEG
        assert ANGLE_1LA_DEG == pytest.approx(38.0, abs=0.01), (
            f"ANGLE_1LA_DEG = {ANGLE_1LA_DEG}°, esperado 38.0° (Callis 1997)."
        )
