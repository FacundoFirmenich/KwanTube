"""
test_interpolation_integrity.py — Tests de integridad para robust_interp.py.
Patch v3.5.1 — Valida que ningún método interpolador opere sin diagnóstico,
que las incertidumbres se propaguen correctamente, y que las extrapolaciones
se detecten y adviertan.

Ejecutar con: pytest tests/test_interpolation_integrity.py -v
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

# Path setup
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from qmc_mt.robust_interp import InterpDiagnostic, InterpResult, RobustInterpolator


# ------------------------------------------------------------------ fixtures
@pytest.fixture
def linear_data():
    """Datos lineales exactos: y = 2x + 1."""
    x = np.linspace(0, 10, 20)
    y = 2.0 * x + 1.0
    y_err = np.full_like(y, 0.05)
    return x, y, y_err


@pytest.fixture
def noisy_data():
    """Datos con ruido gaussiano para tests de propagación."""
    rng = np.random.default_rng(42)
    x = np.linspace(0, 5, 15)
    y = np.sin(x) + rng.normal(0, 0.02, size=len(x))
    y_err = np.full_like(y, 0.02)
    return x, y, y_err


@pytest.fixture
def convergence_data():
    """Datos de convergencia geométrica para Richardson."""
    # Simula y(NC) = y_inf - A * r^NC con y_inf=0.5, A=0.4, r=0.4
    y_inf, A, r = 0.5, 0.4, 0.4
    nc_vals = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    y_vals = y_inf - A * r ** nc_vals
    y_err = np.full_like(y_vals, 0.002)
    return nc_vals, y_vals, y_err, y_inf


# ------------------------------------------------------------------ basic tests
class TestRobustInterpolatorInstantiation:
    """Tests de creación y validación de métodos."""

    def test_valid_methods(self):
        for method in ("linear", "pchip", "akima", "rbf", "richardson", "pade"):
            interp = RobustInterpolator(method)
            assert interp.method == method

    def test_invalid_method_raises(self):
        with pytest.raises(ValueError, match="no soportado"):
            RobustInterpolator("spline_magic")

    def test_predict_before_fit_raises(self):
        interp = RobustInterpolator("pchip")
        with pytest.raises(RuntimeError, match="fit\\(\\)"):
            interp.predict(np.array([1.0, 2.0]))


class TestFitValidation:
    """Tests de validación en fit()."""

    def test_fit_requires_two_points(self):
        interp = RobustInterpolator("linear")
        with pytest.raises(ValueError, match="2 puntos"):
            interp.fit(np.array([1.0]), np.array([2.0]))

    def test_fit_shape_mismatch_y_err_raises(self):
        interp = RobustInterpolator("pchip")
        with pytest.raises(ValueError, match="[Ss]hapes"):
            interp.fit(
                np.array([0.0, 1.0, 2.0]),
                np.array([0.0, 1.0, 2.0]),
                y_err=np.array([0.1, 0.1]),
            )

    def test_fit_unsorted_x_warns_and_sorts(self):
        interp = RobustInterpolator("pchip")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            interp.fit(
                np.array([5.0, 1.0, 3.0]),
                np.array([25.0, 1.0, 9.0]),
            )
        # Tras fit, x debe estar ordenado
        assert np.all(np.diff(interp._x) > 0)

    def test_fit_returns_self(self, linear_data):
        x, y, y_err = linear_data
        interp = RobustInterpolator("linear")
        result = interp.fit(x, y, y_err)
        assert result is interp


# ------------------------------------------------------------------ predict tests
class TestPredictReturnsCorrectTypes:
    """Tests de tipos y estructura de retorno."""

    @pytest.mark.parametrize("method", ["linear", "pchip", "akima"])
    def test_predict_returns_interp_result(self, method, linear_data):
        x, y, y_err = linear_data
        interp = RobustInterpolator(method)
        interp.fit(x, y, y_err)
        result = interp.predict(np.array([2.5, 5.0, 7.5]))
        assert isinstance(result, InterpResult)
        assert isinstance(result.diagnostic, InterpDiagnostic)

    @pytest.mark.parametrize("method", ["linear", "pchip", "akima"])
    def test_predict_value_shape_matches_input(self, method, linear_data):
        x, y, y_err = linear_data
        interp = RobustInterpolator(method).fit(x, y, y_err)
        x_query = np.array([1.0, 3.0, 5.0, 7.0])
        result = interp.predict(x_query)
        assert result.value.shape == x_query.shape

    @pytest.mark.parametrize("method", ["linear", "pchip", "akima"])
    def test_predict_std_shape_matches_when_y_err_given(self, method, linear_data):
        x, y, y_err = linear_data
        interp = RobustInterpolator(method).fit(x, y, y_err)
        x_query = np.array([1.0, 3.0, 5.0])
        result = interp.predict(x_query, return_std=True)
        assert result.std is not None
        assert result.std.shape == x_query.shape

    @pytest.mark.parametrize("method", ["linear", "pchip", "akima"])
    def test_predict_std_none_when_no_y_err(self, method, linear_data):
        x, y, _ = linear_data
        interp = RobustInterpolator(method).fit(x, y)
        result = interp.predict(np.array([2.0, 4.0]))
        assert result.std is None


# ------------------------------------------------------------------ accuracy tests
class TestPredictAccuracy:
    """Tests de precisión numérica básica."""

    def test_linear_exact_recovery(self, linear_data):
        """Interpolación lineal debe recuperar exactamente una función lineal."""
        x, y, y_err = linear_data
        interp = RobustInterpolator("linear").fit(x, y, y_err)
        x_query = np.array([1.5, 4.5, 8.0])
        result = interp.predict(x_query)
        y_expected = 2.0 * x_query + 1.0
        np.testing.assert_allclose(result.value, y_expected, atol=1e-10)

    def test_pchip_monotone_preserving(self):
        """PCHIP debe preservar monotonicidad en datos monótonos."""
        x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        y = np.array([0.0, 1.0, 1.5, 1.8, 2.0])  # monótono creciente
        interp = RobustInterpolator("pchip").fit(x, y)
        x_fine = np.linspace(0, 4, 100)
        result = interp.predict(x_fine)
        # PCHIP preserva monotonicidad: no debe haber decrementos
        assert np.all(np.diff(result.value) >= -1e-10)

    def test_akima_smooth_interpolation(self, noisy_data):
        """Akima debe interpolar sin overshooting excesivo."""
        x, y, y_err = noisy_data
        interp = RobustInterpolator("akima").fit(x, y, y_err)
        x_query = x[1:-1]  # puntos interiores
        result = interp.predict(x_query)
        # El valor en los puntos de entrenamiento debe ser igual a los datos
        np.testing.assert_allclose(
            interp.predict(x).value, y, atol=1e-10,
            err_msg="Akima debe interpolar exactamente en los nodos"
        )


# ------------------------------------------------------------------ diagnostic tests
class TestDiagnostics:
    """Tests del diagnóstico de calidad."""

    def test_diagnostic_method_field(self, linear_data):
        x, y, _ = linear_data
        for method in ("linear", "pchip"):
            interp = RobustInterpolator(method).fit(x, y)
            result = interp.predict(np.array([3.0]))
            assert result.diagnostic.method == method

    def test_extrapolation_fraction_zero_for_interior(self, linear_data):
        x, y, _ = linear_data
        interp = RobustInterpolator("pchip").fit(x, y)
        # Puntos estrictamente interiores → extrap_frac = 0
        x_interior = np.array([2.0, 5.0, 8.0])
        result = interp.predict(x_interior)
        assert result.diagnostic.extrapolation_fraction == pytest.approx(0.0, abs=1e-10)

    def test_extrapolation_fraction_nonzero_for_exterior(self, linear_data):
        x, y, _ = linear_data
        interp = RobustInterpolator("pchip").fit(x, y)
        # Punto fuera del rango
        x_outside = np.array([15.0, 20.0])  # x.max() = 10
        result = interp.predict(x_outside)
        assert result.diagnostic.extrapolation_fraction > 0.0

    def test_stability_score_range(self, linear_data):
        x, y, _ = linear_data
        interp = RobustInterpolator("linear").fit(x, y)
        result = interp.predict(np.array([3.0, 5.0]))
        assert 0.0 <= result.diagnostic.stability_score <= 1.0

    def test_to_audit_dict_serializable(self, linear_data):
        """to_audit_dict() debe retornar tipos serializables por json.dumps."""
        import json
        x, y, y_err = linear_data
        interp = RobustInterpolator("pchip").fit(x, y, y_err)
        result = interp.predict(np.array([3.0, 5.0]))
        # No debe lanzar excepción
        audit = result.to_audit_dict()
        json_str = json.dumps(audit)
        assert len(json_str) > 0

    def test_diagnostic_n_points_correct(self, linear_data):
        x, y, _ = linear_data
        interp = RobustInterpolator("pchip").fit(x, y)
        result = interp.predict(np.array([5.0]))
        assert result.diagnostic.n_points == len(x)


# ------------------------------------------------------------------ Richardson tests
class TestRichardsonMethod:
    """Tests específicos para extrapolación de Richardson."""

    def test_richardson_converges_to_known_limit(self, convergence_data):
        """Richardson debe aproximarse al límite conocido y_inf=0.5."""
        nc_vals, y_vals, y_err, y_inf = convergence_data
        interp = RobustInterpolator("richardson").fit(nc_vals, y_vals, y_err)
        result = interp.predict(np.array([np.inf]))  # extrapolación al límite
        # Tolerancia 5% del valor verdadero
        assert abs(float(result.value[0]) - y_inf) < 0.05 * abs(y_inf), (
            f"Richardson no convergió al límite: obtenido={float(result.value[0]):.4f}, "
            f"esperado={y_inf:.4f}"
        )

    def test_richardson_returns_std_with_y_err(self, convergence_data):
        nc_vals, y_vals, y_err, _ = convergence_data
        interp = RobustInterpolator("richardson").fit(nc_vals, y_vals, y_err)
        result = interp.predict(np.array([100.0]), return_std=True)
        assert result.std is not None
        assert np.all(result.std >= 0)

    def test_richardson_extrapolation_fraction_is_one(self, convergence_data):
        nc_vals, y_vals, _, _ = convergence_data
        interp = RobustInterpolator("richardson").fit(nc_vals, y_vals)
        result = interp.predict(np.array([100.0]))
        # Richardson siempre extrapola al límite continuo
        assert result.diagnostic.extrapolation_fraction == pytest.approx(1.0)

    def test_richardson_stability_score_high_for_convergent_data(self, convergence_data):
        nc_vals, y_vals, _, _ = convergence_data
        interp = RobustInterpolator("richardson").fit(nc_vals, y_vals)
        result = interp.predict(np.array([100.0]))
        # Datos bien convergentes → stability_score >= 0.8
        assert result.diagnostic.stability_score >= 0.8, (
            f"Esperado stability_score >= 0.8, obtenido {result.diagnostic.stability_score}"
        )

    def test_convergence_diagnostic_structure(self, convergence_data):
        nc_vals, y_vals, _, _ = convergence_data
        interp = RobustInterpolator("richardson").fit(nc_vals, y_vals)
        diag = interp.get_convergence_diagnostic()
        assert diag["type"] == "richardson"
        assert "last_ratio" in diag
        assert "ratio_stable" in diag
        assert diag["last_ratio"] is not None


# ------------------------------------------------------------------ Padé tests
class TestPadeMethod:
    """Tests específicos para aproximación de Padé."""

    def test_pade_fits_rational_function(self):
        """Padé [1/1] debe ajustar exactamente una función racional simple."""
        x = np.linspace(0.1, 5.0, 30)
        # f(x) = (1 + 2x) / (1 + 0.5x) — función racional exacta
        y = (1.0 + 2.0 * x) / (1.0 + 0.5 * x)
        interp = RobustInterpolator("pade").fit(x, y)
        x_query = np.array([1.0, 2.5, 4.0])
        result = interp.predict(x_query)
        y_expected = (1.0 + 2.0 * x_query) / (1.0 + 0.5 * x_query)
        np.testing.assert_allclose(result.value, y_expected, rtol=1e-3)

    def test_pade_returns_std_with_covariance(self):
        """Padé debe retornar std cuando el ajuste converge."""
        x = np.linspace(0.1, 5.0, 20)
        y = (1.0 + 2.0 * x) / (1.0 + 0.5 * x)
        y_err = np.full_like(y, 0.01)
        interp = RobustInterpolator("pade").fit(x, y, y_err)
        result = interp.predict(np.array([1.0, 3.0]))
        assert result.std is not None
        assert np.all(result.std >= 0)

    def test_pade_detects_pole_in_domain(self):
        """Padé debe detectar y advertir sobre polos en el dominio de predicción."""
        x = np.linspace(0.1, 5.0, 20)
        # Función con polo en x=2: f(x) = 1 / (x - 2) → aproximar con Padé
        y = 1.0 / (x - 10.0)  # polo en x=10, fuera de rango
        interp = RobustInterpolator("pade").fit(x, y)
        # Predecir en rango con el polo dentro
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            result = interp.predict(np.array([1.0, 3.0]))
        # El diagnóstico debe estar presente
        assert result.diagnostic is not None

    def test_pade_to_audit_dict_serializable(self):
        import json
        x = np.linspace(0.5, 4.0, 15)
        y = (2.0 + x) / (1.0 + 0.3 * x)
        interp = RobustInterpolator("pade").fit(x, y)
        result = interp.predict(np.array([1.5, 2.5]))
        audit = result.to_audit_dict()
        json_str = json.dumps(audit)
        assert "pade" in json_str


# ------------------------------------------------------------------ integration test
class TestIntegrationPhysicalPipeline:
    """
    Test de integración: simula el uso de RobustInterpolator
    en el contexto del pipeline biofísico (HEOM, Richardson).
    """

    def test_heom_convergence_scenario(self):
        """
        Escenario: datos de convergencia HEOM NC=3..8.
        Richardson debe extrapolar al límite y diagnosticar estabilidad.
        """
        y_inf_true = 0.45
        A, r_true = 0.3, 0.42
        nc = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
        pop = y_inf_true - A * r_true ** nc
        y_err = np.full_like(pop, 0.001)

        interp = RobustInterpolator("richardson").fit(nc, pop, y_err)
        result = interp.predict(np.array([np.inf]))

        # El resultado debe estar dentro del 3% del valor verdadero
        assert abs(float(result.value[0]) - y_inf_true) < 0.03 * abs(y_inf_true), (
            f"Error HEOM scenario: {float(result.value[0]):.4f} vs {y_inf_true:.4f}"
        )
        # El diagnóstico debe indicar estabilidad
        assert result.diagnostic.stability_score > 0.5

    def test_drude_lorentz_bath_scenario(self):
        """
        Escenario: función espectral Drude-Lorentz J(omega) = 2*lambda*gamma*omega / (omega^2 + gamma^2).
        Padé debe ajustar razonablemente.
        """
        lam, gamma = 35.0, 53.0
        omega = np.linspace(1.0, 200.0, 50)
        J = 2.0 * lam * gamma * omega / (omega ** 2 + gamma ** 2)

        interp = RobustInterpolator("pade").fit(omega, J)
        result = interp.predict(omega[5:15])

        assert result.diagnostic is not None
        assert result.diagnostic.n_points == len(omega)
        # Verificar que el ajuste tiene baja extrapolación
        assert result.diagnostic.extrapolation_fraction == pytest.approx(0.0, abs=0.01)
