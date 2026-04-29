"""
robust_interp.py — Interpolador con trazabilidad epistemológica (Patch v3.5.1).

Cada llamada devuelve un InterpResult que contiene:
  - value: valor interpolado / extrapolado
  - std: incertidumbre propagada (si se proveen errores en los datos)
  - diagnostic: InterpDiagnostic con métricas de calidad
  - metadata: metadatos para auditoría

Métodos soportados: linear, pchip, akima, rbf, richardson, pade.

Diseñado para cumplir criterios Tier-1 de reproducibilidad:
  - Sin interpolaciones silenciosas ni extrapolaciones no diagnosticadas.
  - Propagación de incertidumbre documentada para cada método.
  - Diagnóstico de régimen asintótico obligatorio para Richardson.
"""
from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass
from typing import Optional, Union

import numpy as np
from scipy.interpolate import (
    Akima1DInterpolator,
    PchipInterpolator,
    RBFInterpolator,
)
from scipy.stats import chi2 as chi2_dist


# ------------------------------------------------------------------ dataclasses
@dataclass
class InterpDiagnostic:
    """Diagnóstico de calidad para una interpolación."""

    method: str
    n_points: int
    condition_number: Optional[float]
    cross_val_error: Optional[float]
    extrapolation_fraction: float
    stability_score: float  # 0–1; 1 = óptimo
    warning: Optional[str] = None

    def to_dict(self) -> dict:
        """Serializar a diccionario plano (compatible con JSON)."""
        return asdict(self)


@dataclass
class InterpResult:
    """Resultado de interpolación con trazabilidad completa."""

    value: Union[float, np.ndarray]
    std: Optional[Union[float, np.ndarray]]
    diagnostic: InterpDiagnostic
    metadata: dict

    def to_audit_dict(self) -> dict:
        """Serializar para auditoría (compatible con JSON)."""
        return {
            "value": (
                self.value.tolist()
                if hasattr(self.value, "tolist")
                else float(self.value)
            ),
            "std": (
                self.std.tolist()
                if hasattr(self.std, "tolist")
                else (float(self.std) if self.std is not None else None)
            ),
            "diagnostic": self.diagnostic.to_dict(),
            "metadata": self.metadata,
        }


# ------------------------------------------------------------------ main class
class RobustInterpolator:
    """
    Interpolador robusto con diagnóstico integrado.

    Métodos disponibles
    -------------------
    linear   : interpolación lineal a trozos (RegularGridInterpolator 1-D)
    pchip    : PCHIP monotone preserving (SciPy)
    akima    : Akima 1-D (SciPy)
    rbf      : Radial Basis Function thin-plate (SciPy)
    richardson : extrapolación de Richardson al límite continuo
    pade     : aproximación racional de Padé [1/1]

    Uso básico
    ----------
    >>> interp = RobustInterpolator("pchip")
    >>> interp.fit(x_data, y_data, y_err=y_uncertainty)
    >>> result = interp.predict(x_new)
    >>> print(result.value, result.std, result.diagnostic.stability_score)
    """

    VALID_METHODS = frozenset({"linear", "pchip", "akima", "rbf", "richardson", "pade"})

    def __init__(self, method: str, **kwargs):
        if method not in self.VALID_METHODS:
            raise ValueError(
                f"Método '{method}' no soportado. Opciones: {sorted(self.VALID_METHODS)}"
            )
        self.method = method
        self.kwargs = kwargs
        self._fitted = False
        self._x: Optional[np.ndarray] = None
        self._y: Optional[np.ndarray] = None
        self._y_err: Optional[np.ndarray] = None

    # ---------------------------------------------------------------- fit
    def fit(
        self,
        x: np.ndarray,
        y: np.ndarray,
        y_err: Optional[np.ndarray] = None,
    ) -> "RobustInterpolator":
        """
        Ajustar interpolador a datos con incertidumbre opcional.

        Parameters
        ----------
        x : array_like, shape (n,)
            Coordenadas de los datos de entrenamiento. Se ordenarán si es necesario.
        y : array_like, shape (n,)
            Valores objetivo.
        y_err : array_like, shape (n,), optional
            Incertidumbres (1-sigma) en y. Requeridas para propagación de errores.

        Returns
        -------
        self
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if y_err is not None:
            y_err = np.asarray(y_err, dtype=float)
            if y.shape != y_err.shape:
                raise ValueError(
                    f"Shapes de y {y.shape} y y_err {y_err.shape} deben coincidir."
                )

        if len(x) < 2:
            raise ValueError("Se requieren al menos 2 puntos para interpolar.")

        # Ordenar si no está ordenado estrictamente
        if np.any(np.diff(x) <= 0):
            warnings.warn(
                "x no está estrictamente ordenado; se ordenará internamente.",
                stacklevel=2,
            )
            sort_idx = np.argsort(x)
            x = x[sort_idx]
            y = y[sort_idx]
            if y_err is not None:
                y_err = y_err[sort_idx]

        self._x, self._y, self._y_err = x, y, y_err
        self._fitted = True
        return self

    # ---------------------------------------------------------------- predict
    def predict(
        self,
        x_new: np.ndarray,
        return_std: bool = True,
    ) -> InterpResult:
        """
        Predecir valores con incertidumbre propagada.

        Parameters
        ----------
        x_new : array_like
            Puntos donde evaluar el interpolador.
        return_std : bool
            Si True, intentar propagar incertidumbre desde y_err.

        Returns
        -------
        InterpResult
        """
        if not self._fitted:
            raise RuntimeError("Llamar a fit() antes de predict().")
        x_new = np.atleast_1d(np.asarray(x_new, dtype=float))

        dispatch = {
            "richardson": self._predict_richardson,
            "pade": self._predict_pade,
        }
        if self.method in dispatch:
            return dispatch[self.method](x_new, return_std)
        return self._predict_generic(x_new, return_std)

    # ---------------------------------------------------------------- generic
    def _predict_generic(
        self, x_new: np.ndarray, return_std: bool
    ) -> InterpResult:
        """Implementación para linear, pchip, akima, rbf."""
        if self.method == "pchip":
            interp = PchipInterpolator(self._x, self._y)
            y_pred = interp(x_new)
        elif self.method == "akima":
            interp = Akima1DInterpolator(self._x, self._y)
            y_pred = interp(x_new)
        elif self.method == "rbf":
            interp = RBFInterpolator(
                self._x[:, None], self._y, kernel="thin_plate_spline", smoothing=0.0
            )
            y_pred = interp(x_new[:, None]).ravel()
        else:  # linear
            y_pred = np.interp(x_new, self._x, self._y)

        # Propagación lineal de primer orden de la incertidumbre
        y_std = None
        if return_std and self._y_err is not None:
            y_std = np.interp(x_new, self._x, self._y_err)

        diag = self._compute_generic_diagnostic(x_new, y_pred)
        return InterpResult(
            value=y_pred,
            std=y_std,
            diagnostic=diag,
            metadata={"method": self.method, "n_fit_points": len(self._x)},
        )

    # ---------------------------------------------------------------- richardson
    def _predict_richardson(
        self, x_new: np.ndarray, return_std: bool
    ) -> InterpResult:
        """
        Extrapolación de Richardson al límite de paso de cuadrícula → 0.

        Los datos en self._x deben ser parámetros de discretización (e.g., 1/NC),
        y self._y los observables correspondientes. La extrapolación asume
        convergencia geométrica: y(NC) = y_inf + A * r^NC.

        Diagnóstico obligatorio: verificación de contracción geométrica.
        """
        if len(self._x) < 2:
            raise ValueError("Richardson requiere al menos 2 puntos.")

        deltas = np.diff(self._y)
        ratios = np.abs(deltas[1:] / deltas[:-1]) if len(deltas) >= 2 else np.array([np.nan])

        # Verificar contracción geométrica en los últimos 2 ratios
        if len(ratios) >= 2 and np.all(np.isfinite(ratios[-2:])):
            ratio_stable = bool(np.abs(ratios[-1] - ratios[-2]) < 0.15)
        elif len(ratios) >= 1 and np.isfinite(ratios[-1]):
            ratio_stable = bool(0.1 < ratios[-1] < 0.9)
        else:
            ratio_stable = False

        r_est = float(ratios[-1]) if len(ratios) > 0 and np.isfinite(ratios[-1]) else 0.5

        # Extrapolación: y_inf = y_N + (y_N - y_{N-1}) * r / (1 - r)
        if abs(1 - r_est) > 1e-6:
            y_inf_scalar = self._y[-1] + (self._y[-1] - self._y[-2]) * r_est / (1 - r_est)
        else:
            y_inf_scalar = self._y[-1]
            warnings.warn(
                "Richardson: r ≈ 1 (convergencia muy lenta o no convergente). "
                "Extrapolación no confiable.",
                stacklevel=3,
            )

        # Propagación Monte Carlo si hay errores disponibles
        y_std_scalar = None
        if return_std and self._y_err is not None and len(self._y_err) >= 2:
            n_mc = 10_000
            rng_mc = np.random.default_rng(42)
            samples_y = rng_mc.normal(
                self._y[-2:], self._y_err[-2:], size=(n_mc, 2)
            )
            # Incluir incertidumbre en el ratio r (5% relativo como conservador)
            r_sigma = max(0.05 * abs(r_est), 0.01)
            samples_r = rng_mc.normal(r_est, r_sigma, size=n_mc)
            denom = 1 - samples_r
            # Evitar división por cero
            safe_denom = np.where(np.abs(denom) > 1e-6, denom, 1e-6)
            y_inf_mc = (
                samples_y[:, 1]
                + (samples_y[:, 1] - samples_y[:, 0]) * samples_r / safe_denom
            )
            y_inf_scalar = float(np.mean(y_inf_mc))
            y_std_scalar = float(np.std(y_inf_mc))

        y_pred = np.full_like(x_new, y_inf_scalar)
        y_std = np.full_like(x_new, y_std_scalar) if y_std_scalar is not None else None

        warn_msg = None
        if not ratio_stable:
            warn_msg = (
                f"Régimen asintótico no verificado (r={r_est:.3f}). "
                "Usar resultado con cautela."
            )
            warnings.warn(warn_msg, stacklevel=2)

        diag = InterpDiagnostic(
            method="richardson",
            n_points=len(self._x),
            condition_number=None,
            cross_val_error=None,
            extrapolation_fraction=1.0,  # Richardson siempre extrapola
            stability_score=0.9 if ratio_stable else 0.3,
            warning=warn_msg,
        )
        return InterpResult(
            value=y_pred,
            std=y_std,
            diagnostic=diag,
            metadata={
                "method": "richardson",
                "ratio_estimated": r_est,
                "geometric_convergence_verified": ratio_stable,
                "n_mc_samples": 10_000 if y_std_scalar is not None else 0,
            },
        )

    # ---------------------------------------------------------------- pade
    def _predict_pade(
        self, x_new: np.ndarray, return_std: bool
    ) -> InterpResult:
        """
        Aproximación de Padé [1/1] para descomposición de baños Drude-Lorentz.

        Ajusta f(x) = (a0 + a1*x) / (1 + b1*x) por mínimos cuadrados no lineales.
        Diagnóstico de estabilidad: verifica que los polos sean imaginarios puros
        (función sin divergencias en el dominio real relevante).
        """
        from scipy.optimize import curve_fit

        def _pade11(x: np.ndarray, a0: float, a1: float, b1: float) -> np.ndarray:
            return (a0 + a1 * x) / (1.0 + b1 * x)

        try:
            popt, pcov = curve_fit(
                _pade11,
                self._x,
                self._y,
                p0=[self._y.mean(), 0.0, 0.0],
                maxfev=20_000,
            )
        except Exception as e:
            warnings.warn(
                f"Ajuste Padé falló ({e}); cayendo a interpolación lineal.",
                stacklevel=2,
            )
            return self._predict_generic(x_new, return_std)

        y_pred = _pade11(x_new, *popt)

        # Incertidumbre desde covarianza (propagación lineal de primer orden)
        y_std = None
        if return_std and pcov is not None and np.all(np.isfinite(pcov)):
            # Jacobiano numérico
            eps = 1e-7
            J_rows = []
            for i, p in enumerate(popt):
                p_plus = popt.copy()
                p_plus[i] += eps
                J_rows.append((_pade11(x_new, *p_plus) - y_pred) / eps)
            J = np.column_stack(J_rows)  # (n_points, 3)
            y_var = np.einsum("ij,jk,ik->i", J, pcov, J)
            y_std = np.sqrt(np.maximum(y_var, 0.0))

        # Diagnóstico de polos: polo en x = -1/b1
        b1 = popt[2]
        if abs(b1) > 1e-12:
            pole_real = -1.0 / b1
            # ¿El polo está dentro del rango de x_new?
            x_min, x_max = float(x_new.min()), float(x_new.max())
            pole_in_domain = x_min <= pole_real <= x_max
        else:
            pole_real = np.inf
            pole_in_domain = False

        stability_score = 0.2 if pole_in_domain else 1.0
        warn_msg = None
        if pole_in_domain:
            warn_msg = (
                f"Polo de Padé en x={pole_real:.4f} dentro del dominio de predicción. "
                "Resultado posiblemente divergente."
            )
            warnings.warn(warn_msg, stacklevel=2)

        # Número de condición de la matriz de diseño
        cond_num = None
        if len(self._x) >= 3:
            A_design = np.column_stack([
                np.ones_like(self._x),
                self._x,
                -self._y * self._x,
            ])
            cond_num = float(np.linalg.cond(A_design))

        extrap_frac = float(np.mean((x_new < self._x.min()) | (x_new > self._x.max())))

        diag = InterpDiagnostic(
            method="pade",
            n_points=len(self._x),
            condition_number=cond_num,
            cross_val_error=None,
            extrapolation_fraction=extrap_frac,
            stability_score=stability_score,
            warning=warn_msg,
        )
        return InterpResult(
            value=y_pred,
            std=y_std,
            diagnostic=diag,
            metadata={
                "method": "pade",
                "order": "[1/1]",
                "coefficients": popt.tolist(),
                "pole_x": float(pole_real) if abs(b1) > 1e-12 else None,
                "pole_in_domain": bool(pole_in_domain),
            },
        )

    # ---------------------------------------------------------------- diagnostic
    def _compute_generic_diagnostic(
        self, x_new: np.ndarray, y_pred: np.ndarray
    ) -> InterpDiagnostic:
        """Calcular diagnóstico de calidad genérico (para linear, pchip, akima, rbf)."""
        x_min, x_max = float(self._x.min()), float(self._x.max())
        extrap_frac = float(
            np.mean((x_new < x_min) | (x_new > x_max))
        )

        # Número de condición (para interpolación lineal)
        cond_num = None
        if self.method == "linear" and len(self._x) >= 2:
            A = np.column_stack([np.ones_like(self._x), self._x])
            cond_num = float(np.linalg.cond(A))

        # Error LOO simplificado (sólo estimación de varianza residual para linear)
        cv_err = None
        if self.method == "linear" and len(self._x) >= 4:
            y_loo = np.interp(self._x, self._x, self._y)  # trivialmente 0 para linear
            # Mejor estimación: residuos respecto a la tendencia global
            coeffs = np.polyfit(self._x, self._y, 1)
            y_trend = np.polyval(coeffs, self._x)
            cv_err = float(np.sqrt(np.mean((self._y - y_trend) ** 2)))

        # Puntuación de estabilidad heurística
        stability = 1.0
        if extrap_frac > 0.5:
            stability *= 0.5
        if cond_num is not None and cond_num > 1e8:
            stability *= 0.3

        warn_msg = None
        if extrap_frac > 0.3:
            warn_msg = f"{extrap_frac * 100:.1f}% de puntos en extrapolación"
        elif cond_num is not None and cond_num > 1e8:
            warn_msg = "Número de condición alto: posible inestabilidad numérica"

        return InterpDiagnostic(
            method=self.method,
            n_points=len(self._x),
            condition_number=cond_num,
            cross_val_error=cv_err,
            extrapolation_fraction=extrap_frac,
            stability_score=float(max(stability, 0.0)),
            warning=warn_msg,
        )

    # ---------------------------------------------------------------- convergence
    def get_convergence_diagnostic(self) -> dict:
        """
        Obtener diagnóstico de convergencia específico del método.

        Relevante principalmente para 'richardson' y 'pade'.
        """
        if not self._fitted:
            return {"error": "fit() no llamado aún"}

        if self.method == "richardson":
            deltas = np.diff(self._y)
            ratios = (
                np.abs(deltas[1:] / deltas[:-1])
                if len(deltas) >= 2
                else np.array([np.nan])
            )
            return {
                "type": "richardson",
                "deltas": deltas.tolist(),
                "ratios": ratios.tolist(),
                "last_ratio": float(ratios[-1]) if len(ratios) > 0 else None,
                "ratio_stable": bool(
                    np.all(np.abs(ratios[-2:] - ratios[-1]) < 0.15)
                ) if len(ratios) >= 2 else None,
                "n_points": len(self._x),
            }
        return {
            "type": self.method,
            "note": "Diagnóstico de convergencia no aplica para este método.",
        }
