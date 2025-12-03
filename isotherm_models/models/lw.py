"""
Leitner-Wolkoff (LW) isotherm models for moisture sorption.

This module implements LW I and LW II models that use the Lambert W function
for predicting moisture sorption isotherms in cellulose-based materials.
"""

import numpy as np
from .base import IsothermModel, NumericType
from ..utils import lambert_w, cardanos_formula_cubic_root, safe_log


class LWI(IsothermModel):
    """Leitner-Wolkoff I isotherm model."""

    @property
    def model_name(self) -> str:
        return "LW I"

    @property
    def description(self) -> str:
        return "Leitner-Wolkoff I model using Lambert W function"

    def _compute_xi(self, omega: NumericType, theta: NumericType) -> NumericType:
        c = self.params["c"]
        kappa_val = self.kappa(theta)
        omega_fsp_val = self.omega_fsp(theta)
        return (omega / kappa_val / (1 - (omega / omega_fsp_val) ** self.m)) ** (1 / c)

    def compute_water_activity(
        self, omega: NumericType, theta: NumericType
    ) -> NumericType:
        xi = self._compute_xi(omega, theta)
        return 1 - np.exp(-lambert_w(xi))

    def _compute_v(self, aw: NumericType, theta: NumericType) -> NumericType:
        c = self.params["c"]
        kappa_val = self.kappa(theta)
        omega_fsp_val = self.omega_fsp(theta)
        aw_safe = np.clip(aw, 1e-10, 1 - 1e-10)
        return omega_fsp_val / (
            3 * kappa_val * (-safe_log(1 - aw_safe) / (1 - aw_safe)) ** c
        )

    def compute_moisture_ratio(
        self, aw: NumericType, theta: NumericType
    ) -> NumericType:
        v_val = self._compute_v(aw, theta)
        omega_fsp_val = self.omega_fsp(theta)
        return omega_fsp_val * cardanos_formula_cubic_root(v_val)

    def _compute_h_function(self, aw: NumericType) -> NumericType:
        aw_safe = np.clip(aw, 1e-10, 1 - 1e-10)
        log_term = safe_log(1 - aw_safe)
        return -(1 - aw_safe) / aw_safe * log_term / (1 - log_term)


class LWII(IsothermModel):
    """Leitner-Wolkoff II isotherm model with additional parameter b."""

    @property
    def model_name(self) -> str:
        return "LW II"

    @property
    def description(self) -> str:
        return "Leitner-Wolkoff II model with extended parameterization"

    def _compute_xi(self, omega: NumericType, theta: NumericType) -> NumericType:
        c = self.params["c"]
        kappa_val = self.kappa(theta)
        omega_fsp_val = self.omega_fsp(theta)
        return (omega / kappa_val / (1 - (omega / omega_fsp_val) ** self.m)) ** (1 / c)

    def compute_water_activity(
        self, omega: NumericType, theta: NumericType
    ) -> NumericType:
        b = self.params["b"]
        xi = self._compute_xi(omega, theta)
        return (1 - np.exp(-lambert_w(xi ** (1 / b)))) ** b

    def _compute_v(self, aw: NumericType, theta: NumericType) -> NumericType:
        c = self.params["c"]
        b = self.params["b"]
        kappa_val = self.kappa(theta)
        omega_fsp_val = self.omega_fsp(theta)
        aw_safe = np.clip(aw, 1e-10, 1 - 1e-10)
        aw_b = aw_safe ** (1 / b)
        aw_b_safe = np.clip(aw_b, 1e-10, 1 - 1e-10)
        return omega_fsp_val / (
            3
            * kappa_val
            * ((-safe_log(1 - aw_b_safe) / (1 - aw_b_safe)) ** b) ** c
        )

    def compute_moisture_ratio(
        self, aw: NumericType, theta: NumericType
    ) -> NumericType:
        v_val = self._compute_v(aw, theta)
        omega_fsp_val = self.omega_fsp(theta)
        return omega_fsp_val * cardanos_formula_cubic_root(v_val)

    def _compute_h_function(self, aw: NumericType) -> NumericType:
        b = self.params["b"]
        aw_safe = np.clip(aw, 1e-10, 1 - 1e-10)
        aw_b = aw_safe ** (1 / b)
        aw_b_safe = np.clip(aw_b, 1e-10, 1 - 1e-10)
        log_term = safe_log(1 - aw_b_safe)
        return -(1 - aw_b_safe) / aw_b_safe * log_term / (1 - log_term)
