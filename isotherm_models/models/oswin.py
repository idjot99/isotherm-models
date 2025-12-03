"""
Oswin isotherm models for moisture sorption.

This module implements Oswin I and Oswin II models for predicting
moisture sorption isotherms in cellulose-based materials.
"""

import numpy as np
from .base import IsothermModel, NumericType
from ..utils import cardanos_formula_cubic_root, safe_log


class OswinI(IsothermModel):
    """Oswin I isotherm model."""

    @property
    def model_name(self) -> str:
        return "Oswin I"

    @property
    def description(self) -> str:
        return "Oswin I model for moisture sorption"

    def _compute_xi(self, omega: NumericType, theta: NumericType) -> NumericType:
        c = self.params["c"]
        kappa_val = self.kappa(theta)
        omega_fsp_val = self.omega_fsp(theta)
        return (omega / kappa_val / (1 - (omega / omega_fsp_val) ** self.m)) ** (1 / c)

    def compute_water_activity(
        self, omega: NumericType, theta: NumericType
    ) -> NumericType:
        xi = self._compute_xi(omega, theta)
        return xi / (1 + xi)

    def _compute_v(self, aw: NumericType, theta: NumericType) -> NumericType:
        c = self.params["c"]
        kappa_val = self.kappa(theta)
        omega_fsp_val = self.omega_fsp(theta)
        aw_safe = np.clip(aw, 1e-10, 1 - 1e-10)
        return omega_fsp_val / (3 * kappa_val * (aw_safe / (1 - aw_safe)) ** c)

    def compute_moisture_ratio(
        self, aw: NumericType, theta: NumericType
    ) -> NumericType:
        v_val = self._compute_v(aw, theta)
        omega_fsp_val = self.omega_fsp(theta)
        return omega_fsp_val * cardanos_formula_cubic_root(v_val)

    def _compute_h_function(self, aw: NumericType) -> NumericType:
        return 1 - aw


class OswinII(IsothermModel):
    """Oswin II isotherm model with additional parameter b."""

    @property
    def model_name(self) -> str:
        return "Oswin II"

    @property
    def description(self) -> str:
        return "Oswin II model with extended parameterization"

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
        return ((xi ** (1 / b)) / (1 + xi ** (1 / b))) ** b

    def _compute_v(self, aw: NumericType, theta: NumericType) -> NumericType:
        c = self.params["c"]
        b = self.params["b"]
        kappa_val = self.kappa(theta)
        omega_fsp_val = self.omega_fsp(theta)
        aw_safe = np.clip(aw, 1e-10, 1 - 1e-10)
        return omega_fsp_val / (
            3 * kappa_val * ((aw_safe ** (1 / b) / (1 - aw_safe ** (1 / b))) ** b) ** c
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
        return 1 - aw_safe ** (1 / b)
