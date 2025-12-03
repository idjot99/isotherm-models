"""
GAB (Guggenheim-Anderson-de Boer) isotherm models.

This module implements GAB models for predicting moisture sorption isotherms.
The GAB model is widely used for food and biological materials.
"""

import numpy as np
from typing import Union
from .base import IsothermModel, NumericType
from ..utils import R_SPECIFIC_GAS


class GABI(IsothermModel):
    """
    GAB I isotherm model.

    This variant uses temperature-dependent C and K parameters.
    """

    @property
    def model_name(self) -> str:
        return "GAB I"

    @property
    def description(self) -> str:
        return "GAB I model with temperature-dependent parameters"

    def K(self, theta: NumericType) -> NumericType:
        """Compute temperature-dependent K parameter."""
        K0 = self.params["K0"]
        qK = self.params["qK"]
        return K0 * np.exp(qK / R_SPECIFIC_GAS / theta)

    def C(self, theta: NumericType) -> NumericType:
        """Compute temperature-dependent C parameter."""
        C0 = self.params["C0"]
        qC = self.params["qC"]
        return C0 * np.exp(qC / R_SPECIFIC_GAS / theta)

    def compute_water_activity(
        self, omega: NumericType, theta: NumericType
    ) -> NumericType:
        """Compute water activity for GAB I model."""
        C_val = self.C(theta)
        K_val = self.K(theta)
        omega_fsp_val = self.omega_fsp(theta)
        y = omega / omega_fsp_val

        # Solve quadratic equation for aw
        A = (C_val - 2) * K_val * (y - 1) + (C_val - 1) * K_val ** 2 - 1
        B = 2 * (C_val - 1) * K_val ** 2 * y

        sqrt_term = np.sqrt(A ** 2 + B * y)
        numerator = sqrt_term + A

        denominator = B

        # Avoid division by zero
        denominator_safe = np.where(np.abs(denominator) < 1e-10, 1e-10, denominator)

        return numerator / denominator_safe

    def compute_moisture_ratio(
        self, aw: NumericType, theta: NumericType
    ) -> NumericType:
        """Compute moisture ratio from water activity (inverse function)."""
        C_val = self.C(theta)
        K_val = self.K(theta)
        omega_fsp_val = self.omega_fsp(theta)

        return (
            omega_fsp_val
            * aw
            * (1 - K_val)
            / (1 - K_val * aw)
            * (1 - K_val + C_val * K_val)
            / (1 - K_val * aw + C_val * K_val * aw)
        )

    def compute_net_heat_of_sorption(
        self, omega: NumericType, theta: NumericType
    ) -> NumericType:
        """Compute net heat of sorption for GAB I model."""
        C0 = self.params["C0"]
        qC = self.params["qC"]
        K0 = self.params["K0"]
        qK = self.params["qK"]

        C_val = self.C(theta)
        K_val = self.K(theta)
        omega_fsp_val = self.omega_fsp(theta)

        y = omega / omega_fsp_val

        A = (C_val - 2) * K_val * (y - 1) + (C_val - 1) * K_val ** 2 - 1
        B = 2 * (C_val - 1) * K_val ** 2 * y

        dy = y / omega_fsp_val * self.omega_ref / self.theta1

        rdA = (
            -C_val * K_val * qC * (K_val + y - 1)
            - K_val * qK * ((C_val - 2) * (y - 1) + 2 * (C_val - 1) * K_val)
            + (C_val - 2) * K_val * dy * R_SPECIFIC_GAS * theta ** 2
        )

        rdB = (
            -2 * C_val * qC * K_val ** 2 * y
            - 4 * (C_val - 1) * K_val ** 2 * qK * y
            + 2 * (C_val - 1) * K_val ** 2 * dy * R_SPECIFIC_GAS * theta ** 2
        )

        sqrt_term = np.sqrt(A ** 2 + 2 * B * y)

        Hiso = (
            (A * rdA + rdB * y + B * dy * R_SPECIFIC_GAS * theta ** 2)
            / ((sqrt_term + A) * sqrt_term)
            + rdA / (sqrt_term + A)
            - rdB / B
        )

        return Hiso

    def _compute_h_function(self, aw: NumericType) -> NumericType:
        """Not used for GAB I model (uses custom heat calculation)."""
        return np.ones_like(aw)  # Placeholder, not actually used


class GAB(IsothermModel):
    """
    Standard GAB (Guggenheim-Anderson-de Boer) isotherm model.

    This is the classical formulation used widely in food science.

    References
    ----------
    Furmaniak, Terzyk and Gauden (2007). Journal of Food Engineering 82:528-535
    """

    @property
    def model_name(self) -> str:
        return "GAB"

    @property
    def description(self) -> str:
        return "Standard GAB model (Guggenheim-Anderson-de Boer)"

    def K(self, theta: NumericType) -> NumericType:
        """Compute temperature-dependent K parameter."""
        K0 = self.params["K0"]
        qK = self.params["qK"]
        return K0 * np.exp(qK / R_SPECIFIC_GAS / theta)

    def C(self, theta: NumericType) -> NumericType:
        """Compute temperature-dependent C parameter."""
        C0 = self.params["C0"]
        qC = self.params["qC"]
        return C0 * np.exp(qC / R_SPECIFIC_GAS / theta)

    def omega_fsp_gab(self, theta: NumericType) -> NumericType:
        """Compute GAB-specific fiber saturation point."""
        mGAB = self.params["mGAB"]
        C_val = self.C(theta)
        K_val = self.K(theta)
        return mGAB * C_val * K_val / (1 - K_val) / (1 + (C_val - 1) * K_val)

    def compute_water_activity(
        self, omega: NumericType, theta: NumericType
    ) -> NumericType:
        """Compute water activity for standard GAB model."""
        mGAB = self.params["mGAB"]
        C_val = self.C(theta)
        K_val = self.K(theta)

        term = C_val ** 2 * (mGAB - omega) ** 2 + 4 * C_val * mGAB * omega
        sqrt_term = np.sqrt(term)

        numerator = sqrt_term - C_val * (mGAB - omega) - 2 * omega
        denominator = 2 * (C_val - 1) * K_val * omega

        # Avoid division by zero
        denominator_safe = np.where(
            np.abs(denominator) < 1e-10, 1e-10, denominator
        )

        return numerator / denominator_safe

    def compute_moisture_ratio(
        self, aw: NumericType, theta: NumericType
    ) -> NumericType:
        """Compute moisture ratio from water activity."""
        mGAB = self.params["mGAB"]
        C_val = self.C(theta)
        K_val = self.K(theta)

        return (
            mGAB
            * C_val
            * K_val
            * aw
            / (1 - K_val * aw)
            / (1 - K_val * aw + C_val * K_val * aw)
        )

    def compute_net_heat_of_sorption(
        self, omega: NumericType, theta: NumericType
    ) -> NumericType:
        """Compute net heat of sorption for standard GAB model."""
        C0 = self.params["C0"]
        qC = self.params["qC"]
        K0 = self.params["K0"]
        qK = self.params["qK"]
        mGAB = self.params["mGAB"]

        C_val = self.C(theta)
        K_val = self.K(theta)
        aw = self.compute_water_activity(omega, theta)

        term = (1 - K_val * aw) ** 2
        denominator = 1 + (C_val - 1) * K_val ** 2 * aw ** 2

        return term / denominator * qC + qK

    def _compute_h_function(self, aw: NumericType) -> NumericType:
        """Not used for GAB model (uses custom heat calculation)."""
        return np.ones_like(aw)  # Placeholder

    # Override for GAB-specific omega_fsp
    def omega_fsp(self, theta: NumericType) -> NumericType:
        """Use GAB-specific omega_fsp calculation."""
        return self.omega_fsp_gab(theta)
