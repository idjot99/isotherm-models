"""
Henderson isotherm models for moisture sorption.

This module implements Henderson I and Henderson II models for predicting
moisture sorption isotherms in cellulose-based materials.
"""

import numpy as np
from typing import Union

from .base import IsothermModel, NumericType
from ..utils import cardanos_formula_cubic_root, safe_log


class HendersonI(IsothermModel):
    """
    Henderson I isotherm model.

    This model uses a modified Henderson equation to describe moisture
    sorption behavior in cellulose materials.

    Parameters
    ----------
    c : float
        Model parameter c
    kappa_inf : float
        Asymptotic value of kappa at high temperature
    theta0 : float
        Reference temperature parameter [K]
    n : float
        Power law exponent
    omega_ref : float
        Reference moisture ratio value

    References
    ----------
    Tryding et al. (2022). DOI: 10.1080/07373937.2022.2084104
    """

    @property
    def model_name(self) -> str:
        """Return the model name."""
        return "Henderson I"

    @property
    def description(self) -> str:
        """Return model description."""
        return "Henderson I model for moisture sorption in cellulose materials"

    def _compute_xi(self, omega: NumericType, theta: NumericType) -> NumericType:
        """
        Compute the intermediate variable xi.

        Parameters
        ----------
        omega : float or np.ndarray
            Moisture ratio [-]
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            xi value
        """
        c = self.params["c"]
        kappa_val = self.kappa(theta)
        omega_fsp_val = self.omega_fsp(theta)

        return (omega / kappa_val / (1 - (omega / omega_fsp_val) ** self.m)) ** (1 / c)

    def compute_water_activity(
        self, omega: NumericType, theta: NumericType
    ) -> NumericType:
        """
        Compute water activity using Henderson I equation.

        Parameters
        ----------
        omega : float or np.ndarray
            Moisture ratio [-]
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            Water activity [-]

        Notes
        -----
        aw = 1 - exp(-xi)
        """
        xi = self._compute_xi(omega, theta)
        return 1 - np.exp(-xi)

    def _compute_v(self, aw: NumericType, theta: NumericType) -> NumericType:
        """
        Compute the v parameter for inverse calculation.

        Parameters
        ----------
        aw : float or np.ndarray
            Water activity [-]
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            v parameter value
        """
        c = self.params["c"]
        kappa_val = self.kappa(theta)
        omega_fsp_val = self.omega_fsp(theta)

        return omega_fsp_val / (3 * kappa_val * (-safe_log(1 - aw)) ** c)

    def compute_moisture_ratio(
        self, aw: NumericType, theta: NumericType
    ) -> NumericType:
        """
        Compute moisture ratio from water activity (inverse function).

        Parameters
        ----------
        aw : float or np.ndarray
            Water activity [-]
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            Moisture ratio [-]

        Notes
        -----
        Uses Cardano's formula to solve the cubic equation for m=3.
        """
        v_val = self._compute_v(aw, theta)
        omega_fsp_val = self.omega_fsp(theta)

        return omega_fsp_val * cardanos_formula_cubic_root(v_val)

    def _compute_h_function(self, aw: NumericType) -> NumericType:
        """
        Compute the h function for Henderson I model.

        Parameters
        ----------
        aw : float or np.ndarray
            Water activity [-]

        Returns
        -------
        float or np.ndarray
            h function value

        Notes
        -----
        h(aw) = -(1-aw)/aw * log(1-aw)
        """
        # Avoid division by zero and log(0)
        aw_safe = np.clip(aw, 1e-10, 1 - 1e-10)
        return -(1 - aw_safe) / aw_safe * safe_log(1 - aw_safe)


class HendersonII(IsothermModel):
    """
    Henderson II isotherm model with additional parameter b.

    This is an extended version of Henderson I with an additional parameter
    for improved fitting flexibility.

    Parameters
    ----------
    c : float
        Model parameter c
    kappa_inf : float
        Asymptotic value of kappa at high temperature
    theta0 : float
        Reference temperature parameter [K]
    n : float
        Power law exponent
    b : float
        Additional model parameter
    omega_ref : float
        Reference moisture ratio value

    References
    ----------
    Tryding et al. (2022). DOI: 10.1080/07373937.2022.2084104
    """

    @property
    def model_name(self) -> str:
        """Return the model name."""
        return "Henderson II"

    @property
    def description(self) -> str:
        """Return model description."""
        return "Henderson II model with extended parameterization"

    def _compute_xi(self, omega: NumericType, theta: NumericType) -> NumericType:
        """
        Compute the intermediate variable xi.

        Parameters
        ----------
        omega : float or np.ndarray
            Moisture ratio [-]
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            xi value
        """
        c = self.params["c"]
        kappa_val = self.kappa(theta)
        omega_fsp_val = self.omega_fsp(theta)

        return (omega / kappa_val / (1 - (omega / omega_fsp_val) ** self.m)) ** (1 / c)

    def compute_water_activity(
        self, omega: NumericType, theta: NumericType
    ) -> NumericType:
        """
        Compute water activity using Henderson II equation.

        Parameters
        ----------
        omega : float or np.ndarray
            Moisture ratio [-]
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            Water activity [-]

        Notes
        -----
        aw = (1 - exp(-xi^(1/b)))^b
        """
        b = self.params["b"]
        xi = self._compute_xi(omega, theta)
        return (1 - np.exp(-(xi ** (1 / b)))) ** b

    def _compute_v(self, aw: NumericType, theta: NumericType) -> NumericType:
        """
        Compute the v parameter for inverse calculation.

        Parameters
        ----------
        aw : float or np.ndarray
            Water activity [-]
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            v parameter value
        """
        c = self.params["c"]
        b = self.params["b"]
        kappa_val = self.kappa(theta)
        omega_fsp_val = self.omega_fsp(theta)

        aw_safe = np.clip(aw, 1e-10, 1 - 1e-10)
        return omega_fsp_val / (
            3 * kappa_val * ((-safe_log(1 - aw_safe ** (1 / b))) ** b) ** c
        )

    def compute_moisture_ratio(
        self, aw: NumericType, theta: NumericType
    ) -> NumericType:
        """
        Compute moisture ratio from water activity (inverse function).

        Parameters
        ----------
        aw : float or np.ndarray
            Water activity [-]
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            Moisture ratio [-]

        Notes
        -----
        Uses Cardano's formula to solve the cubic equation for m=3.
        """
        v_val = self._compute_v(aw, theta)
        omega_fsp_val = self.omega_fsp(theta)

        return omega_fsp_val * cardanos_formula_cubic_root(v_val)

    def _compute_h_function(self, aw: NumericType) -> NumericType:
        """
        Compute the h function for Henderson II model.

        Parameters
        ----------
        aw : float or np.ndarray
            Water activity [-]

        Returns
        -------
        float or np.ndarray
            h function value

        Notes
        -----
        h(aw, b) = -(1-aw^(1/b))/aw^(1/b) * log(1-aw^(1/b))
        """
        b = self.params["b"]
        aw_safe = np.clip(aw, 1e-10, 1 - 1e-10)
        aw_b = aw_safe ** (1 / b)
        aw_b_safe = np.clip(aw_b, 1e-10, 1 - 1e-10)

        return -(1 - aw_b_safe) / aw_b_safe * safe_log(1 - aw_b_safe)
