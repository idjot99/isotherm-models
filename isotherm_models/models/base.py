"""
Base class for isotherm models.

This module provides the abstract base class that all isotherm models inherit from.
It defines the common interface and shared functionality for moisture sorption calculations.
"""

from abc import ABC, abstractmethod
from typing import Dict, Optional, Union
import numpy as np

from ..utils import (
    R_SPECIFIC_GAS,
    M_DIMENSIONLESS,
    THETA_REF,
    THETA1,
    validate_positive,
    validate_moisture_ratio,
)

NumericType = Union[float, np.ndarray]


class IsothermModel(ABC):
    """
    Abstract base class for all isotherm models.

    This class defines the interface that all specific isotherm models must implement.
    It provides common functionality for temperature-dependent calculations and
    parameter management.

    Parameters
    ----------
    params : Dict[str, float]
        Dictionary of model parameters. Required keys depend on the specific model.
    omega_ref : float
        Reference moisture ratio value (omega at reference conditions)
    theta_ref : float, optional
        Reference temperature in Kelvin (default: 296.15 K)
    theta1 : float, optional
        Temperature parameter for omega_fsp calculation (default: 310 K)
    m : int, optional
        Dimensionless constant (default: 3)

    Attributes
    ----------
    model_name : str
        Name of the isotherm model
    description : str
        Description of the model

    References
    ----------
    Tryding et al. (2022). A full-range moisture sorption model for
    cellulose-based materials yielding consistent net isosteric heat of sorption.
    Drying Technology. DOI: 10.1080/07373937.2022.2084104
    """

    def __init__(
        self,
        params: Dict[str, float],
        omega_ref: float,
        theta_ref: float = THETA_REF,
        theta1: float = THETA1,
        m: int = M_DIMENSIONLESS,
    ):
        """Initialize the isotherm model with parameters."""
        self.params = params
        self.omega_ref = omega_ref
        self.theta_ref = theta_ref
        self.theta1 = theta1
        self.m = m

        # Validate parameters
        validate_positive(omega_ref, "omega_ref")
        validate_positive(theta_ref, "theta_ref")
        validate_positive(theta1, "theta1")

    @property
    @abstractmethod
    def model_name(self) -> str:
        """Return the name of the model."""
        pass

    @property
    @abstractmethod
    def description(self) -> str:
        """Return a description of the model."""
        pass

    @abstractmethod
    def compute_water_activity(
        self, omega: NumericType, theta: NumericType
    ) -> NumericType:
        """
        Compute water activity as a function of moisture ratio and temperature.

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
        """
        pass

    @abstractmethod
    def compute_moisture_ratio(
        self, aw: NumericType, theta: NumericType
    ) -> NumericType:
        """
        Compute moisture ratio as a function of water activity and temperature.

        This is the inverse of compute_water_activity.

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
        """
        pass

    def kappa(self, theta: NumericType) -> NumericType:
        """
        Compute the temperature-dependent parameter kappa.

        Parameters
        ----------
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            kappa parameter value

        Notes
        -----
        This function implements Equation (X) from the reference paper.
        kappa = kappa_inf * exp(1/(n+1) * (theta0/theta)^(n+1))
        """
        kappa_inf = self.params["kappa_inf"]
        theta0 = self.params["theta0"]
        n = self.params["n"]

        return kappa_inf * np.exp((1 / (n + 1)) * (theta0 / theta) ** (n + 1))

    def omega_fsp(self, theta: NumericType) -> NumericType:
        """
        Compute the fiber saturation point as a function of temperature.

        Parameters
        ----------
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            Fiber saturation point omega_FSP [-]

        Notes
        -----
        Based on Alexandersson's temperature dependency model.
        omega_FSP = omega_ref * (1 + (theta_ref - theta) / theta1)

        References
        ----------
        Alexandersson (referenced in main paper)
        """
        return self.omega_ref * (1 + (self.theta_ref - theta) / self.theta1)

    def hiso_max(self, omega: NumericType, theta: NumericType) -> NumericType:
        """
        Compute the maximum isosteric heat of sorption.

        Parameters
        ----------
        omega : float or np.ndarray
            Moisture ratio [-]
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            Maximum heat of sorption [kJ/kg]

        Notes
        -----
        Hiso_max = (1/c) * R * theta0 * (theta0/theta)^n
        """
        c = self.params["c"]
        theta0 = self.params["theta0"]
        n = self.params["n"]

        return (1 / c) * R_SPECIFIC_GAS * theta0 * (theta0 / theta) ** n

    def eta(self, omega: NumericType, theta: NumericType) -> NumericType:
        """
        Compute the eta correction factor for heat of sorption.

        Parameters
        ----------
        omega : float or np.ndarray
            Moisture ratio [-]
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            Eta correction factor [-]

        Notes
        -----
        This accounts for the temperature dependence of omega_FSP in the
        heat of sorption calculation.
        """
        theta0 = self.params["theta0"]
        n = self.params["n"]

        omega_fsp_val = self.omega_fsp(theta)
        term1 = omega / omega_fsp_val

        return 1 + (
            self.m
            * (term1 ** self.m)
            / (1 - term1 ** self.m)
            * (self.omega_ref / omega_fsp_val)
            * (theta / theta0) ** (n + 2)
            * (theta0 / self.theta1)
        )

    def compute_net_heat_of_sorption(
        self, omega: NumericType, theta: NumericType
    ) -> NumericType:
        """
        Compute the net isosteric heat of sorption.

        Parameters
        ----------
        omega : float or np.ndarray
            Moisture ratio [-]
        theta : float or np.ndarray
            Temperature [K]

        Returns
        -------
        float or np.ndarray
            Net heat of sorption [kJ/kg]

        Notes
        -----
        This default implementation uses:
        H_iso = H_iso_max * h(aw) * eta

        Specific models may override this if they use different formulations.
        """
        aw = self.compute_water_activity(omega, theta)
        h_val = self._compute_h_function(aw)

        return self.hiso_max(omega, theta) * h_val * self.eta(omega, theta)

    @abstractmethod
    def _compute_h_function(self, aw: NumericType) -> NumericType:
        """
        Compute the model-specific h function.

        This is an internal method used in heat of sorption calculations.

        Parameters
        ----------
        aw : float or np.ndarray
            Water activity [-]

        Returns
        -------
        float or np.ndarray
            h function value [-]
        """
        pass

    def compute_sorption_isotherm(
        self, theta: NumericType, n_points: int = 300
    ) -> Dict[str, np.ndarray]:
        """
        Compute a complete sorption isotherm at the given temperature.

        Parameters
        ----------
        theta : float
            Temperature [K]
        n_points : int, optional
            Number of points in the isotherm (default: 300)

        Returns
        -------
        dict
            Dictionary with keys:
            - 'omega': moisture ratio values
            - 'aw': water activity values
            - 'Hiso': net heat of sorption values

        Examples
        --------
        >>> model = some_isotherm_model
        >>> isotherm = model.compute_sorption_isotherm(298.15)
        >>> plt.plot(isotherm['omega'], isotherm['aw'])
        """
        omega_fsp_val = self.omega_fsp(theta)
        omega_vals = np.linspace(0, omega_fsp_val * 0.99, n_points)

        aw_vals = self.compute_water_activity(omega_vals, theta)
        hiso_vals = self.compute_net_heat_of_sorption(omega_vals, theta)

        return {"omega": omega_vals, "aw": aw_vals, "Hiso": hiso_vals, "theta": theta}

    def __repr__(self) -> str:
        """Return string representation of the model."""
        return f"{self.__class__.__name__}(model_name='{self.model_name}')"
