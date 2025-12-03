"""
Mathematical utilities for isotherm calculations.

This module provides mathematical functions used across different isotherm models,
including the Lambert W function and common transformations.
"""

import numpy as np
from scipy.special import lambertw
from typing import Union

# Type alias for numeric types
NumericType = Union[float, np.ndarray]


def lambert_w(x: NumericType) -> NumericType:
    """
    Compute the Lambert W function (principal branch).

    The Lambert W function is the inverse of f(w) = w * exp(w).
    Used in Leitner-Wolkoff (LW) models for water activity calculations.

    Parameters
    ----------
    x : float or np.ndarray
        Input value(s), must be >= -1/e for real solutions

    Returns
    -------
    float or np.ndarray
        The real part of W(x)

    References
    ----------
    Corless et al. (1996). On the Lambert W function.
    Advances in Computational Mathematics 5:329-359

    Notes
    -----
    Uses scipy's lambertw function which implements efficient algorithms
    for computing the Lambert W function.

    Examples
    --------
    >>> lambert_w(0)
    0.0
    >>> lambert_w(np.e)
    1.0
    """
    return lambertw(x).real


def cardanos_formula_cubic_root(v: NumericType) -> NumericType:
    """
    Solve the depressed cubic equation x³ + 3v*x - v = 0 using Cardano's formula.

    This is used to invert the moisture content relationship for models with m=3.
    The equation arises from omega = omegaFSP * (x³ - v*x) where x is the solution.

    Parameters
    ----------
    v : float or np.ndarray
        The parameter v from the depressed cubic equation

    Returns
    -------
    float or np.ndarray
        The real root of the cubic equation

    Notes
    -----
    The formula used is:
    x = (0.5 + sqrt(0.25 + v³))^(1/3) - v / (0.5 + sqrt(0.25 + v³))^(1/3)

    This provides the real root for all real values of v.

    References
    ----------
    Based on Cardano's formula for solving cubic equations.

    Examples
    --------
    >>> cardanos_formula_cubic_root(0.1)
    0.316...
    """
    sqrt_term = np.sqrt(0.25 + v ** 3)
    cube_root_term = (0.5 + sqrt_term) ** (1 / 3)
    return cube_root_term - v / cube_root_term


def validate_water_activity(aw: NumericType, name: str = "water activity") -> None:
    """
    Validate that water activity values are in the valid range [0, 1].

    Parameters
    ----------
    aw : float or np.ndarray
        Water activity value(s) to validate
    name : str, optional
        Name of the variable for error message (default: "water activity")

    Raises
    ------
    ValueError
        If any water activity value is outside [0, 1]

    Examples
    --------
    >>> validate_water_activity(0.5)  # No error
    >>> validate_water_activity(1.5)  # Raises ValueError
    """
    if np.any(aw < 0) or np.any(aw > 1):
        raise ValueError(f"{name} must be between 0 and 1, got {aw}")


def validate_positive(value: NumericType, name: str) -> None:
    """
    Validate that a value is positive.

    Parameters
    ----------
    value : float or np.ndarray
        Value(s) to validate
    name : str
        Name of the variable for error message

    Raises
    ------
    ValueError
        If any value is <= 0

    Examples
    --------
    >>> validate_positive(1.0, "temperature")  # No error
    >>> validate_positive(-1.0, "temperature")  # Raises ValueError
    """
    if np.any(value <= 0):
        raise ValueError(f"{name} must be positive, got {value}")


def validate_moisture_ratio(omega: NumericType, omega_fsp: NumericType) -> None:
    """
    Validate that moisture ratio is in the valid range [0, omega_FSP].

    Parameters
    ----------
    omega : float or np.ndarray
        Moisture ratio value(s) to validate
    omega_fsp : float or np.ndarray
        Fiber saturation point value(s)

    Raises
    ------
    ValueError
        If any omega value is outside [0, omega_fsp]

    Examples
    --------
    >>> validate_moisture_ratio(0.1, 0.3)  # No error
    >>> validate_moisture_ratio(0.5, 0.3)  # Raises ValueError
    """
    if np.any(omega < 0):
        raise ValueError(f"Moisture ratio must be non-negative, got {omega}")
    if np.any(omega > omega_fsp):
        raise ValueError(
            f"Moisture ratio must not exceed omega_FSP={omega_fsp}, got {omega}"
        )


def safe_log(x: NumericType, epsilon: float = 1e-10) -> NumericType:
    """
    Compute natural logarithm with protection against log(0).

    Parameters
    ----------
    x : float or np.ndarray
        Input value(s)
    epsilon : float, optional
        Small value to replace zero (default: 1e-10)

    Returns
    -------
    float or np.ndarray
        Natural logarithm of max(x, epsilon)

    Notes
    -----
    This function replaces values <= 0 with epsilon before computing
    the logarithm to avoid numerical errors.

    Examples
    --------
    >>> safe_log(np.e)
    1.0
    >>> safe_log(0)  # Returns log(epsilon) instead of -inf
    -23.025...
    """
    x_safe = np.maximum(x, epsilon)
    return np.log(x_safe)


def safe_divide(numerator: NumericType, denominator: NumericType,
                epsilon: float = 1e-10) -> NumericType:
    """
    Perform division with protection against division by zero.

    Parameters
    ----------
    numerator : float or np.ndarray
        Numerator value(s)
    denominator : float or np.ndarray
        Denominator value(s)
    epsilon : float, optional
        Small value to replace zero denominator (default: 1e-10)

    Returns
    -------
    float or np.ndarray
        Result of numerator / max(denominator, epsilon)

    Examples
    --------
    >>> safe_divide(1.0, 2.0)
    0.5
    >>> safe_divide(1.0, 0.0)  # Returns 1.0 / epsilon instead of inf
    10000000000.0
    """
    denominator_safe = np.where(
        np.abs(denominator) < epsilon,
        epsilon,
        denominator
    )
    return numerator / denominator_safe
