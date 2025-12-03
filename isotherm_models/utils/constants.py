"""
Physical constants and default parameters for isotherm models.

This module contains all physical constants, reference values, and model
parameters used in the isotherm calculations.

References:
    Tryding et al. (2022). A full-range moisture sorption model for
    cellulose-based materials yielding consistent net isosteric heat of sorption.
    Drying Technology. DOI: 10.1080/07373937.2022.2084104

    Leuk et al. (2016). Drying Technology 34(5):563-573
"""

# Physical constants
R_SPECIFIC_GAS = 461.5e-3  # Specific gas constant for water vapor [kJ/kg/K]

# Dimensionless constants
M_DIMENSIONLESS = 3  # Dimensionless constant used in eta calculation

# Reference temperatures
THETA_REF = 296.15  # Reference temperature [K]
THETA1 = 310        # Temperature parameter for omegaFSP [K]

# Default temperature values
DEFAULT_TEMPERATURES_C = [25, 80]  # Default isotherm temperatures [°C]


def celsius_to_kelvin(temp_c):
    """
    Convert temperature from Celsius to Kelvin.

    Parameters
    ----------
    temp_c : float or array-like
        Temperature in Celsius

    Returns
    -------
    float or array-like
        Temperature in Kelvin
    """
    return temp_c + 273.15


def kelvin_to_celsius(temp_k):
    """
    Convert temperature from Kelvin to Celsius.

    Parameters
    ----------
    temp_k : float or array-like
        Temperature in Kelvin

    Returns
    -------
    float or array-like
        Temperature in Celsius
    """
    return temp_k - 273.15


# Model parameters calibrated to bleached fiber data
# Source: Leuk et al. Drying Technology 34(5):563-573, 2016
MODEL_PARAMETERS = {
    "Henderson I": {
        "params": [0.5786, 0.0629, 377.8934, 4.7100],
        "param_names": ["c", "kappa_inf", "theta0", "n"],
        "omega_ref": 2.49,
        "description": "Henderson I model for moisture sorption"
    },
    "Henderson II": {
        "params": [0.4979, 0.0519, 375.3296, 4.7748, 1.4827],
        "param_names": ["c", "kappa_inf", "theta0", "n", "b"],
        "omega_ref": 2.49,
        "description": "Henderson II model with additional parameter b"
    },
    "Oswin I": {
        "params": [0.4306, 0.0541, 364.0613, 6.0171],
        "param_names": ["c", "kappa_inf", "theta0", "n"],
        "omega_ref": 2.49,
        "description": "Oswin I model for moisture sorption"
    },
    "Oswin II": {
        "params": [0.5837, 0.0694, 373.8213, 5.0128, 0.5539],
        "param_names": ["c", "kappa_inf", "theta0", "n", "b"],
        "omega_ref": 2.49,
        "description": "Oswin II model with additional parameter b"
    },
    "LW I": {
        "params": [0.3644, 0.0491, 359.4227, 6.7272],
        "param_names": ["c", "kappa_inf", "theta0", "n"],
        "omega_ref": 2.49,
        "description": "Leitner-Wolkoff I model"
    },
    "LW II": {
        "params": [0.5125, 0.0617, 370.8969, 5.2068, 0.6089],
        "param_names": ["c", "kappa_inf", "theta0", "n", "b"],
        "omega_ref": 2.49,
        "description": "Leitner-Wolkoff II model with additional parameter b"
    },
    "GAB I": {
        "params": [34.9718, -86.0959, 1.0098, -4.1838],
        "param_names": ["C0", "qC", "K0", "qK"],
        "omega_ref": 2.49,
        "description": "GAB I model (modified Guggenheim-Anderson-de Boer)"
    },
    "GAB": {
        "params": [0.4273, 468.9266, 0.1011, 297.7358, 0.0578],
        "param_names": ["C0", "qC", "K0", "qK", "mGAB"],
        "omega_ref": None,
        "description": "GAB model (Guggenheim-Anderson-de Boer) - standard form"
    }
}


# Plotting configuration
PLOT_COLORS = {
    "Henderson I": 'r',
    "Henderson II": 'r',
    "Oswin I": 'g',
    "Oswin II": 'g',
    "LW I": 'b',
    "LW II": 'b',
    "GAB I": 'k',
    "GAB": 'k'
}

PLOT_STYLES = {
    "T1": '-',   # Solid line for first temperature
    "T2": ':'    # Dotted line for second temperature
}
