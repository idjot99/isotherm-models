"""
Isotherm Models Package

A comprehensive package for modeling moisture sorption in cellulose-based materials.

This package implements various isotherm models including Henderson, Oswin,
Leitner-Wolkoff (LW), and GAB models, providing tools for computing water activity,
moisture ratio, and net isosteric heat of sorption.

References
----------
Tryding et al. (2022). A full-range moisture sorption model for
cellulose-based materials yielding consistent net isosteric heat of sorption.
Drying Technology. DOI: 10.1080/07373937.2022.2084104
"""

__version__ = "1.0.0"

from . import models, utils, plotting
from .models import (
    IsothermModel,
    HendersonI,
    HendersonII,
    OswinI,
    OswinII,
    LWI,
    LWII,
    GABI,
    GAB,
    MODEL_REGISTRY,
)
from .utils import MODEL_PARAMETERS, celsius_to_kelvin
from .plotting import plot_water_activity, plot_heat_of_sorption, plot_combined


def create_model(model_name: str, **kwargs):
    """
    Factory function to create isotherm models by name.

    Parameters
    ----------
    model_name : str
        Name of the model (e.g., "Henderson I", "Oswin II", "GAB")
    **kwargs : dict
        Additional keyword arguments to pass to the model constructor.
        If not provided, default parameters from MODEL_PARAMETERS will be used.

    Returns
    -------
    IsothermModel
        Instance of the requested model

    Examples
    --------
    >>> model = create_model("Henderson I")
    >>> aw = model.compute_water_activity(0.1, 298.15)

    >>> model = create_model("Oswin II", params={"c": 0.6, "kappa_inf": 0.07, ...})
    """
    if model_name not in MODEL_REGISTRY:
        available = ", ".join(MODEL_REGISTRY.keys())
        raise ValueError(
            f"Unknown model '{model_name}'. Available models: {available}"
        )

    # Get default parameters if not provided
    if "params" not in kwargs and model_name in MODEL_PARAMETERS:
        param_data = MODEL_PARAMETERS[model_name]
        params_dict = dict(zip(param_data["param_names"], param_data["params"]))
        kwargs["params"] = params_dict

        if "omega_ref" not in kwargs and param_data["omega_ref"] is not None:
            kwargs["omega_ref"] = param_data["omega_ref"]

    model_class = MODEL_REGISTRY[model_name]
    return model_class(**kwargs)


__all__ = [
    "__version__",
    "models",
    "utils",
    "plotting",
    "IsothermModel",
    "HendersonI",
    "HendersonII",
    "OswinI",
    "OswinII",
    "LWI",
    "LWII",
    "GABI",
    "GAB",
    "MODEL_REGISTRY",
    "create_model",
    "plot_water_activity",
    "plot_heat_of_sorption",
    "plot_combined",
]
