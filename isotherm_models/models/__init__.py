"""Isotherm models for moisture sorption in cellulose-based materials."""

from .base import IsothermModel
from .henderson import HendersonI, HendersonII
from .oswin import OswinI, OswinII
from .lw import LWI, LWII
from .gab import GABI, GAB

# Model registry for easy access by name
MODEL_REGISTRY = {
    "Henderson I": HendersonI,
    "Henderson II": HendersonII,
    "Oswin I": OswinI,
    "Oswin II": OswinII,
    "LW I": LWI,
    "LW II": LWII,
    "GAB I": GABI,
    "GAB": GAB,
}

__all__ = [
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
]
