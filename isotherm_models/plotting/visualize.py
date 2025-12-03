"""
Visualization tools for isotherm models.

This module provides functions for plotting moisture sorption isotherms
and heat of sorption curves.
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import List, Optional, Dict
from ..models import IsothermModel
from ..utils import PLOT_COLORS, PLOT_STYLES


def plot_water_activity(
    models: List[IsothermModel],
    temperatures_k: List[float],
    omega_range: Optional[tuple] = None,
    show_inverse: bool = True,
    figsize: tuple = (10, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plot water activity vs moisture ratio isotherms.

    Parameters
    ----------
    models : List[IsothermModel]
        List of isotherm models to plot
    temperatures_k : List[float]
        List of temperatures in Kelvin
    omega_range : tuple, optional
        Range for omega axis (min, max). If None, uses (0, 0.3)
    show_inverse : bool, optional
        If True, also plot the inverse relationship (omega vs aw)
    figsize : tuple, optional
        Figure size (width, height) in inches
    save_path : str, optional
        If provided, save the figure to this path

    Returns
    -------
    plt.Figure
        Matplotlib figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlabel('Moisture ratio, ω [-]', fontsize=14)
    ax.set_ylabel('Water activity, a_ω [-]', fontsize=14)

    if omega_range is None:
        omega_range = (0, 0.3)

    for model in models:
        color = PLOT_COLORS.get(model.model_name, 'k')

        for i, temp in enumerate(temperatures_k):
            line_style = PLOT_STYLES["T1"] if i == 0 else PLOT_STYLES["T2"]

            # Compute isotherm
            omega_fsp = model.omega_fsp(temp)
            omega_vals = np.linspace(0, min(omega_fsp * 0.99, omega_range[1]), 300)
            aw_vals = model.compute_water_activity(omega_vals, temp)

            # Plot main isotherm
            label = f'{model.model_name}, T={temp-273.15:.0f}°C' if i == 0 else None
            ax.plot(omega_vals, aw_vals, line_style, color=color, linewidth=2, label=label)

            # Plot inverse if requested
            if show_inverse:
                aw_inv = np.linspace(0.01, 0.99, 300)
                omega_inv = model.compute_moisture_ratio(aw_inv, temp)
                ax.plot(omega_inv, aw_inv, 'w--', linewidth=1, alpha=0.7)

    ax.set_xlim(omega_range)
    ax.set_ylim(0, 1)
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig


def plot_heat_of_sorption(
    models: List[IsothermModel],
    temperatures_k: List[float],
    omega_range: Optional[tuple] = None,
    hiso_range: Optional[tuple] = None,
    figsize: tuple = (10, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plot net isosteric heat of sorption vs moisture ratio.

    Parameters
    ----------
    models : List[IsothermModel]
        List of isotherm models to plot
    temperatures_k : List[float]
        List of temperatures in Kelvin
    omega_range : tuple, optional
        Range for omega axis (min, max). If None, uses (0, 0.45)
    hiso_range : tuple, optional
        Range for heat axis (min, max). If None, uses (0, 1700)
    figsize : tuple, optional
        Figure size (width, height) in inches
    save_path : str, optional
        If provided, save the figure to this path

    Returns
    -------
    plt.Figure
        Matplotlib figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlabel('Moisture ratio, ω [-]', fontsize=14)
    ax.set_ylabel('Net heat of sorption, ΔH_iso [kJ/kg]', fontsize=14)

    if omega_range is None:
        omega_range = (0, 0.45)
    if hiso_range is None:
        hiso_range = (0, 1700)

    for model in models:
        color = PLOT_COLORS.get(model.model_name, 'k')

        for i, temp in enumerate(temperatures_k):
            line_style = PLOT_STYLES["T1"] if i == 0 else PLOT_STYLES["T2"]

            # Compute heat of sorption
            omega_fsp = model.omega_fsp(temp)
            omega_vals = np.linspace(0, min(omega_fsp * 0.99, omega_range[1]), 300)
            hiso_vals = model.compute_net_heat_of_sorption(omega_vals, temp)

            # Plot
            label = f'{model.model_name}, T={temp-273.15:.0f}°C' if i == 0 else None
            ax.plot(omega_vals, hiso_vals, line_style, color=color, linewidth=2, label=label)

    ax.set_xlim(omega_range)
    ax.set_ylim(hiso_range)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig


def plot_combined(
    models: List[IsothermModel],
    temperatures_k: List[float],
    save_dir: Optional[str] = None,
) -> Dict[str, plt.Figure]:
    """
    Create both water activity and heat of sorption plots.

    Parameters
    ----------
    models : List[IsothermModel]
        List of isotherm models to plot
    temperatures_k : List[float]
        List of temperatures in Kelvin
    save_dir : str, optional
        If provided, save figures to this directory

    Returns
    -------
    Dict[str, plt.Figure]
        Dictionary with keys 'water_activity' and 'heat_of_sorption'
    """
    import os

    fig_aw = plot_water_activity(
        models,
        temperatures_k,
        save_path=os.path.join(save_dir, 'water_activity.png') if save_dir else None
    )

    fig_hiso = plot_heat_of_sorption(
        models,
        temperatures_k,
        save_path=os.path.join(save_dir, 'heat_of_sorption.png') if save_dir else None
    )

    return {'water_activity': fig_aw, 'heat_of_sorption': fig_hiso}
