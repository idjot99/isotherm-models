"""
Integration tests for the isotherm models package.

These tests verify that different components work together correctly.
"""

import pytest
import numpy as np
from isotherm_models import create_model, MODEL_REGISTRY, celsius_to_kelvin
from isotherm_models.plotting import plot_water_activity, plot_heat_of_sorption
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for testing
import matplotlib.pyplot as plt


def test_all_models_can_be_created():
    """Test that all registered models can be created with default parameters."""
    for model_name in MODEL_REGISTRY.keys():
        model = create_model(model_name)
        assert model is not None
        assert model.model_name == model_name


def test_all_models_compute_water_activity():
    """Test that all models can compute water activity."""
    theta = celsius_to_kelvin(25)
    omega = 0.1

    for model_name in MODEL_REGISTRY.keys():
        model = create_model(model_name)
        aw = model.compute_water_activity(omega, theta)

        # Check valid range
        assert 0 <= aw <= 1, f"{model_name} produced invalid aw: {aw}"


def test_all_models_compute_moisture_ratio():
    """Test that all models can compute moisture ratio (inverse)."""
    theta = celsius_to_kelvin(25)
    aw = 0.5

    for model_name in MODEL_REGISTRY.keys():
        model = create_model(model_name)
        omega = model.compute_moisture_ratio(aw, theta)

        # Check valid range
        omega_fsp = model.omega_fsp(theta)
        assert 0 <= omega <= omega_fsp, f"{model_name} produced invalid omega: {omega}"


def test_all_models_compute_heat_of_sorption():
    """Test that all models can compute heat of sorption."""
    theta = celsius_to_kelvin(25)
    omega = 0.1

    for model_name in MODEL_REGISTRY.keys():
        model = create_model(model_name)
        hiso = model.compute_net_heat_of_sorption(omega, theta)

        # Heat should be positive and reasonable
        assert hiso > 0, f"{model_name} produced negative heat: {hiso}"
        assert hiso < 5000, f"{model_name} produced unreasonably large heat: {hiso}"


def test_all_models_compute_complete_isotherm():
    """Test that all models can compute complete isotherms."""
    theta = celsius_to_kelvin(25)

    for model_name in MODEL_REGISTRY.keys():
        model = create_model(model_name)
        isotherm = model.compute_sorption_isotherm(theta, n_points=50)

        assert 'omega' in isotherm
        assert 'aw' in isotherm
        assert 'Hiso' in isotherm
        assert len(isotherm['omega']) == 50


def test_plotting_water_activity():
    """Test that water activity plotting works."""
    models = [
        create_model("Henderson I"),
        create_model("Oswin I"),
    ]
    temperatures_k = [celsius_to_kelvin(25), celsius_to_kelvin(80)]

    fig = plot_water_activity(models, temperatures_k)
    assert fig is not None
    plt.close(fig)


def test_plotting_heat_of_sorption():
    """Test that heat of sorption plotting works."""
    models = [
        create_model("Henderson I"),
        create_model("LW I"),
    ]
    temperatures_k = [celsius_to_kelvin(25), celsius_to_kelvin(80)]

    fig = plot_heat_of_sorption(models, temperatures_k)
    assert fig is not None
    plt.close(fig)


def test_inverse_consistency_all_models():
    """Test forward and inverse consistency for all models."""
    theta = celsius_to_kelvin(25)
    omega_test = 0.15

    for model_name in MODEL_REGISTRY.keys():
        model = create_model(model_name)

        # Forward
        aw = model.compute_water_activity(omega_test, theta)

        # Inverse
        omega_computed = model.compute_moisture_ratio(aw, theta)

        # Should match within tolerance
        assert np.isclose(omega_test, omega_computed, rtol=0.01), \
            f"{model_name} failed inverse consistency test"


def test_temperature_dependency_all_models():
    """Test that temperature affects results for all models."""
    omega = 0.1
    theta1 = celsius_to_kelvin(25)
    theta2 = celsius_to_kelvin(80)

    for model_name in MODEL_REGISTRY.keys():
        model = create_model(model_name)

        aw1 = model.compute_water_activity(omega, theta1)
        aw2 = model.compute_water_activity(omega, theta2)

        # Results should be different at different temperatures
        assert not np.isclose(aw1, aw2, rtol=0.05), \
            f"{model_name} shows no temperature dependency"


def test_omega_fsp_temperature_dependency():
    """Test that fiber saturation point varies with temperature."""
    theta1 = celsius_to_kelvin(25)
    theta2 = celsius_to_kelvin(80)

    for model_name in MODEL_REGISTRY.keys():
        if model_name == "GAB":  # GAB uses different omega_fsp calculation
            continue

        model = create_model(model_name)

        omega_fsp1 = model.omega_fsp(theta1)
        omega_fsp2 = model.omega_fsp(theta2)

        # Should decrease with increasing temperature
        assert omega_fsp1 > omega_fsp2, \
            f"{model_name} omega_FSP doesn't decrease with temperature"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
