"""
Unit tests for Henderson models.
"""

import pytest
import numpy as np
from isotherm_models import HendersonI, HendersonII, create_model


class TestHendersonI:
    """Test suite for Henderson I model."""

    @pytest.fixture
    def model(self):
        """Create a Henderson I model instance with default parameters."""
        return create_model("Henderson I")

    def test_model_creation(self, model):
        """Test that model can be created successfully."""
        assert model is not None
        assert model.model_name == "Henderson I"
        assert "Henderson" in model.description

    def test_water_activity_range(self, model):
        """Test that water activity is in valid range [0, 1]."""
        theta = 298.15  # 25°C
        omega_vals = np.linspace(0, 0.2, 50)
        aw_vals = model.compute_water_activity(omega_vals, theta)

        assert np.all(aw_vals >= 0)
        assert np.all(aw_vals <= 1)

    def test_water_activity_monotonic(self, model):
        """Test that water activity increases monotonically with omega."""
        theta = 298.15
        omega_vals = np.linspace(0.01, 0.2, 50)
        aw_vals = model.compute_water_activity(omega_vals, theta)

        # Check that aw increases
        assert np.all(np.diff(aw_vals) > 0)

    def test_inverse_consistency(self, model):
        """Test that inverse function is consistent."""
        theta = 298.15
        omega_original = 0.1

        # Forward: omega -> aw
        aw = model.compute_water_activity(omega_original, theta)

        # Inverse: aw -> omega
        omega_computed = model.compute_moisture_ratio(aw, theta)

        # Should be approximately equal
        assert np.isclose(omega_original, omega_computed, rtol=1e-3)

    def test_temperature_dependency(self, model):
        """Test that water activity changes with temperature."""
        omega = 0.1
        theta1 = 298.15  # 25°C
        theta2 = 353.15  # 80°C

        aw1 = model.compute_water_activity(omega, theta1)
        aw2 = model.compute_water_activity(omega, theta2)

        # Water activity should be different at different temperatures
        assert not np.isclose(aw1, aw2, rtol=0.01)

    def test_heat_of_sorption_positive(self, model):
        """Test that heat of sorption is positive."""
        theta = 298.15
        omega_vals = np.linspace(0.01, 0.2, 50)
        hiso_vals = model.compute_net_heat_of_sorption(omega_vals, theta)

        assert np.all(hiso_vals > 0)

    def test_heat_of_sorption_decreases(self, model):
        """Test that heat of sorption generally decreases with omega."""
        theta = 298.15
        omega_vals = np.linspace(0.01, 0.15, 50)
        hiso_vals = model.compute_net_heat_of_sorption(omega_vals, theta)

        # Should generally decrease (allow some small fluctuations)
        assert hiso_vals[0] > hiso_vals[-1]

    def test_compute_isotherm(self, model):
        """Test complete isotherm computation."""
        theta = 298.15
        isotherm = model.compute_sorption_isotherm(theta, n_points=100)

        assert 'omega' in isotherm
        assert 'aw' in isotherm
        assert 'Hiso' in isotherm
        assert len(isotherm['omega']) == 100
        assert len(isotherm['aw']) == 100
        assert len(isotherm['Hiso']) == 100


class TestHendersonII:
    """Test suite for Henderson II model."""

    @pytest.fixture
    def model(self):
        """Create a Henderson II model instance."""
        return create_model("Henderson II")

    def test_model_creation(self, model):
        """Test that model can be created successfully."""
        assert model is not None
        assert model.model_name == "Henderson II"

    def test_water_activity_range(self, model):
        """Test that water activity is in valid range."""
        theta = 298.15
        omega_vals = np.linspace(0, 0.2, 50)
        aw_vals = model.compute_water_activity(omega_vals, theta)

        assert np.all(aw_vals >= 0)
        assert np.all(aw_vals <= 1)

    def test_inverse_consistency(self, model):
        """Test that inverse function is consistent."""
        theta = 298.15
        omega_original = 0.1

        aw = model.compute_water_activity(omega_original, theta)
        omega_computed = model.compute_moisture_ratio(aw, theta)

        assert np.isclose(omega_original, omega_computed, rtol=1e-3)

    def test_b_parameter_effect(self):
        """Test that parameter b affects the model."""
        theta = 298.15
        omega = 0.1

        model1 = create_model("Henderson I")
        model2 = create_model("Henderson II")

        aw1 = model1.compute_water_activity(omega, theta)
        aw2 = model2.compute_water_activity(omega, theta)

        # Henderson II (with b parameter) should give different result
        assert not np.isclose(aw1, aw2, rtol=0.01)


def test_parameter_validation():
    """Test that invalid parameters raise errors."""
    # Should work with valid parameters
    model = create_model("Henderson I")
    assert model is not None

    # Should fail with negative omega_ref
    with pytest.raises(ValueError):
        HendersonI(
            params={"c": 0.5, "kappa_inf": 0.06, "theta0": 377, "n": 4.7},
            omega_ref=-1.0
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
