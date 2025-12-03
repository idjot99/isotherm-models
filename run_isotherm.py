#!/usr/bin/env python3
"""
Command-line interface for isotherm models.

This script allows users to compute and visualize moisture sorption isotherms
for various models from the command line.
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from isotherm_models import (
    create_model,
    MODEL_REGISTRY,
    plot_water_activity,
    plot_heat_of_sorption,
    celsius_to_kelvin,
)


def main():
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        description="Compute and visualize moisture sorption isotherms",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot Henderson I model at 25°C and 80°C
  python run_isotherm.py --model "Henderson I" --temperatures 25 80

  # Plot multiple models
  python run_isotherm.py --model "Henderson I" "Oswin I" "LW I" --temperatures 25 80

  # Save plots to directory
  python run_isotherm.py --model "GAB" --temperatures 25 80 --output plots/

  # List available models
  python run_isotherm.py --list-models
        """
    )

    parser.add_argument(
        "--model",
        "-m",
        nargs="+",
        help="Model name(s) to use. Use --list-models to see available models.",
    )

    parser.add_argument(
        "--temperatures",
        "-t",
        nargs="+",
        type=float,
        default=[25, 80],
        help="Temperatures in Celsius (default: 25 80)",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="Output directory for saving plots. If not provided, plots are shown interactively.",
    )

    parser.add_argument(
        "--list-models",
        action="store_true",
        help="List all available models and exit",
    )

    parser.add_argument(
        "--no-inverse",
        action="store_true",
        help="Do not plot inverse relationship on water activity plot",
    )

    parser.add_argument(
        "--format",
        choices=["png", "pdf", "svg"],
        default="png",
        help="Output format for saved figures (default: png)",
    )

    args = parser.parse_args()

    # List models if requested
    if args.list_models:
        print("\nAvailable isotherm models:")
        print("-" * 60)
        for name, model_class in MODEL_REGISTRY.items():
            # Create temporary instance to get description
            try:
                model = create_model(name)
                print(f"  {name:20s} - {model.description}")
            except Exception as e:
                print(f"  {name:20s} - (Error loading: {e})")
        print("-" * 60)
        return 0

    # Validate that model is provided
    if not args.model:
        parser.error("--model is required unless --list-models is used")

    # Validate model names
    for model_name in args.model:
        if model_name not in MODEL_REGISTRY:
            print(f"Error: Unknown model '{model_name}'", file=sys.stderr)
            print("Use --list-models to see available models", file=sys.stderr)
            return 1

    # Convert temperatures to Kelvin
    temperatures_k = [celsius_to_kelvin(t) for t in args.temperatures]

    # Create models
    print(f"\nCreating {len(args.model)} model(s)...")
    models = []
    for model_name in args.model:
        try:
            model = create_model(model_name)
            models.append(model)
            print(f"  ✓ {model_name}")
        except Exception as e:
            print(f"  ✗ {model_name}: {e}", file=sys.stderr)
            return 1

    # Create output directory if needed
    if args.output:
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"\nSaving plots to: {output_dir}")

    # Plot water activity
    print("\nGenerating water activity plot...")
    fig_aw = plot_water_activity(
        models,
        temperatures_k,
        show_inverse=not args.no_inverse,
        save_path=str(output_dir / f"water_activity.{args.format}") if args.output else None,
    )

    # Plot heat of sorption
    print("Generating heat of sorption plot...")
    fig_hiso = plot_heat_of_sorption(
        models,
        temperatures_k,
        save_path=str(output_dir / f"heat_of_sorption.{args.format}") if args.output else None,
    )

    # Show plots if not saving
    if not args.output:
        print("\nDisplaying plots... (close windows to exit)")
        plt.show()
    else:
        print(f"\n✓ Plots saved successfully!")

    return 0


if __name__ == "__main__":
    sys.exit(main())
