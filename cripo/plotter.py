"""Plot utilities for CRIPO cross-section outputs."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt


SERIES = [
    ("EL.COH.", "elastic_coherent"),
    ("EL.INC.", "elastic_incoherent"),
    ("INEL.INC.", "inelastic_incoherent"),
    ("INEL.COH.", "inelastic_coherent"),
    ("ABSORPTION", "absorption"),
    ("TOTAL", "total"),
]

NEUTRON_WAVELENGTH_ANGSTROM_EV = 0.286014351


def read_dat(dat_path: str | Path) -> dict[str, list[float]]:
    """Read a legacy CRIPO `.DAT` file into column-oriented Python lists."""
    path = Path(dat_path)
    data = {
        "energy_ev": [],
        "elastic_coherent": [],
        "elastic_incoherent": [],
        "inelastic_incoherent": [],
        "inelastic_coherent": [],
        "absorption": [],
        "total": [],
    }
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) != 7:
                continue
            try:
                values = [float(part) for part in parts]
            except ValueError:
                continue
            data["energy_ev"].append(values[0])
            data["elastic_coherent"].append(values[1])
            data["elastic_incoherent"].append(values[2])
            data["inelastic_incoherent"].append(values[3])
            data["inelastic_coherent"].append(values[4])
            data["absorption"].append(values[5])
            data["total"].append(values[6])
    if not data["energy_ev"]:
        raise ValueError(f"No numeric cross-section rows found in {path}")
    return data


def energy_to_wavelength_angstrom(energy_ev: float) -> float:
    """Convert neutron energy in eV to wavelength in Angstrom."""
    return NEUTRON_WAVELENGTH_ANGSTROM_EV / math.sqrt(energy_ev)


def plot_cross_section_data(
    data: dict[str, list[float]],
    output_path: str | Path = "CRIPO_cross_sections.png",
    include_partials: bool = True,
    title: str = "CRIPO Cross Sections",
    x_axis: str = "energy",
) -> Path:
    """Plot cross sections from an in-memory data dictionary.

    The data structure is the same one returned by ``CripoModel`` and
    ``read_dat``. The plot is saved to ``output_path`` and that path is
    returned.
    """
    output = Path(output_path)
    if x_axis not in {"energy", "wavelength"}:
        raise ValueError("x_axis must be 'energy' or 'wavelength'")

    if x_axis == "energy":
        x_values = data["energy_ev"]
        xlabel = "Energy (eV)"
    else:
        x_values = [energy_to_wavelength_angstrom(value) for value in data["energy_ev"]]
        xlabel = "Neutron wavelength (Angstrom)"

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x_values, data["total"], label="TOTAL", linewidth=2.2, color="#111827")
    if include_partials:
        colors = ["#2563eb", "#059669", "#dc2626", "#7c3aed", "#d97706"]
        for color, (label, key) in zip(colors, SERIES[:-1]):
            ax.plot(x_values, data[key], label=label, linewidth=1.5, color=color)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Cross section (barn)")
    ax.set_title(title)
    ax.grid(True, which="both", linestyle=":", alpha=0.35)
    if x_axis == "wavelength":
        ax.invert_xaxis()
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    return fig


def plot_cross_sections(
    dat_path: str | Path = "CRIPOOUT.DAT",
    output_path: str | Path = "CRIPO_cross_sections.png",
    include_partials: bool = True,
    x_axis: str = "energy",
) -> Path:
    """Read a `.DAT` file and generate a saved cross-section plot from it."""
    data = read_dat(dat_path)
    return plot_cross_section_data(
        data, output_path=output_path, include_partials=include_partials, x_axis=x_axis
    )


def main() -> None:
    """Command-line entry point retained for direct module execution."""
    parser = argparse.ArgumentParser(
        description="Plot CRIPO cross sections from a .DAT output file."
    )
    parser.add_argument(
        "dat_path",
        nargs="?",
        default="CRIPOOUT.DAT",
        help="Path to the CRIPO .DAT file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="CRIPO_cross_sections.png",
        help="Path to the output image file.",
    )
    parser.add_argument(
        "--total-only",
        action="store_true",
        help="Plot only the total cross section.",
    )
    parser.add_argument(
        "--x-axis",
        choices=["energy", "wavelength"],
        default="energy",
        help="Choose whether the horizontal axis is energy or neutron wavelength.",
    )
    args = parser.parse_args()

    output = plot_cross_sections(
        dat_path=args.dat_path,
        output_path=args.output,
        include_partials=not args.total_only,
        x_axis=args.x_axis,
    )
    print(f"Saved plot to {output}")


if __name__ == "__main__":
    main()
