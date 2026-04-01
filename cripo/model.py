"""High-level OO interface for CRIPO calculations."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import pandas as pd

from .driver import CripoInput, run_cripo
from .plotter import energy_to_wavelength_angstrom, plot_cross_section_data


class CripoModel:
    """High-level object-oriented interface for running CRIPO in Python.

    The constructor takes the same physical inputs as the historical CLI. Use
    :meth:`run` to execute the calculation, :meth:`get_cross_section_data` or
    :meth:`to_dataframe` to inspect the results, and :meth:`plot_xs` to create
    plots directly from the in-memory data.
    """

    def __init__(
        self,
        name: str,
        sigma_incoherent: float,
        b_coherent_fm: float,
        sigma_thermal: float,
        debye_temperature: float,
        temperature: float,
        cell_parameter_a_angstrom: float,
        cell_parameter_c_angstrom: float,
        atomic_mass_amu: float,
        crystal_type: int,
        lowest_lethargy_exponent: int,
        highest_lethargy_exponent: int,
        points_per_lethargy_decade: int,
        electrons_z: int,
        neutron_electron_length_fm: float = 0.0013,
    ) -> None:
        """Store the CRIPO inputs for a later calculation run."""
        self.inputs = CripoInput(
            name=name,
            sigma_incoherent=sigma_incoherent,
            b_coherent_fm=b_coherent_fm,
            sigma_thermal=sigma_thermal,
            debye_temperature=debye_temperature,
            temperature=temperature,
            cell_parameter_a_angstrom=cell_parameter_a_angstrom,
            cell_parameter_c_angstrom=cell_parameter_c_angstrom,
            atomic_mass_amu=atomic_mass_amu,
            crystal_type=crystal_type,
            lowest_lethargy_exponent=lowest_lethargy_exponent,
            highest_lethargy_exponent=highest_lethargy_exponent,
            points_per_lethargy_decade=points_per_lethargy_decade,
            electrons_z=electrons_z,
            neutron_electron_length_fm=neutron_electron_length_fm,
        )
        self.result: dict[str, object] | None = None
        self.run()

    def run(self, output_dir: str | Path = ".") -> "CripoModel":
        """Execute the calculation and cache the full result on the instance."""
        self.result = run_cripo(self.inputs, output_dir=output_dir)
        return self

    def get_cross_section_data(self) -> dict[str, list[float]]:
        """Return the computed cross-section components as plain Python lists."""
        if self.result is None:
            raise RuntimeError("Run the calculation first with .run().")
        rows = self.result["rows"]
        data = {
            "energy_ev": [],
            "lethargy_mantissa": [],
            "elastic_coherent": [],
            "elastic_incoherent": [],
            "inelastic_incoherent": [],
            "inelastic_coherent": [],
            "absorption": [],
            "total": [],
        }
        for row in rows:
            for key in data:
                data[key].append(row[key])
        return data

    def to_dataframe(self) -> pd.DataFrame:
        """Return the cross-section data as a pandas ``DataFrame``."""
        data = self.get_cross_section_data()
        frame = pd.DataFrame(data)
        frame["wavelength_angstrom"] = frame["energy_ev"].map(
            energy_to_wavelength_angstrom
        )
        ordered_columns = [
            "energy_ev",
            "wavelength_angstrom",
            "lethargy_mantissa",
            "elastic_coherent",
            "elastic_incoherent",
            "inelastic_incoherent",
            "inelastic_coherent",
            "absorption",
            "total",
        ]
        return frame[ordered_columns]

    def plot_xs(
        self,
        output_path: str | Path | None = None,
        include_partials: bool = True,
        x_axis: str = "energy",
    ) -> Path:
        """Plot the computed cross sections to an image file.

        Set ``x_axis="wavelength"`` to use neutron wavelength in Angstrom
        instead of energy in eV on the horizontal axis.
        """
        data = self.get_cross_section_data()
        if output_path is None:
            safe_name = self.inputs.name.strip().replace(" ", "_") or "CRIPO"
            suffix = "wavelength" if x_axis == "wavelength" else "energy"
            output_path = f"{safe_name}_cross_sections_{suffix}.png"
        return plot_cross_section_data(
            data,
            output_path=output_path,
            include_partials=include_partials,
            title=f"{self.inputs.name} Cross Sections",
            x_axis=x_axis,
        )

    def summary(self) -> dict[str, object]:
        """Return a compact summary of the inputs, outputs, and run metadata."""
        if self.result is None:
            raise RuntimeError("Run the calculation first with .run().")
        return {
            "inputs": asdict(self.inputs),
            "fdw": self.result["fdw"],
            "idim": self.result["idim"],
            "elim": self.result["elim"],
            "output_files": self.result["output_files"],
        }
