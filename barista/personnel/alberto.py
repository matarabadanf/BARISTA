#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(
    description="""Takes a excited state optimization in ORCA and plots the energy of each state at each step of the opitmizaiton and the actual root.""",
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="ORCA optimization file"
)
parser.add_argument("-o", type=str, default=True, help="Output image filename")
parser.add_argument(
    "-en",
    type=float,
    default=0,
    help="Reference energy of the ground state optimized energy",
)
parser.add_argument("-u", type=str, default="eV", help="Energy units")

# display help message if there is no arguments

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional, Union, Tuple


class Alberto:
    """
    Takes an excited state optimization in ORCA and plots the energy of each state at each
    step of the optimization and the actual root.

    This class parses ORCA output files from TD-DFT/TDA calculations to track and
    visualize the energy evolution of electronic states during geometry optimization.

    Attributes
    ----------
    filename : str
        Name of the ORCA output file.
    reference_energy : float
        Value of the reference energy (typically ground state energy) in Hartree.
    units : str
        Units for the energy values in plots. Default is "eV".
    content_list : List[str]
        Content of the ORCA output file as a list of lines.
    state_indexes : List[int]
        List of state indexes found in the file.
    number_of_states : int
        Total number of excited states found.
    number_of_steps : int
        Number of optimization steps.
    cis_array : np.ndarray
        Array containing SCF energies at each step.
    energy_array : np.ndarray
        2D array containing energy values for each state at each step.
    bulk_energy_array : np.ndarray
        Copy of energy_array before reference shift.
    curr_energy : np.ndarray
        Energy of the current root at each optimization step.
    _starting_root : int
        Initial root used for the optimization.
    _final_root : int
        Final root after optimization.

    Methods
    -------
    _initialize() -> None
        Initialize all data structures and parse the output file.
    _get_file_content() -> None
        Open and read the output file.
    _get_state_indexes() -> None
        Extract state indexes from the output file.
    _get_cis_energies() -> None
        Extract SCF energies from the output file.
    _get_abs_energy_array() -> None
        Calculate absolute energies for all states.
    _get_relative_energy() -> None
        Calculate energies relative to the reference energy.
    _get_current_energy() -> None
        Extract the current energy surface at each step.
    set_units(units: str) -> None
        Change the energy units used for plotting.
    _rescale_units(prev_units: str) -> None
        Rescale energy values according to the selected units.
    set_imagename(image_name: str) -> None
        Set the output image name.
    _prepare_plot() -> Tuple[plt.Figure, plt.Axes]
        Prepare the energy plot.
    plot() -> None
        Display the energy plot.
    generate_image(output_image_name: str) -> None
        Save the energy plot to an image file.
    """

    def __init__(
        self,
        filename: str,
        reference_energy: float = 0.0,
        units: str = "eV",
    ) -> None:
        """
        Initialize the Alberto class with the specified parameters.

        Parameters
        ----------
        filename : str
            Path to the ORCA output file to analyze.
        reference_energy : float, optional
            Reference energy value in Hartree, by default 0.0.
            Used to calculate relative energies.
        units : str, optional
            Energy units for plotting ("eV" or "Hartree"), by default "eV".
        """
        self.filename = filename
        self.reference_energy = reference_energy
        self.units = units
        self.content_list: List[str] = []
        self._starting_root: Optional[int] = None
        self._final_root: Optional[int] = None

        self._initialize()

    def _initialize(self) -> None:
        """
        Initialize all data structures and run the parsing methods in sequence.
        
        This method coordinates the full initialization process by calling all
        the necessary parsing and processing methods in the correct order.
        """
        self._get_file_content()
        self._get_state_indexes()
        self._get_cis_energies()
        self._get_abs_energy_array()
        self._get_relative_energy()
        self._rescale_units(self.units)
        self._get_current_energy()

    def _get_file_content(self) -> None:
        """
        Open the ORCA output file and store its content as a list of strings.

        Raises
        ------
        ValueError
            If there is no filename specified.
        FileNotFoundError
            If the file does not exist.
        """
        if self.filename is None:
            raise ValueError("No result file was specified")

        try:
            with open(self.filename, "r") as file:
                self.content_list = file.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"File {self.filename:s} does not exist.")

    def _get_state_indexes(self) -> None:
        """
        Extract state indexes from the ORCA output file.
        
        Finds all states mentioned in the TD-DFT/TDA section and stores their indices.
        Sets the number_of_states attribute based on discovered states.
        
        Raises
        ------
        ValueError
            If no result file was found.
        """
        if self.filename is None:
            raise ValueError("No result file was found")

        counter = 0
        self.state_indexes: List[int] = []

        for line in self.content_list:
            if "TD-DFT/TDA EXCITED STATES (SINGLETS)" in line:
                counter += 1
                if counter >= 2:
                    break
            elif (
                "STATE" in line
                and "TD-DFT/TDA EXCITED STATES (SINGLETS)" not in line
                and "EXCITED STATE GRADIENT DONE" not in line
            ):
                self.state_indexes.append(
                    int(line.strip().split()[1].replace(":", ""))
                )

        self.number_of_states: int = len(self.state_indexes)

    def _get_cis_energies(self) -> None:
        """
        Extract SCF energies from the ORCA output file.
        
        Collects all SCF energies from the file and stores them in cis_array.
        If no reference energy was provided, the first SCF energy is used.
        Sets the number_of_steps attribute based on collected energies.
        """
        cis_energies: List[float] = []

        for line in self.content_list:
            if "E(SCF)" in line:
                cis_energies.append(float(line.strip().split()[2]))

        if self.reference_energy == 0:
            self.reference_energy = cis_energies[0]

        self.cis_array: np.ndarray = np.array(cis_energies).reshape(-1)
        #     - self.reference_energy
        # )

        self.number_of_steps: int = len(self.cis_array)

    def _get_abs_energy_array(self) -> None:
        """
        Parse the ORCA output file to obtain the absolute energies in a.u.

        The procedure performed is the following:
         - Extract the energy increase from the current s0 to the excited states
         - Set the first row of the energy array to the SCF energy
         - Add the energy of the excited states to the SCF in each state

        This returns the total absolute energy of each of the states
        (typically in the test energies between -562.25 and -561.85 a.u.).
        """
        total_energy_list: List[List[float]] = [[] for i in range(self.number_of_states + 1)]

        for line in self.content_list:
            if (
                "STATE" in line
                and "TD-DFT/TDA EXCITED STATES (SINGLETS)" not in line
                and "EXCITED STATE GRADIENT DONE" not in line
            ):
                index = int(line.strip().split()[1].replace(":", "")) - 1
                energy = float(line.strip().split()[3].replace(":", ""))
                total_energy_list[index + 1].append(energy)

        self.energy_array: np.ndarray = np.zeros(
            [self.number_of_states + 1, self.number_of_steps]
        )

        self.energy_array[0] = np.copy(self.cis_array)

        for pes_index, state in enumerate(total_energy_list):
            for index, energy in enumerate(state):
                self.energy_array[pes_index, index] = energy

        self.energy_array[1:] += self.energy_array[0]

        self.bulk_energy_array: np.ndarray = np.copy(self.energy_array)

    def _get_relative_energy(self) -> None:
        """
        Calculate energies relative to the reference energy.
        
        Subtracts the reference energy from all energy values to obtain
        relative energies that are more suitable for visualization.
        """
        if self.reference_energy == 0:
            self.reference_energy = self.cis_array[0]

        self.energy_array = (
            np.copy(self.bulk_energy_array) - self.reference_energy
        )

    def _get_current_energy(self) -> None:
        """
        Extract the current energy surface at each optimization step.
        
        Identifies which electronic state is being followed during the optimization
        at each step and extracts its energy. Also determines the starting and final 
        roots used in the optimization.
        """
        curr_energy_index: np.ndarray = np.zeros(self.number_of_steps, dtype=int)
        counter = 0

        for line in self.content_list:
            if "DE(CIS) =" in line:
                curr_energy_index[counter] = int(
                    line.strip().split()[-1].replace(")", "")
                )
                counter += 1
        self._starting_root = curr_energy_index[0]
        self._final_root = curr_energy_index[-1]

        self.curr_energy: np.ndarray = np.zeros(self.number_of_steps)

        for i, index in enumerate(curr_energy_index):
            self.curr_energy[i] = self.energy_array[index][i]

    def set_units(self, units: str) -> None:
        """
        Change the energy units used for plotting.
        
        Parameters
        ----------
        units : str
            Energy units to use ("eV" or "Hartree").
            
        Raises
        ------
        ValueError
            If an unsupported unit is specified.
        """
        prev_units = self.units

        if units in ["eV", "Hartree"]:
            self.units = units
        else:
            raise ValueError('Units can only be "eV" or "Hartree".')

        self._rescale_units(prev_units)

    def _rescale_units(self, prev_units: str) -> None:
        """
        Rescale energy values according to the selected units.
        
        Converts energy values between eV and Hartree using the conversion
        factor 27.2114 eV/Hartree.
        
        Parameters
        ----------
        prev_units : str
            Previous units used for energy values.
        """
        if (
            self.units == prev_units
            and self.units == "eV"
            or self.units != prev_units
            and self.units == "eV"
        ):
            self.cis_array *= 27.2114
            self.energy_array *= 27.2114

        elif (
            self.units == prev_units
            and self.units == "Hartree"
            or self.units != prev_units
            and self.units == "Hartree"
        ):
            self.cis_array /= 27.2114
            self.energy_array /= 27.2114

    def set_imagename(self, image_name: str) -> None:
        """
        Sets the output image name for saving plots.

        Parameters
        ----------
        image_name : str
            Image name to set.
        """
        self.output_image: str = image_name

    def _prepare_plot(self) -> Tuple[plt.Figure, plt.Axes]:
        """
        Prepare the energy plot showing all states and the current surface.
        
        Creates a figure with plots for all electronic states and marks the
        current surface being followed during optimization.
        
        Returns
        -------
        Tuple[plt.Figure, plt.Axes]
            The created figure and axes objects.
        """
        fig, ax = plt.subplots()

        x_range = np.arange(0, self.number_of_steps, 1)
        for surface, state in enumerate(self.energy_array):
            ax.plot(x_range, state, label=f"$S{surface}$")

        ax.scatter(
            x_range,
            self.curr_energy,
            label="Current surface",
            c="rebeccapurple",
            marker="x",
        )

        ax.set_ylabel(f"Energy difference / {self.units}")
        ax.set_xlabel("Step")
        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        fig.set_size_inches(8, 4)

        return fig, ax

    def plot(self) -> None:
        """
        Display the energy plot showing all states and the current surface.
        
        Creates and displays a plot showing the energy evolution of all states
        during the optimization, highlighting the current surface being followed.
        """
        fig, ax = self._prepare_plot()

    def generate_image(self, output_image_name: str) -> None:
        """
        Save the energy plot to an image file.
        
        Parameters
        ----------
        output_image_name : str
            Name of the output image file.
        """
        self.plot()
        plt.savefig(output_image_name, dpi=300)

    @property
    def starting_root(self) -> int:
        """
        Get the initial root used for the optimization.
        
        Returns
        -------
        int
            Index of the starting root.
        """
        return np.copy(self._starting_root)

    @property
    def final_root(self) -> int:
        """
        Get the final root after optimization.
        
        Returns
        -------
        int
            Index of the final root.
        """
        return np.copy(self._final_root)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    a = Alberto(args.f, args.en, args.u)

    print(a.starting_root, a.final_root)

    if args.o is not True:
        a.generate_image(args.o)

    # tests a = Alberto("tests/xanthine_opt_followiroot_6.in.out", 0, "eV")
