#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt

# Parser is used to input via terminal the required arguments
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


class Alberto:
    """
    Takes a excited state optimization in ORCA and plots the energy of each state at each
    step of the opitmizaiton and the actual root.

    Attributes
    ----------
    filename : str
        Name of the orca output file.
    reference_energy : float
        Value of the reference energy.
    units : str
        Units for the plot. Default is eV.


    Methods
    -------
    plot()
        Plots the energies at each step of the optimization.

    save_image(image_name)
        Saves the plot to an image.
    """

    def __init__(
        self,
        filename: str,
        reference_energy: float = 0.0,
        units: str = "eV",
    ) -> None:
        self.filename = filename
        self.reference_energy = reference_energy
        self.units = units
        self.content_list = []

        self._initialize()

    def _initialize(self):
        self._get_file_content()
        self._get_state_indexes()
        self._get_cis_energies()
        self._get_energy_array()
        self._rescale_units(self.units)

    def _get_file_content(self):
        """
        Open the file and store its content as a list

        Raises
        ------
        ValueError
            If there is no filename specified.
        FileNotFoundError
            If the file does not exist.

        Returns
        -------
        None.

        """
        if self.filename is None:
            raise ValueError("No result file was specified")

        try:
            with open(self.filename, "r") as file:
                self.content_list = file.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"File {self.filename:s} does not exist.")

    def _get_state_indexes(self):
        if self.filename is None:
            raise ValueError("No result file was found")

        counter = 0
        self.state_indexes = []

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

        self.number_of_states = len(self.state_indexes)

    def _get_cis_energies(self):
        cis_energies = []

        for line in self.content_list:
            if "E(SCF)" in line:
                cis_energies.append(float(line.strip().split()[2]))

        if self.reference_energy == 0:
            self.reference_energy = cis_energies[0]

        self.cis_array = (
            np.array(cis_energies).reshape(-1) - self.reference_energy
        )

        self.number_of_steps = len(self.cis_array)

    def _get_energy_array(self):
        total_energy_list = [[] for i in range(self.number_of_states)]

        for line in self.content_list:
            if (
                "STATE" in line
                and "TD-DFT/TDA EXCITED STATES (SINGLETS)" not in line
                and "EXCITED STATE GRADIENT DONE" not in line
            ):
                index = int(line.strip().split()[1].replace(":", "")) - 1
                energy = float(line.strip().split()[3].replace(":", ""))
                total_energy_list[index].append(energy)

        self.energy_array = np.zeros(
            [self.number_of_states, self.number_of_steps]
        )

        for pes_index, state in enumerate(total_energy_list):
            for index, energy in enumerate(state):
                self.energy_array[pes_index, index] = energy

    def set_units(self, units):
        prev_units = self.units

        if units in ["eV", "Hartree"]:
            self.units = units
        else:
            raise ValueError('Units can only be "eV" or "Hartree".')

        self._rescale_units(prev_units)

    def _rescale_units(self, prev_units):  # ??? Something is wrong, redo
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

    def set_imagename(self, image_name: str):
        """
        Sets the class parameter Alberto.imagename

        Parameters
        ----------
        image_name : str
            Image name to set.

        Returns
        -------
        None,

        """
        self.output_image = image_name

    def plot(self):
        x_range = np.arange(0, self.number_of_steps, 1)
        for state in self.energy_array:
            plt.plot(x_range, state)
        plt.plot(x_range, self.cis_array)

    def generate_image(self, output_image_name: str):

        self.plot()

        plt.savefig(output_image_name, dpi=300)


def alberto(
    filename: str,
    reference_energy: float = 0,
    output_image: str = True,
    units: str = "eV",
) -> np.array:
    # get the file and content

    file = open(filename, "r")
    cont = file.readlines()
    # grep the lines of interest
    counter = 0
    n_os = []

    for line in cont:
        if "TD-DFT/TDA EXCITED STATES (SINGLETS)" in line:
            counter += 1
            if counter >= 2:
                break
        elif (
            "STATE" in line
            and "TD-DFT/TDA EXCITED STATES (SINGLETS)" not in line
            and "EXCITED STATE GRADIENT DONE" not in line
        ):
            n_os.append(int(line.strip().split()[1].replace(":", "")))

    total_list = [[] for i in range(0, max(n_os))]
    cis_energies = []

    for line in cont:
        if (
            "STATE" in line
            and "TD-DFT/TDA EXCITED STATES (SINGLETS)" not in line
            and "EXCITED STATE GRADIENT DONE" not in line
        ):
            index = int(line.strip().split()[1].replace(":", "")) - 1
            energy = float(line.strip().split()[3].replace(":", ""))
            total_list[index].append(energy)
        elif "E(SCF)" in line:
            cis_energies.append(float(line.strip().split()[2]))

    # print(total_list)

    cis_array = np.array(cis_energies) - reference_energy
    total_arrays = [np.array(l) for l in total_list]

    if units == "eV":
        cis_array *= 27.2114
        for i in total_arrays:
            i *= 27.2114
    else:
        units = "Hartree"

    x = np.arange(1, len(total_list[0]) + 1, 1)

    plt.plot(x, cis_array, label="Ground state")

    # obtain the state of the actual root
    curr_energy = []
    curr_energy_index = []
    for line in cont:
        if "DE(CIS) =" in line:
            curr_energy_index.append(
                float(line.strip().split()[-1].replace(")", ""))
            )

    for i, root in enumerate(curr_energy_index):
        curr_energy.append(total_arrays[int(root) - 1][i])

    if output_image is True:
        output_image = filename

    for i in range(len(total_arrays)):
        plt.plot(x, cis_array + total_arrays[i], label="Root %i" % (i + 1))

    plt.scatter(
        x,
        cis_array + curr_energy,
        label="Active root",
        marker="x",
        c="rebeccapurple",
    )
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.xlabel("Step")
    plt.ylabel("Energy / %s" % units)
    plt.title(output_image.replace(".in.out", ""))
    plt.savefig(
        output_image.replace(".in.out", ".jpg"), dpi=250, bbox_inches="tight"
    )
    # plt.show()

    to_return = [np.array(cis_array + array) for array in total_arrays]

    to_return.append(cis_array + curr_energy)

    return to_return


if __name__ == "__main__":
    # if len(sys.argv) == 1:
    #     parser.print_help(sys.stderr)
    #     sys.exit(1)

    # args = parser.parse_args()

    a = Alberto("upgrades/beta-carboline_opt_followiroot_3.in.out")
