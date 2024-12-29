#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
import numpy as np

np.set_printoptions(linewidth=np.inf)

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Exctract the final geometry of an optimization run in ORCA.""",
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="Input file name"
)
parser.add_argument("-o", type=str, default=None, help="Output file name")

parser.add_argument("--no_energies", default=True, action='store_false', help="Dont include energies in the xyzfile as comment.")


class Laura:

    def __init__(self, orca_output_filename):
        self.orca_output_filename = orca_output_filename
        self._file_content = None
        self._gradient_calculation = False
        self.ground_state_energy = None
        self.energies = None
        self._coordinates_dataframe = None

        self._initialize()

    @classmethod
    def from_file(cls, orca_output_filename):
        instance = cls(orca_output_filename)

        return instance

    def _initialize(self):
        self._get_file_content()
        self._get_nroots()
        self._is_gradient_calculation()
        self._get_gs_energy()
        self._get_excitation_energies()
        self._get_coordinates_section()
        self._extract_coordinates_to_dataframe()

    def _get_file_content(self):
        try:
            with open(self.orca_output_filename, "r") as f:
                self._file_content = f.readlines()

        except FileNotFoundError:
            raise FileNotFoundError(f"Could not find file: {self.orca_output_filename}")

    def _get_nroots(self):
        for line in self._file_content:
            if "Number of roots to be determined" in line:
                self.n_roots = int(line.strip().split()[-1])
                break

    def _is_gradient_calculation(self):
        for line in self._file_content:
            if 'SCF GRADIENT' in line:
                self._gradient_calculation = True
                break

    def _get_gs_energy(self):

        if self._gradient_calculation:
            for line in self._file_content:
                if 'Total Energy       : ' in line:
                    self.ground_state_energy = float(line.strip().split()[-4])

        else:
            for line in self._file_content:
                if "E(SCF)  =" in line:
                    print(line)
                    self.ground_state_energy = float(line.strip().split()[-2])

    def _get_excitation_energies(self):
        self.energies = np.zeros(self.n_roots + 1)
        self.energies[0] = self.ground_state_energy

        for line in self._file_content:
            if (
                "STATE" in line
                and "TD-DFT/TDA EXCITED STATES (SINGLETS)" not in line
                and "EXCITED STATE GRADIENT DONE" not in line
                and "NATURAL" not in line
            ):

                index = int(line.replace(':', ' ').strip().split()[1])
                self.energies[index] = float(line.strip().split()[3]) + self.ground_state_energy

    def _get_coordinates_section(self):

        for index, line in enumerate(self._file_content):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                self.starting_index = index + 2
            elif "CARTESIAN COORDINATES (A.U.)" in line:
                self.end_index = index - 2

        return self.starting_index, self.end_index

    def _extract_coordinates_to_dataframe(self):

        positions = np.array(
            [
                np.array(vector.strip().split()[1:]).astype(float)
                for vector in self._file_content[self.starting_index:self.end_index]
            ]
        )

        labels = [
            np.array(vector.strip().split()[0])
            for vector in self._file_content[self.starting_index:self.end_index]
        ]

        df = pd.DataFrame(
            {
                "Atom labels": labels,
                "x": positions[:, 0],
                "y": positions[:, 1],
                "z": positions[:, 2],
            }
        )

        self._coordinates_dataframe = df

    def get_coordinates_dataframe(self):

        return self._coordinates_dataframe.copy()

    def generate_xyzfile(self, filename=None, print_energies=True):
        if filename is None:
            filename = (
                self.orca_output_filename.replace(".in", "").replace(".out", "") + ".xyz"
            )

        with open(filename, "w") as out_file:
            out_file.write(str(len(self._file_content[self.starting_index:self.end_index])) + "\n")

            if print_energies:
                out_file.write(str(self.energies).replace(']', '').replace('[', '') + "\n")

            else:
                out_file.write('Imagine a world where atoms were very small kittens\n')

            for i in self._file_content[self.starting_index:self.end_index]:
                out_file.write(i.strip() + "\n")


if __name__ == "__main__":

    # display help message if there is no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    a = Laura.from_file(args.f)
    a.generate_xyzfile(args.o, args.no_energies)
