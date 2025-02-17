#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
import numpy as np

np.set_printoptions(linewidth=np.inf)

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Exctract the final geometry of an optimization run in Gaussian. Extract vibrational frequencies if availiable""",
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="Input file name"
)

parser.add_argument("-o", type=str, default=None, help="Output file name. Default is Input_filename.xyz")

parser.add_argument("--frequencies", action='store_true', help="Extract frequencies to Input_filename_freq.dat")


class Dani:

    atomic_dict = {
        1:'H',
        6:'C',
        7:'N',
        8:'O',
    }

    def __init__(self, gaussian_output_filename):
        self.gaussian_output_filename = gaussian_output_filename
        self._file_content = None

        self._initialize()

    @classmethod
    def from_file(cls, gaussian_output_filename):
        instance = cls(gaussian_output_filename)

        return instance

    def _initialize(self):
        self._get_file_content()
        self._get_coordinates_section()

    def _get_file_content(self):
        try:
            with open(self.gaussian_output_filename, "r") as f:
                self._file_content = f.readlines()

        except FileNotFoundError:
            raise FileNotFoundError(f"Could not find file: {self.gaussian_output_filename}")

    def _get_coordinates_section(self):

        for index, line in enumerate(self._file_content[::-1]):
            if "Coordinates (Angstroms)" in line:
                self.starting_index = len(self._file_content) - index + 2
                self.end_index = self.starting_index + 300
                break
        for newindex, line in enumerate(self._file_content[self.starting_index:self.end_index]):
            if '---------------------' in line:
                self.end_index = self.starting_index + newindex
                break

        return self.starting_index, self.end_index


    def generate_xyzfile(self, filename=None):
        if filename is None:
            filename = (
                self.gaussian_output_filename.replace(".log", "").replace(".out", "") + ".xyz"
            )

        with open(filename, "w") as out_file:
            out_file.write(str(len(self._file_content[self.starting_index:self.end_index])) + "\n")

            out_file.write('Imagine a world where atoms were very small kittens\n')

            for atom in a._file_content[self.starting_index:self.end_index]:
                list_to_print = [a.atomic_dict[int(atom.strip().split()[1])]] + atom.strip().split()[-3:] + ['\n']
                out_file.write(' '.join(list_to_print))

    def extract_vibrational_modes(self, filename=None):
        intensities = []
        freqs = []
        for line in self._file_content:
            if 'Frequencies' in line:
                freqs += line.strip().split()[-3:]
            elif 'IR Inten' in line:
                intensities += line.strip().split()[-3:]


        if filename is None:
            filename = (
                self.gaussian_output_filename.replace(".log", "").replace(".out", "") + "_frequencies.dat"
        )

        with open(filename, "w") as out_file:
            for mode, frequency in enumerate(freqs):
                out_file.write(' '.join([str(mode+1), str(frequency), str(intensities[mode]),'\n']))



if __name__ == "__main__":

    # display help message if there is no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    a = Dani.from_file(args.f)
    print('coordinates', a._get_coordinates_section())
    i = a._get_coordinates_section()[0]
    f = a._get_coordinates_section()[1]
    print(i, f)
    for atom in a._file_content[i:f]:
        print([a.atomic_dict[int(atom.strip().split()[1])]] + atom.strip().split()[-3:])
    a.generate_xyzfile(args.o)
    if args.frequencies:
        a.extract_vibrational_modes()
