#!/usr/bin/env python3

import numpy as np
import pandas as pd


class Jeremy:

    def __init__(self, xyzfile):
        self.xyzfile = xyzfile
        self._initialize()

    def _initialize(self):
        self._read_xyzfile()
        self._extract_atoms()

    def _read_xyzfile(self):
        try:
            with open(self.xyzfile, "r") as f:
                self._file_content_list = f.readlines()[2:]

            self._n_atoms = len(self._file_content_list) - 2

        except FileNotFoundError:
            raise FileNotFoundError(
                f'"{self.xyzfile}" geometry file was not found'
            )

    def _extract_atoms(self):
        self._atom_list = []
        for atom in self._file_content_list:
            self._atom_list.append(atom.strip().split()[0])

    def _extract_positions(self):
        self.position_array = np.zeros([self._n_atoms, 3], dtype=float)

    def _build_xyz_dataframe(self):
        pass


if __name__ == "__main__":
    J = Jeremy("upgrades/asdf.xyz")
