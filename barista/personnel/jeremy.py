#!/usr/bin/env python3

from functools import cached_property

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike


class Jeremy:

    # =========================================================================
    #     Magic Methods
    # =========================================================================

    def __init__(self, xyzfile, bond_thresh: float = 1.6):
        self.xyzfile = xyzfile

        self._position_array = None

        self._custom_connectivity = None
        self.bond_thresh = bond_thresh

        self._read_xyzfile()
        self._extract_atoms()
        self._extract_positions()
        self._extend_labels()

    # =========================================================================
    #     Public interface methods
    # =========================================================================

    def override_connectivity_matrix(self, connectivity_matrix: ArrayLike):
        """
        Override connectivity matrix to generate internals from a reference.

        Redefines the internal coordinates based in the reference connectivity.
        Useful when comparing two structures, being able to determine their
        'likeness' by obtaining the same internals of both and their values.

        Parameters
        ----------
        connectivity_matrix : ArrayLike
            Connectivity matrix of the reference system.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if connectivity_matrix.shape != self.connectivity_matrix.shape:
            raise ValueError("The dimensions of both molecules are not equal")

    def clear_cache(self):

        cached_properties = [
            "xyz_df",
            "connectivity_matrix",
            "distance_matrix",
            "r_vector_matrix",
            "bond_list",
            "angle_list",
            "dihedral_list",
            "bond_values",
            "angle_values",
            "dihedral_values",
        ]

        for prop in cached_properties:
            if prop in self.__dict__:
                del self.__dict__[prop]

    # =========================================================================
    #     Property methods
    # =========================================================================
    @property
    def position_array(self) -> ArrayLike:
        """
        Return a copy of the xyz cooridnates array.

        Returns
        -------
        ArrayLike

        """
        return np.copy(self._position_array)

    @property
    def atom_labels(self):
        return self._atom_labels.copy()

    @property
    def internal_list(self) -> list:
        """
        Return copy of the internal coordinate list.

        Returns
        -------
        list

        """
        return self.bond_list + self.angle_list + self.dihedral_list

    @property
    def internal_values(self) -> list:
        """
        Return copy of the list containing internal coordinate values.

        Returns
        -------
        list

        """
        return self.bond_values + self.angle_values + self.dihedral_values

    @cached_property
    def connectivity_matrix(self):
        """
        Build connectivity matrix.

        If below the threshold, will fill with 1. If no connection below the
        threshold, it will search for the nearest atom and assign the value 10
        in order to leave no atom unconnected.


        Parameters
        ----------
        bond_thresh : float, optional
            Bond threshold. The default is 1.6.

        Returns
        -------
        connectivity_matrix : ArrayLike
            Connectivity matrix of the molecule.

        """
        if self._custom_connectivity is not None:
            return np.copy(self._custom_connectivity)

        _connectivity_matrix = np.zeros([self._n_atoms, self._n_atoms])

        for i in range(0, self._n_atoms):
            for j in range(i, self._n_atoms):
                r_ij = np.linalg.norm(
                    self._position_array[i] - self._position_array[j]
                )
                if r_ij < self.bond_thresh and i != j:
                    _connectivity_matrix[i][j] = _connectivity_matrix[j][i] = 1

        for i, row in enumerate(_connectivity_matrix):
            if sum(row) == 0:

                distances = sorted(self.distance_matrix[i], reverse=False)

                j = np.where(self.distance_matrix[i] == distances[1])

                _connectivity_matrix[i][j[0][0]] = _connectivity_matrix[
                    j[0][0]
                ][i] = 10
        return np.copy(_connectivity_matrix)

    @cached_property
    def xyz_df(self) -> pd.DataFrame:
        """
        Return copy of the cartesian dataframe.

        Returns
        -------
        pd.DataFrame

        """
        position_data = {
            "Atom": self.atom_labels,
            "x": self.position_array[:, 0],
            "y": self.position_array[:, 1],
            "z": self.position_array[:, 2],
        }
        _xyz_df = pd.DataFrame.from_dict(position_data)

        return _xyz_df.copy()

    @cached_property
    def distance_matrix(self):
        _distance_matrix = np.zeros([self._n_atoms, self._n_atoms])
        for i in range(0, self._n_atoms):
            for j in range(i, self._n_atoms):
                r_ij = np.linalg.norm(
                    self._position_array[i] - self._position_array[j]
                )
                _distance_matrix[i][j] = _distance_matrix[j][i] = r_ij
        return np.copy(_distance_matrix)

    @cached_property
    def r_vector_matrix(self):
        _r_vector_matrix = []

        for i in range(self._n_atoms):
            append_list = []
            for j in range(self._n_atoms):
                append_list.append(
                    -(self._position_array[j] - self._position_array[i])
                )
            _r_vector_matrix.append(append_list)

        return np.copy(np.array(_r_vector_matrix))

    # Build internals

    @cached_property
    def bond_list(self):

        _bond_list = []

        for i in range(self._n_atoms):
            for j in range(i, self._n_atoms):
                if (
                    self.connectivity_matrix[i][j] == 1
                    or self.connectivity_matrix[i][j] == 10
                ):
                    _bond_list.append([i, j])

        return _bond_list.copy()

    @cached_property
    def angle_list(self) -> list:
        """
        Build the angle list.

        Done by adding an atom to the end of the bonds forwards and backwards.
        If the angle is new it is stored. If not, discarded.

        Returns
        -------
        list
            Angle list with atomic indices.

        """
        _angle_list = []
        angle_raw_list = []  # placeholder list for bonds that can be repeated

        for _ in range(2):  # we iterate twice, once forward and once backwards
            self.bond_list = [
                bond[::-1] for bond in self.bond_list
            ]  # reverse bond indexes
            for bond in self.bond_list:
                for atom in range(len(self.bond_list)):
                    if (
                        self.connectivity_matrix[bond[-1]][atom] in [1, 10]
                        and atom not in bond
                    ):  # Check connectivity
                        angle_raw_list.append([bond[0], bond[1], atom])

        for angle in angle_raw_list:
            # check uniqueness (avoid repetition of the same angle indexes flipped)
            if angle not in _angle_list and angle[::-1] not in _angle_list:
                _angle_list.append(angle)

        _angle_list.sort(key=lambda x: x[0])  # sort the list by the first index

        return _angle_list.copy()

    @cached_property
    def dihedral_list(self):
        """Build the dihedral list.

        It is done by adding an atom to the end of the angles forwards and
        backwards. If the dihedral is new it is stored. If not, discarded

        """
        dihedral_raw_list = []
        _dihedral_list = []

        for _ in range(2):  # same as before, two iterations
            self.angle_list = [
                angle[::-1] for angle in self.angle_list
            ]  # reverse angle indexes
            for angle in self.angle_list:
                for atom in range(self._n_atoms):
                    if (
                        self.connectivity_matrix[angle[-1]][atom] == 1
                        and atom not in angle
                    ):  # Check connectivity
                        dihedral_raw_list.append(
                            [angle[0], angle[1], angle[2], atom]
                        )

        for dihedral in dihedral_raw_list:
            if (
                dihedral not in _dihedral_list
                and dihedral[::-1] not in _dihedral_list
            ):  # check uniqueness
                _dihedral_list.append(dihedral)

        _dihedral_list.sort(
            key=lambda x: x[0]
        )  # sort the list by the first index

        return _dihedral_list

    # Calculate internals

    @cached_property
    def bond_values(self):
        _bond_values = []

        for bond in self.bond_list:
            i, j = bond
            _bond_values.append(self.distance_matrix[i][j])

        return _bond_values.copy()

    @cached_property
    def angle_values(self):
        """Calculate the angles.

        Uses theta = arccos((r_ba · r_bc)/(||r_ba||*||r_bc||)) * 180/pi.

        """
        angle_values = []

        for angle in self.angle_list:
            a, b, c = angle
            angle_values.append(
                float(
                    np.arccos(
                        (
                            self.r_vector_matrix[b][a]
                            @ self.r_vector_matrix[b][c]
                        )
                        / (
                            self.distance_matrix[a][b]
                            * self.distance_matrix[b][c]
                        )
                    )
                    * 180
                    / np.pi
                )
            )

        return angle_values.copy()

    @cached_property
    def dihedral_values(self):
        """Calculate the dihedral angles.

        Done with the atan(sin_phi, cos_phi) formula.

        """
        dih_val = []

        for dihedral in self.dihedral_list:
            a, b, c, d = dihedral
            t = np.cross(self.r_vector_matrix[a][b], self.r_vector_matrix[b][c])
            u = np.cross(self.r_vector_matrix[b][c], self.r_vector_matrix[c][d])
            v = np.cross(t, u)

            cos_phi = t @ u / (np.linalg.norm(t) * np.linalg.norm(u))
            sin_phi = (self.r_vector_matrix[b][c] @ v) / (
                np.linalg.norm(t)
                * np.linalg.norm(u)
                * np.linalg.norm(self.r_vector_matrix[b][c])
            )

            phi = np.arctan2(sin_phi, cos_phi) * 180 / np.pi

            dih_val.append(phi)

        return dih_val.copy()

    # =========================================================================
    #     Internal computation preparation methods
    # =========================================================================

    def _read_xyzfile(self):
        try:
            with open(self.xyzfile, "r") as f:
                self._file_content_list = f.readlines()[2:]
            self._n_atoms = len(self._file_content_list)

        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f'"{self.xyzfile}" geometry file was not found'
            ) from exc

    def _extract_atoms(self):
        self._atom_labels = []
        for atom in self._file_content_list:
            self._atom_labels.append(atom.strip().split()[0].capitalize())

    def _extract_positions(self):
        self._position_array = np.zeros([self._n_atoms, 3], dtype=float)

        for i, atom in enumerate(self._file_content_list):
            self._position_array[i] = np.array(
                atom.strip().split()[1:], dtype=float
            )

    def _extend_labels(self):
        counter_labels = list(set(self._atom_labels))
        counter = np.zeros(len(counter_labels))

        extended_labels = []

        for label in self._atom_labels:
            counter[counter_labels.index(label)] += 1
            extended_labels.append(
                label + str(int(counter[counter_labels.index(label)]))
            )

        self._atom_labels = extended_labels


if __name__ == "__main__":
    co2 = Jeremy("updates/asdf.xyz")
    co2_distorted = Jeremy("updates/jkl.xyz")
