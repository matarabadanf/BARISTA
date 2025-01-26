#!/usr/bin/env python3

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike


class Jeremy:

    # =========================================================================
    #     Magic Methods
    # =========================================================================

    def __init__(self, xyzfile):
        self.xyzfile = xyzfile

        self._position_array = None

        self._connectivity_matrix = None

        self._bond_list = None
        self._angle_list = None
        self._dihedral_list = None

        self._bond_values = None
        self._angle_values = None
        self._dihedral_values = None

        self._read_xyzfile()
        self._extract_atoms()
        self._extract_positions()
        self._extend_labels()
        self._build_xyz_dataframe()
        self._build_distance_matrix()
        self._calculate_r_vector_matrix()

    # Public interface methods
    def build_connectivity_matrix(self, bond_thresh: float = 1.6):
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
        self._connectivity_matrix = np.zeros([self._n_atoms, self._n_atoms])

        for i in range(0, self._n_atoms):
            for j in range(i, self._n_atoms):
                r_ij = np.linalg.norm(
                    self._position_array[i] - self._position_array[j]
                )
                if r_ij < bond_thresh and i != j:
                    self._connectivity_matrix[i][j] = self._connectivity_matrix[
                        j
                    ][i] = 1

        for i, row in enumerate(self._connectivity_matrix):
            if sum(row) == 0:

                distances = sorted(self.distance_matrix[i], reverse=False)

                j = np.where(self.distance_matrix[i] == distances[1])

                self._connectivity_matrix[i][j[0][0]] = (
                    self._connectivity_matrix[j[0][0]][i]
                ) = 10

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
        if connectivity_matrix.shape != self._connectivity_matrix.shape:
            raise ValueError("The dimensions of both molecules are not equal")

        self._connectivity_matrix = connectivity_matrix
        self._get_bond_list()
        self._build_angles()
        self._build_dihedrals()

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
    def connectivity_matrix(self) -> ArrayLike:
        """
        Return a copy of the system connectivity matrix.

        Returns
        -------
        ArrayLike
            DESCRIPTION.

        """
        if self._connectivity_matrix is None:
            self.build_connectivity_matrix()

        return np.copy(self._connectivity_matrix)

    @property
    def xyz_df(self) -> pd.DataFrame:
        """
        Return copy of the cartesian dataframe.

        Returns
        -------
        pd.DataFrame

        """
        if self._xyz_df is None:
            self._build_xyz_dataframe()

        return self._xyz_df.copy()

    @property
    def bond_list(self) -> list:
        """
        Return copy of the bond list.

        Returns
        -------
        list

        """
        if self._bond_list is None:
            self._get_bond_list()
        return self._bond_list.copy()

    @property
    def angle_list(self) -> list:
        """
        Return copy of the angle list.

        Returns
        -------
        list

        """
        if self._angle_list is None:
            self._build_angles()
        return self._angle_list.copy()

    @property
    def dihedral_list(self) -> list:
        """
        Return copy of the dihedral list.

        Returns
        -------
        list

        """
        if self._dihedral_list is None:
            self._build_dihedrals()
        return self._dihedral_list.copy()

    @property
    def internal_list(self) -> list:
        """
        Return copy of the internal coordinate list.

        Returns
        -------
        list

        """
        return self._build_internal_coordinates().copy()

    @property
    def internal_values(self) -> list:
        """
        Return copy of the list containing internal coordinate values.

        Returns
        -------
        list

        """
        return np.array(self._calculate_internals().copy())

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
        self.atom_labels = []
        for atom in self._file_content_list:
            self.atom_labels.append(atom.strip().split()[0].capitalize())

    def _extract_positions(self):
        self._position_array = np.zeros([self._n_atoms, 3], dtype=float)

        for i, atom in enumerate(self._file_content_list):
            self._position_array[i] = np.array(
                atom.strip().split()[1:], dtype=float
            )

    def _build_distance_matrix(self):
        self.distance_matrix = np.zeros([self._n_atoms, self._n_atoms])
        for i in range(0, self._n_atoms):
            for j in range(i, self._n_atoms):
                r_ij = np.linalg.norm(
                    self._position_array[i] - self._position_array[j]
                )
                self.distance_matrix[i][j] = self.distance_matrix[j][i] = r_ij

    def _calculate_r_vector_matrix(self):
        self.r_vector_matrix = []

        for i in range(self._n_atoms):
            append_list = []
            for j in range(self._n_atoms):
                append_list.append(
                    -(self._position_array[j] - self._position_array[i])
                )
            self.r_vector_matrix.append(append_list)

    def _extend_labels(self):
        counter_labels = list(set(self.atom_labels))
        counter = np.zeros(len(counter_labels))

        extended_labels = []

        for label in self.atom_labels:
            counter[counter_labels.index(label)] += 1
            extended_labels.append(
                label + str(int(counter[counter_labels.index(label)]))
            )

        self.atom_labels = extended_labels

    def _build_xyz_dataframe(self):
        position_data = {
            "Atom": self.atom_labels,
            "x": self._position_array[:, 0],
            "y": self._position_array[:, 1],
            "z": self._position_array[:, 2],
        }
        self._xyz_df = pd.DataFrame.from_dict(position_data)

    # =========================================================================
    #     Computations
    # =========================================================================

    def _build_internal_coordinates(self):

        self._calculate_r_vector_matrix()

        return self.bond_list + self.angle_list + self.dihedral_list

    def _calculate_internals(self):
        return (
            self._calculate_bonds()
            + self._calculate_angles()
            + self._calculate_dihedrals()
        )

    # =========================================================================
    # Build internals
    # =========================================================================

    def _get_bond_list(self):

        if self._connectivity_matrix is None:
            self.build_connectivity_matrix()

        self._bond_list = []

        for i in range(self._n_atoms):
            for j in range(i, self._n_atoms):
                if (
                    self._connectivity_matrix[i][j] == 1
                    or self._connectivity_matrix[i][j] == 10
                ):
                    self._bond_list.append([i, j])

    def _build_angles(self):
        """Build the angle list.

        Done by adding an atom to the end of the bonds forwards and backwards.
        If the angle is new it is stored. If not, discarded

        """
        if self._bond_list is None:
            self._get_bond_list()

        self._angle_list = []
        angle_raw_list = []  # placeholder list for bonds that can be repeated

        for _ in range(2):  # we iterate twice, once forward and once backwards
            self._bond_list = [
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
            if (
                angle not in self._angle_list
                and angle[::-1] not in self._angle_list
            ):
                self._angle_list.append(angle)

        self._angle_list.sort(
            key=lambda x: x[0]
        )  # sort the list by the first index

    def _build_dihedrals(self):
        """Build the dihedral list.

        It is done by adding an atom to the end of the angles forwards and
        backwards. If the dihedral is new it is stored. If not, discarded

        """
        dihedral_raw_list = []
        self._dihedral_list = []

        for _ in range(2):  # same as before, two iterations
            self._angle_list = [
                angle[::-1] for angle in self._angle_list
            ]  # reverse angle indexes
            for angle in self._angle_list:
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
                dihedral not in self._dihedral_list
                and dihedral[::-1] not in self._dihedral_list
            ):  # check uniqueness
                self._dihedral_list.append(dihedral)

        self._dihedral_list.sort(
            key=lambda x: x[0]
        )  # sort the list by the first index

        return self._dihedral_list

    # =========================================================================
    #  Calculate internals
    # =========================================================================

    def _calculate_bonds(self):
        bond_values = []

        for bond in self.bond_list:
            i, j = bond
            bond_values.append(self.distance_matrix[i][j])

        self._bond_values = bond_values

        return self._bond_values.copy()

    def _calculate_angles(self):
        """Calculate the angles.

        Uses theta = arccos((r_ba · r_bc)/(||r_ba||*||r_bc||)) * 180/pi.

        """
        if self._angle_values is None:
            self._build_internal_coordinates()

        angle_values = []

        for angle in self._angle_list:
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

        self._angle_values = angle_values

        return self._angle_values.copy()

    def _calculate_dihedrals(self):
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

        self._dihedral_values = dih_val

        return self._dihedral_values.copy()


if __name__ == "__main__":
    co2 = Jeremy("updates/asdf.xyz")
    co2_distorted = Jeremy("updates/jkl.xyz")
