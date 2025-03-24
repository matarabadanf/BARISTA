#!/usr/bin/env python3
import inspect
from functools import cached_property, lru_cache
import argparse
import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
import sys

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Extract atoms from an xyz file with a center point and a radius."""
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="Original XYZ file."
)

parser.add_argument("-o", type=str, default='None', help="Output xyz filename.")

parser.add_argument(
    "-r",
    type=float,
    default=1.,
    help="Sphere radius.",
)

parser.add_argument("-p", type=str, default="0 0 0", help="Central sphere point.")

class Jeremy:

    # =========================================================================
    #     Special Methods
    # =========================================================================

    def __init__(self, xyzfile, bond_thresh: float = 1.6):
        self.xyzfile = xyzfile

        self._position_array = None

        self._custom_connectivity = None
        self.bond_thresh = bond_thresh

        self._read_xyzfile()
        self._extract_atoms()
        self._extract_positions()
        # self._extend_labels()

    # =========================================================================
    #     Public interface methods
    # =========================================================================

    def override_connectivity_matrix(self, cm: ArrayLike) -> None:
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
        if cm.shape != self.connectivity_matrix.shape:
            raise ValueError("The dimensions of both molecules are not equal")

        self._custom_connectivity = cm

        self.clear_cache()

        _ = self.dihedral_list
        _ = self.angle_list
        _ = self.bond_list

    def clear_cache(self) -> None:
        """
        Clear cached properties.

        Returns
        -------
        None

        """
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
            "bond_list",
            "angle_list",
            "dihedral_list",
            "dihedral_values",
            "internal_list",
            "internal_values",
        ]

        for prop in cached_properties:
            if prop in self.__dict__:
                del self.__dict__[prop]

    def rmsdiff(self, reference):
        coords_reference = reference._position_array
        self_coords = self._position_array

        sq_diff = (coords_reference - self_coords)**2
        mean_square_difference = np.mean(sq_diff)
        
        return np.sqrt(mean_square_difference)

    def compare_internals(self, reference:'Jeremy') -> np.ndarray:

        # we will enforce the overwrite of the connectivity matriux in order to 
        # ensure that the internal coordinate types remain consistent

        reference.override_connectivity_matrix(self.connectivity_matrix)

        # separate coordinates in lengths and angles

        n_bonds = len([1 for i in self.internal_list if len(i) == 2])

        print(f'The number of bonds of these molecules are: {n_bonds}')

        bond_diff = self.internal_values[:n_bonds] - reference.internal_values[:n_bonds]
        print(bond_diff)

        angle_pairs = [np.array([a, b]) for a, b in zip(self.internal_values[n_bonds:], reference.internal_values[n_bonds:])]

        print(angle_pairs)

        angle_diff = np.zeros(len(angle_pairs))

        for index, pair in enumerate(angle_pairs):
            if np.all(pair < 0) or np.all(pair > 0):

                angle_diff[index] = abs(abs(max(pair)) - abs(min(pair))) % 360
                print(pair, angle_diff[index])
            else:
                pos = max(pair)
                neg = min(pair) + 360

                print(f'Postive is {pos:5.3f}, negative is {neg:5.3f}')
                angle_diff[index] = abs(abs(min([pos,neg])) - abs(max([pos,neg]))) % 360 if neg > 180 else 0 
                print(f'discrepant pair: {pair}, {angle_diff[index]:5.3f}')

        # for angle, value in enumerate(angle_diff):
        #     print(f'Angle deviation of {angle+n_bonds}, {str(self.internal_list[n_bonds+angle]):>16s} is: {value:6.3f}')
 
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
    def atom_labels(self) -> list:
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
        return np.array(
            self.bond_values + self.angle_values + self.dihedral_values
        )

    @cached_property
    def connectivity_matrix(self) -> ArrayLike:
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
    def distance_matrix(self) -> ArrayLike:
        _distance_matrix = np.zeros([self._n_atoms, self._n_atoms])
        for i in range(0, self._n_atoms):
            for j in range(i, self._n_atoms):
                r_ij = np.linalg.norm(
                    self._position_array[i] - self._position_array[j]
                )
                _distance_matrix[i][j] = _distance_matrix[j][i] = r_ij
        return np.copy(_distance_matrix)

    @cached_property
    def r_vector_matrix(self) -> ArrayLike:
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
    def bond_list(self) -> list:

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
                for atom in range(len(self.connectivity_matrix)):
                    if (
                        self.connectivity_matrix[bond[-1]][atom] in [1, 10]
                        and atom not in bond
                    ):  # Check connectivity
                        angle_raw_list.append([bond[0], bond[1], atom])

        for angle in angle_raw_list:
            # check uniqueness (avoid repetition of the same angle indexes flipped)
            if angle not in _angle_list and angle[::-1] not in _angle_list:
                _angle_list.append(angle)

        _angle_list.sort(
            key=lambda x: x[0]
        )  # sort the list by the first index

        return _angle_list.copy()

    @cached_property
    def dihedral_list(self) -> list:
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
    def bond_values(self) -> list:
        _bond_values = []

        for bond in self.bond_list:
            i, j = bond
            _bond_values.append(self.distance_matrix[i][j])

        return _bond_values.copy()

    @cached_property
    def angle_values(self) -> list:
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
    def dihedral_values(self) -> list:
        """Calculate the dihedral angles.

        Done with the atan(sin_phi, cos_phi) formula.

        """
        dih_val = []

        for dihedral in self.dihedral_list:
            a, b, c, d = dihedral
            t = np.cross(
                self.r_vector_matrix[a][b], self.r_vector_matrix[b][c]
            )
            u = np.cross(
                self.r_vector_matrix[b][c], self.r_vector_matrix[c][d]
            )
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

    def _read_xyzfile(self) -> None:
        try:
            with open(self.xyzfile, "r") as f:
                self._file_content_list = f.readlines()[2:]
            self._n_atoms = len(self._file_content_list)

        except FileNotFoundError as exc:
            raise FileNotFoundError(
                f'"{self.xyzfile}" geometry file was not found'
            ) from exc

    def _extract_atoms(self) -> None:
        self._atom_labels = []
        for atom in self._file_content_list:
            self._atom_labels.append(atom.strip().split()[0].capitalize())

    def _extract_positions(self) -> None:
        self._position_array = np.zeros([self._n_atoms, 3], dtype=float)

        for i, atom in enumerate(self._file_content_list):
            self._position_array[i] = np.array(
                atom.strip().split()[1:], dtype=float
            )

    def _extend_labels(self) -> None:
        counter_labels = list(set(self._atom_labels))
        counter = np.zeros(len(counter_labels))

        extended_labels = []

        for label in self._atom_labels:
            counter[counter_labels.index(label)] += 1
            extended_labels.append(
                label + str(int(counter[counter_labels.index(label)]))
            )

        self._atom_labels = extended_labels

MoleculeHandler = Jeremy

if __name__ == "__main__":
    # a = Jeremy("upgrades/xanthine_opt_iroot_1.xyz")
    # b = Jeremy("upgrades/xanthine_FC.xyz")
    pass


class Maxime(Jeremy):

    def extract_radius_atoms(
            self, 
            point: list[float]= [0,0,0],
            r:float=1,
            new_filename='None'
        ) -> list:

        if len(point) == 3:
            point == np.array(point)
        else:
            print(f'The inputed point is {len(point)} dimensional.')
            return 

        print(f'Selected point is {point}. Sphere radius is {r}.')

        preserved_atoms_indices = []
        frontier_atoms_indices = []

        for index, row in self.xyz_df.iterrows():
            atom_coords = np.array([row['x'], row['y'], row['z']])
            # print(f'Coordinates of atom{index} are {atom_coords}. Distance to selected point is: {np.linalg.norm(atom_coords - point)}. Included in sphere -> {np.linalg.norm(atom_coords - point)< r}')
            if np.linalg.norm(atom_coords - point)< r:
                preserved_atoms_indices.append(index)
                # print(f'Coordinates of atom {index} are {atom_coords}. Distance to selected point is: {np.linalg.norm(atom_coords - point)}, r= {r}.')
                
                if np.linalg.norm(atom_coords - point) < r and np.linalg.norm(atom_coords - point) > (r - 2):
                    frontier_atoms_indices.append(index)
                    # print(f'Coordinates of atom {index} are {atom_coords}. Distance to selected point is: {np.linalg.norm(atom_coords - point)}.')
        
        frontier_df = self.xyz_df.iloc[[i for i in frontier_atoms_indices]]

        completion_indices = self.complete_molecules(frontier_df)

        total_atom_indices = preserved_atoms_indices + completion_indices

        print('\n\n\n\nTOTAL ATOMS IN THE CASE IS')
        print(total_atom_indices)
        

        preserved_df = self.xyz_df.iloc[[i for i in total_atom_indices]]

        print(len(preserved_df))

        if new_filename != 'None':
            with open(new_filename, 'w') as f:
                f.write(f'{len(preserved_df)}\n')
                f.write(f'Extracted atoms from {point} with {r} radius\n')
                for index, row in preserved_df.iterrows():
                    f.write(f'{row["Atom"]} {row["x"]} {row["y"]} {row["z"]}\n')
        
        return preserved_df, r, point, frontier_atoms_indices
    
    def complete_molecules(self, frontier_df: pd.DataFrame) -> pd.DataFrame:

        frontier_indices = [i for i, row in frontier_df.iterrows()]

        total_molecule_indices = []
        
        while len(frontier_indices) >= 1: 
            print(len(frontier_indices))

            total_molecule_indices += self.complete_molecule(frontier_indices[0])
            frontier_indices_temp = [i for i in frontier_indices if i not in total_molecule_indices]
            frontier_indices = frontier_indices_temp
            print(len([i for i in frontier_indices if i not in total_molecule_indices]))
        
        print(total_molecule_indices)

        return total_molecule_indices

    def complete_molecule(self, index):

        starting_point = index
    
        already_checked = []
        already_checked.append(starting_point)
        molecule_indices = []
        neighbour_indices = [index for index, value in enumerate(self.connectivity_matrix[starting_point]) if value == 1]
        molecule_indices += neighbour_indices

        not_checked = [i for i in molecule_indices if i not in already_checked]

        while len(not_checked) >= 1: 
            starting_point = not_checked[0]
            already_checked.append(starting_point)
            neighbour_indices = [index for index, value in enumerate(self.connectivity_matrix[starting_point]) if value == 1 and index not in molecule_indices]
            molecule_indices += neighbour_indices
            not_checked = [i for i in molecule_indices if i not in already_checked]
        
        return(molecule_indices)



    # @cached_property
    # def strict_radius(self):
    #     return self.extract_radius_atoms()
    
    # def build_cluster(self):
    #     preserved_df, r, point = self.strict_radius

    #     for index, row in preserved_df.iterrows():
    #         atom_coords = np.array([row['x'], row['y'], row['z']])
    #         print(f'Coordinates of atom{index} are {atom_coords}. Distance to selected point is: {np.linalg.norm(atom_coords - point)}. Included in sphere -> {np.linalg.norm(atom_coords - point)< r}')
    #         if np.linalg.norm(atom_coords - point)< r:
    #             preserved_atoms_indices.append(index)



    
    # def write_to_xyz(self, new_filename:str='None') -> None:


        

if __name__ == "__main__":
    
    # display help message if there is no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    m = Maxime(args.f)
    try:
        refp = [float(i) for i in args.p.strip().split()]
    except:
        exit()
    # new_filename=args.o
    m.extract_radius_atoms(point=refp, r=args.r, new_filename=args.o)

