#!/usr/bin/env python3

from functools import cached_property

import argparse
import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
import sys

# Parser is used to input via terminal the required arguments
parser = argparse.ArgumentParser(
    description="""Extract energies of an ORCA NEB calculation.""",
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="ORCA NEB optmization file."
)
parser.add_argument("-o", type=str, default=None, help="Output file name.")

parser.add_argument(
    "-en",
    type=float,
    default=0,
    help="Reference energy of the ground state optimized energy",
)

parser.add_argument("-u", type=str, default="eV", help="Energy units")

class Stef:
    """
    A class for extracting and processing energies from ORCA NEB calculations.
    
    This class reads an ORCA NEB output file, extracts image energies along the
    reaction path, and provides functionality to save and potentially plot these energies.
    
    Parameters
    ----------
    filename : str
        Path to the ORCA NEB output file
    reference_energy : float, optional
        Reference energy to subtract from all extracted energies, default is 0.0
    units : str, optional
        Energy units to use, default is "eV"
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

        self._initialize()

    def _initialize(self):
        """
        Initialize the object by reading the file content.
        
        This method is called automatically during object initialization.
        
        Returns
        -------
        None
        """
        self._get_file_content()

    def _get_file_content(self):
        """
        Open the file and store its content as a list.
        
        Raises
        ------
        ValueError
            If there is no filename specified.
        FileNotFoundError
            If the file does not exist.

        Returns
        -------
        None
        """
        if self.filename is None:
            raise ValueError("No result file was specified")

        try:
            with open(self.filename, "r") as file:
                self.content_list = file.readlines()
        except FileNotFoundError:
            raise FileNotFoundError(f"File {self.filename:s} does not exist.")


    def _get_image_energies(self) ->  None:
        """
        Extract image energies from the ORCA NEB output file.
        
        This method searches for the PATH SUMMARY section in the output file
        and extracts the image energies from this section. The results are
        stored in the self._image_energies attribute.
        
        Returns
        -------
        None
        """
        start_section = 0 
        end_section = 0 

        for index, line in enumerate(self.content_list):
            if "PATH SUMMARY" in line:
                start_section = index + 5 
        for index, line in enumerate(self.content_list[start_section:]):
            if "Straight line distance between images along the path" in line:
                end_section = start_section + index - 2 

        if start_section == end_section:
            print('No NEB path was found in File.')
            return None

        root_energies = []

        for line in self.content_list[start_section:end_section]:
            root_energies.append(float(line.strip().split()[2]))

        self._image_energies = np.array(root_energies)

        return None
        
    @cached_property
    def image_energies(self) -> NDArray[np.float64]:
        """
        Get the energies of all images along the NEB path.
        
        This property calculates the energies relative to the reference energy
        and converts them to the specified units if necessary. The conversion
        factor 27.2114 is used to convert from Hartree to eV.
        
        Returns
        -------
        NDArray[np.float64]
            Array of energy values for each image
        """
        self._get_image_energies()

        if self.units == 'eV':
            return 27.2114*(np.copy(self._image_energies) - self.reference_energy)

        return np.copy(self._image_energies) - self.reference_energy

    def save_energies(self, filename:str) -> None:
        """
        Save the image energies to a file.
        
        Parameters
        ----------
        filename : str
            Path to the output file
            
        Returns
        -------
        None
        """
        print(self.image_energies)
        with open(filename, 'w') as f:
            for image, energy in enumerate(self.image_energies):
                f.write(f'{image} {energy:8.5f}\n')
                
    def savefig(self, filename):
        filename = filename.replace('.dat', '.png')
        plt.plot([i for i in range(len(self.image_energies))], self.image_energies*27.2114)
        plt.ylabel('$\Delta$E w.r.t. $S_0$ FC')
        plt.savefig(filename, dpi=600)

    def plot(self):
        plt.plot([i for i in range(len(self.image_energies))], self.image_energies*27.2114)
        plt.ylabel('$\Delta$E w.r.t. $S_0$ FC')
        plt.show()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    a = Stef(args.f, args.en, args.u)

    if args.o is not None:
        a.save_energies(args.o)
        a.savefig(args.o)
    a.plot()

#     a = Stef("tests/stef/neb.dat", 0, "eV")
# 
#     print(a.image_energies)
#     a.save_energies('tests/stef_out.dat')
