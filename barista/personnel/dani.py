#!/usr/bin/env python3
import argparse
import sys
import numpy as np
from typing import List, Tuple, Optional

np.set_printoptions(linewidth=np.inf)

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Extract the final geometry of an optimization run in Gaussian. Extract vibrational frequencies if available""",
    epilog="""It should run appropriately with gaussian 16.""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="Input file name"
)

parser.add_argument("-o", type=str, default=None, help="Output file name. Default is Input_filename.xyz")

parser.add_argument("--frequencies", action='store_true', help="Extract frequencies to Input_filename_freq.dat")


class Dani:
    """
    A class for extracting and processing data from Gaussian output files.
    
    This class extracts final geometry coordinates and vibrational frequencies 
    from Gaussian output files, and can generate XYZ format files.
    
    Attributes
    ----------
    gaussian_output_filename : str
        Path to the Gaussian output file
    atomic_dict : Dict[int, str]
        Dictionary mapping atomic numbers to element symbols
    """

    # Mapping of atomic numbers to element symbols
    atomic_dict = {
        1: 'H',
        6: 'C',
        7: 'N',
        8: 'O',
    }

    # =========================================================================
    #     Special Methods
    # =========================================================================
    
    def __init__(self, gaussian_output_filename: str) -> None:
        """
        Initialize the Dani class with a Gaussian output file.
        
        Parameters
        ----------
        gaussian_output_filename : str
            Path to the Gaussian output file
        """
        self.gaussian_output_filename = gaussian_output_filename
        self._file_content: Optional[List[str]] = None
        self.starting_index: Optional[int] = None
        self.end_index: Optional[int] = None

        self._initialize()

    # =========================================================================
    #     Class Methods
    # =========================================================================
    
    @classmethod
    def from_file(cls, gaussian_output_filename: str) -> 'Dani':
        """
        Factory method to create and initialize a Dani instance from a file.
        
        Parameters
        ----------
        gaussian_output_filename : str
            Path to the Gaussian output file
            
        Returns
        -------
        Dani
            An initialized Dani instance
        """
        instance = cls(gaussian_output_filename)
        return instance

    # =========================================================================
    #     Public Methods
    # =========================================================================
    
    def generate_xyzfile(self, filename: Optional[str] = None) -> None:
        """
        Generate an XYZ format file from the extracted coordinates.
        
        Parameters
        ----------
        filename : str, optional
            Output filename, if None the input filename with .xyz extension will be used
            
        Returns
        -------
        None
        """
        if filename is None:
            filename = (
                self.gaussian_output_filename.replace(".log", "").replace(".out", "") + ".xyz"
            )

        with open(filename, "w") as out_file:
            out_file.write(str(len(self._file_content[self.starting_index:self.end_index])) + "\n")
            out_file.write('Imagine a world where atoms were very small kittens\n')

            for atom in self._file_content[self.starting_index:self.end_index]:
                list_to_print = [self.atomic_dict[int(atom.strip().split()[1])]] + atom.strip().split()[-3:] + ['\n']
                out_file.write(' '.join(list_to_print))

    def extract_vibrational_modes(self, filename: Optional[str] = None) -> None:
        """
        Extract vibrational frequencies and IR intensities from the Gaussian output.
        
        Parameters
        ----------
        filename : str, optional
            Output filename, if None the input filename with _frequencies.dat extension will be used
            
        Returns
        -------
        None
        """
        intensities: List[str] = []
        freqs: List[str] = []
        
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
                out_file.write(' '.join([str(mode+1), str(frequency), str(intensities[mode]), '\n']))

    # =========================================================================
    #     Private Methods
    # =========================================================================
    
    def _initialize(self) -> None:
        """
        Initialize the object by reading the file and extracting coordinates.
        
        Returns
        -------
        None
        """
        self._get_file_content()
        self._get_coordinates_section()

    def _get_file_content(self) -> None:
        """
        Read the content of the Gaussian output file.
        
        Returns
        -------
        None
        
        Raises
        ------
        FileNotFoundError
            If the specified file cannot be found
        """
        try:
            with open(self.gaussian_output_filename, "r") as f:
                self._file_content = f.readlines()

        except FileNotFoundError:
            raise FileNotFoundError(f"Could not find file: {self.gaussian_output_filename}")

    def _get_coordinates_section(self) -> Tuple[int, int]:
        """
        Extract the section containing the final coordinates from the Gaussian output file.
        
        This method searches backwards through the file to find the last occurrence of
        the coordinates section, which should contain the final optimized geometry.
        
        Returns
        -------
        tuple
            A tuple containing the start and end indices of the coordinates section
        """
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


if __name__ == "__main__":
    # display help message if there is no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    a = Dani.from_file(args.f)
    # print('coordinates', a._get_coordinates_section())
    # i = a._get_coordinates_section()[0]
    # f = a._get_coordinates_section()[1]
    # print(i, f)
    # for atom in a._file_content[i:f]:
    #     print([a.atomic_dict[int(atom.strip().split()[1])]] + atom.strip().split()[-3:])
    
    a.generate_xyzfile(args.o)
    if args.frequencies:
        a.extract_vibrational_modes()
