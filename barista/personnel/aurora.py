#!/usr/bin/env python3
from dataclasses import dataclass
import argparse
import sys
import pandas as pd
import numpy as np
from numpy.typing import NDArray
from typing import Optional, List, Tuple
import matplotlib.pyplot as plt

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Exctract the potential energy surface path and energies of of a PES branch.
    Takes as input the ground state equilibrium geometry xyz file  and the xyz file of the last point of the PES branch.
    
    \n\nThis script determines the path, highlights, and generates a pandas dataframe that can be later plotted and saved as an image.""",
    epilog=""" It is required that the last branch file names follow the brewer naming convention prefix_initial-state_final-state_identifier.xyz (see documentation). 
    It is also required that the formatting of the xyz file is the same as the one generated by laura (see documentation)""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="FC xyz file with energy as comment."
)

parser.add_argument(
    "-b", type=str, required=True, default=None, help="Last point in the branch."
)

parser.add_argument("-o", type=str, default=None, help="Output image name.")


@dataclass
class AuroraConfiguration:
    """
    Configuration class for Aurora with default parameters.
    
    Parameters
    ----------
    fc_file : str, optional
        Path to the Franck-Condon xyz file with energy as comment, default is None
    last_branch_file : str, optional
        Path to the xyz file of the last point in the branch, default is None
    """
    fc_file: str = None
    last_branch_file: str = None


class Aurora:
    """
    A class for extracting and analyzing potential energy surface (PES) paths.
    
    This class processes xyz files containing energy data to extract PES branch information,
    calculate relative energies, and visualize energy paths along different electronic states.
    
    Parameters
    ----------
    fc_filename : str
        Path to the Franck-Condon xyz file with energy as comment
    last_branch_filename : str
        Path to the xyz file of the last point in the branch
    """
    # =========================================================================
    #     Special Methods
    # =========================================================================

    def __init__(self, fc_filename: str, last_branch_filename: str) -> None:
        self.fc_filename = fc_filename
        self.last_branch_filename = last_branch_filename
        self._pre_identifier_name = None
        self.reference_energy = None
        self._step_identifiers = None
        self._filenames = None
        self._highlighted_path = None
        self._pes_dataFrame = pd.DataFrame({
            'Step': [],
            'State': [],
            'Relative Energy': [],
        })
        self._initialize()

    # =========================================================================
    #     Class Methods
    # =========================================================================

    @classmethod
    def from_files(cls, fc_file: str, last_branch_file: str) -> 'Aurora':
        """
        Factory method to create and initialize an Aurora instance.
        
        Parameters
        ----------
        fc_file : str
            Path to the Franck-Condon xyz file
        last_branch_file : str
            Path to the xyz file of the last point in the branch
            
        Returns
        -------
        Aurora
            An initialized Aurora instance
        """
        instance = cls(fc_file, last_branch_file)
        return instance
    
    # =========================================================================
    #     Public Methods
    # =========================================================================

    def get_relative_energies(self, xyz_file: str) -> NDArray[np.float64]:
        """
        Obtain energy values of an xyz file that has energy values as comment string.
        
        This method reads the xyz file, extracts the energy values from the comment line,
        and calculates relative energies with respect to the reference energy in eV.
        
        Parameters
        ----------
        xyz_file : str
            Path to the xyz file containing energy data
            
        Returns
        -------
        numpy.ndarray
            Array of relative energies in eV
            
        Raises
        ------
        ValueError
            If the reference energy has not been determined
        """
        if self.reference_energy is None:
            raise ValueError('Reference energy has not been determined. Call get_reference_energy first.')

        try:
            with open(xyz_file, 'r') as f:
                cont = f.readlines()
            relative_energies = (np.array([float(energy) for energy in cont[1].strip().split()]) - self.reference_energy) * 27.211

            return relative_energies
        except FileNotFoundError:
            print('The file "%s" does not exist.' % xyz_file)

            return None

    def pes_DataFrame(self) -> pd.DataFrame:
        """
        Get a copy of the PES dataframe.
        
        Returns
        -------
        pandas.DataFrame
            A copy of the dataframe containing all PES data
        """
        return self._pes_dataFrame.copy()

    def plot_branch(self) -> None:
        """
        Plot the PES branch and display the figure.
        
        This method calls _preplot_branch to prepare the plot and then displays it.
        
        Returns
        -------
        None
        """
        self._preplot_branch()
        plt.show()

    def save_branch_plot(self, name: Optional[str] = None) -> None:
        """
        Save the PES branch plot to a file.
        
        Parameters
        ----------
        name : str, optional
            The output filename, if None the last branch filename with jpg extension will be used
            
        Returns
        -------
        None
        """
        self._preplot_branch()
        if name is None:
            image_filename = self.last_branch_filename.replace('xyz', 'jpg')
        else:
            image_filename = str(name)
        plt.savefig(image_filename, dpi=300)

    # =========================================================================
    #     Private Methods
    # =========================================================================
    
    def _initialize(self) -> None:
        """
        Initialize all required data before calculations.
        
        This method calls several private methods to set up the reference energy,
        get name and identifiers, determine filenames, and fill the dataframe.
        
        Returns
        -------
        None
        """
        self._set_reference_energy()
        self._get_name_and_identifiers()
        self._get_filenames()
        self._fill_dataframe()

    def _get_name_and_identifiers(self) -> Tuple[str, NDArray[np.int_]]:
        """
        Extract the prefix name and identifiers from the last branch filename.
        
        This method parses the last branch filename to extract the pre-identifier name
        and the step identifiers that define the path through the PES.
        
        Returns
        -------
        tuple
            Tuple containing the pre-identifier name and step identifiers
        """
        split_name = self.last_branch_filename.replace('.xyz', '')

        for i, item in enumerate(split_name):
            try:
                int(item)
                index = i
                break
            except ValueError:
                pass

        path_information = split_name[index:].replace('ci', '-1_-1_-1').split('_')
        self._pre_identifier_name = split_name[:index - 1]

        path_information = [int(i) for i in path_information if i != '']

        step_identifiers = np.array(path_information).reshape([-1, 3])

        self._step_identifiers = step_identifiers

        return self._pre_identifier_name, self._step_identifiers

    def _get_filenames(self) -> List[str]:
        """
        Generate filenames for all points in the branch.
        
        This method creates a list of filenames for all points in the branch
        based on the pre-identifier name and step identifiers.
        
        Returns
        -------
        list
            List of filenames for all points in the branch
            
        Raises
        ------
        ValueError
            If the general name or path step identifiers have not been determined
        """
        if self._pre_identifier_name is None:
            raise ValueError('The general name was not identified. Call get_name_and_identifiers.')

        if self._step_identifiers is None:
            raise ValueError('Path step identifiers were not determined. Call get_name_and identifiers.')

        a, b, c = self._step_identifiers[0]

        first_filename = ''.join(self._pre_identifier_name) + '_%i_%i_%i.xyz' % (a, b, c)

        filenames = [first_filename]

        for identifier in self._step_identifiers[1:]:
            if max(identifier) == -1:
                filenames.append(filenames[-1].replace('.xyz', '') + '_ci.xyz')
            else:
                a, b, c = identifier
                filenames.append(filenames[-1].replace('.xyz', '') + '_%i_%i_%i.xyz' % (a, b, c))

        self._filenames = filenames

        return self._filenames

    def _set_reference_energy(self) -> float:
        """
        Set the reference energy from the Franck-Condon state.
        
        This method reads the FC file, extracts the energies, sets the reference energy
        to the ground state energy, and initializes the PES dataframe with FC energies.
        
        Returns
        -------
        float
            The reference energy
        """
        with open(self.fc_filename, 'r') as f:
            cont = f.readlines()

        fc_energies = (np.array([float(energy) for energy in cont[1].strip().split()]))

        self.reference_energy = fc_energies[0]

        fc_energies -= self.reference_energy
        fc_energies *= 27.2114
        step_dataframe = pd.DataFrame({
            'State': [i for i in range(len(fc_energies))],
            'Step': np.zeros(len(fc_energies)),
            'Relative Energy': fc_energies
        })

        self._pes_dataFrame = pd.concat([self._pes_dataFrame, step_dataframe])

        return self.reference_energy

    def _fill_dataframe(self) -> None:
        """
        Fill the dataframe with all the energies and steps of the branch.
        
        This method iterates through all filenames, extracts the energies,
        and populates the PES dataframe with the relative energies for each state and step.
        
        Returns
        -------
        None
        """
        for index, filename in enumerate(self._filenames):
            energies = self.get_relative_energies(filename)
            step_array = np.zeros(len(energies)) + index + 1
            step_dataframe = pd.DataFrame({
                'State': [i for i in range(len(energies))],
                'Step': step_array,
                'Relative Energy': energies
            })
            self._pes_dataFrame = pd.concat([self._pes_dataFrame, step_dataframe])

    def _preplot_branch(self) -> plt.Figure:
        """
        Prepare the plot for the PES branch.
        
        This method sets up the plot with appropriate labels, plots all states,
        and highlights the path through the PES.
        
        Returns
        -------
        matplotlib.pyplot
            The prepared plot
        """
        unique_states = self._pes_dataFrame['State'].unique()

        plt.title(self.last_branch_filename)

        x, x_tics = self._get_plot_labels()

        plt.xticks(x, x_tics)
        plt.ylabel('Relative energy / eV')

        for state in unique_states:
            reduced_df = self._pes_dataFrame[self._pes_dataFrame['State'] == state]
            plt.plot(reduced_df['Step'], reduced_df['Relative Energy'], alpha=.5, linestyle=':')
            plt.scatter(reduced_df['Step'], reduced_df['Relative Energy'], alpha=.5)

        plt.plot(*self._get_path_energies(), c='black', linestyle='-.', label='PES path')

        plt.legend()

        return plt

    def _get_path_energies(self) -> Tuple[List[int], List[np.ndarray]]:
        """
        Get the energies along the highlighted path.
        
        This method extracts the energy values for each step along the highlighted path
        through the PES.
        
        Returns
        -------
        tuple
            Tuple containing x-coordinates and corresponding energy values
        """
        highlighted_values = []
        x = []

        for step, actual_state in enumerate(self._get_highlighted_path()):
            highlighted_values.append(
                self._pes_dataFrame.query(
                    "State == %s.0 and Step == %s.0" % (actual_state, step)
                )['Relative Energy'].values)

            x.append(step)

        return x, highlighted_values

    def _get_plot_labels(self) -> Tuple[List[int], List[str]]:
        """
        Generate x-axis labels for the plot.
        
        This method creates appropriate labels for each step in the PES path.
        
        Returns
        -------
        tuple
            Tuple containing x-coordinates and corresponding tick labels
        """
        x_tics = ['FC']

        for i in self._step_identifiers:
            if max(i) == -1:
                x_tics.append('CI')
            else:
                x_tics.append('Minimum\n #%i $S_{%i}$' % (i[-1], i[-2]))

        x = [i for i in range(1 + int(max(self._pes_dataFrame['Step'].unique())))]

        return x, x_tics

    def _get_highlighted_path(self) -> List[int]:
        """
        Get the highlighted path through the PES.
        
        This method determines the sequence of electronic states that form
        the path through the PES.
        
        Returns
        -------
        list
            List of state indices forming the highlighted path
        """
        if self._highlighted_path is None:
            highlighted = [self._step_identifiers[0][0]]

            for identifier in self._step_identifiers:
                highlighted.append(identifier[1])

            for i, highlight in enumerate(highlighted):
                if highlight == -1:
                    highlighted[i] = highlighted[i - 1]

            self._highlighted_path = highlighted

        return self._highlighted_path.copy()


if __name__ == "__main__":
    # display help message if there is no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    a = Aurora.from_files(args.f, args.b)

    if args.o is not None:
        a.save_branch_plot(args.o)
    else:
        a.plot_branch()