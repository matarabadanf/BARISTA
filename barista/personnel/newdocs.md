# Python Code Documentation

This document provides comprehensive documentation for Python scripts.

## Table of Contents

### Command Line Interfaces

* [alberto](#alberto)
* [aurora](#aurora)
* [dani](#dani)
* [emma](#emma)
* [fonsi](#fonsi)
* [javi](#javi)
* [laura](#laura)
* [lorena](#lorena)
* [maxime](#maxime)
* [rodrigo](#rodrigo)
* [stef](#stef)

### Classes

* [alberto.Alberto](#alberto-alberto)
* [aurora.Aurora](#aurora-aurora)
* [aurora.AuroraConfiguration](#aurora-auroraconfiguration)
* [dani.Dani](#dani-dani)
* [emma.Emma](#emma-emma)
* [fonsi.Fonsi](#fonsi-fonsi)
* [javi.Javi](#javi-javi)
* [jeremy.Jeremy](#jeremy-jeremy)
* [jeremy.MoleculeHandler](#jeremy-moleculehandler)
* [laura.Laura](#laura-laura)
* [maxime.Jeremy](#maxime-jeremy)
* [maxime.Maxime](#maxime-maxime)
* [NXSpecReader.NXSpecReader](#nxspecreader-nxspecreader)
* [rodrigo.Rodrigo](#rodrigo-rodrigo)
* [stef.NebExtractor](#stef-nebextractor)

---

# Command Line Interfaces

<!-- Source: alberto.py -->
## CLI Arguments: alberto

Takes a excited state optimization in ORCA and plots the energy of each state at each step of the opitmizaiton and the actual root.

### Usage
```
docs3.py [-h] -f F [-o O] [-en EN] [-u U]
```

### Arguments

#### Options

* `-f` (str) [required]: ORCA optimization file
  * Default: `42`

* `-o` (str): Output image filename
  * Default: `True`

* `-en` (float): Reference energy of the ground state optimized energy
  * Default: `0`

* `-u` (str): Energy units
  * Default: `eV`


---

<!-- Source: aurora.py -->
## CLI Arguments: aurora

Exctract the potential energy surface path and energies of of a PES branch.
    Takes as input the ground state equilibrium geometry xyz file  and the xyz file of the last point of the PES branch.
    
    

This script determines the path, highlights, and generates a pandas dataframe that can be later plotted and saved as an image.

### Usage
```
docs3.py [-h] -f F -b B [-o O]
```

### Arguments

#### Options

* `-f` (str) [required]: FC xyz file with energy as comment.
  * Default: `42`

* `-b` (str) [required]: Last point in the branch.

* `-o` (str): Output image name.


---

<!-- Source: dani.py -->
## CLI Arguments: dani

Extract the final geometry of an optimization run in Gaussian. Extract vibrational frequencies if available

### Usage
```
docs3.py [-h] -f F [-o O] [--frequencies]
```

### Arguments

#### Options

* `-f` (str) [required]: Input file name
  * Default: `42`

* `-o` (str): Output file name. Default is Input_filename.xyz

* `--frequencies`: Extract frequencies to Input_filename_freq.dat
  * Default: `False`


---

<!-- Source: emma.py -->
## CLI Arguments: emma

Takes the ground state geometry of a molecule and compares it to optimized TDDFT geometries run in ORCA.        
It compares the energy, RMSD and final root, plotting the results in an image. 

### Usage
```
docs3.py [-h] -fc FC -ex EX [EX ...] [-O O] [-o O] [-r R] [-md MD] [-i I]
```

### Arguments

#### Options

* `-fc` (str) [required]: Optimized ground state .xyz file
  * Default: `42`

* `-ex` (str) [required]: Excited state optimization .out "filenames" (it is important that the name ends with *_N.in.out)
  * Default: `42`

* `-O` (str): Output file extension (default=".in.out")
  * Default: `42`

* `-o` (str): Output image filename
  * Default: `42`

* `-r` (str): Generate a report table in output file
  * Default: `42`

* `-md` (str): Generate report in MD format (default is True)
  * Default: `True`

* `-i` (bool): Interactive plot mode (default is False)
  * Default: `False`


---

<!-- Source: fonsi.py -->
## CLI Arguments: fonsi

Takes a brewer optimizaiton energies.dat and plots the energy at each step.

### Usage
```
docs3.py [-h] -f F [-o O] [-en EN] [-u U] [--interactive]
```

### Arguments

#### Options

* `-f` (str) [required]: energies.dat file
  * Default: `42`

* `-o` (str): Output image filename

* `-en` (float): Reference energy of the ground state optimized energy
  * Default: `0`

* `-u` (str): Energy units
  * Default: `eV`

* `--interactive`: Open interactive plot.
  * Default: `False`


---

<!-- Source: javi.py -->
## CLI Arguments: javi

Takes the lower and upper state gradients in a CI with the NAC vector and characterizes the type oc CI.        
Results can be plotted interactively using plotly.

### Usage
```
docs3.py [-h] -g0 G0 -g1 G1 -nac NAC [--interactive]
```

### Arguments

#### Options

* `-g0` (str) [required]: Lower state gradient file.

* `-g1` (str) [required]: Higher state gradient file.

* `-nac` (str) [required]: NAC vector file.

* `--interactive`: Open interactive plot.
  * Default: `False`


---

<!-- Source: laura.py -->
## CLI Arguments: laura

Exctract the final geometry of an optimization run in ORCA.

### Usage
```
docs3.py [-h] -f F [-o O] [--no_energies]
```

### Arguments

#### Options

* `-f` (str) [required]: Input file name
  * Default: `42`

* `-o` (str): Output file name

* `--no_energies`: Dont include energies in the xyzfile as comment.
  * Default: `True`


---

<!-- Source: lorena.py -->
## CLI Arguments: lorena

Takes the an xyz file and shakes the geometry

### Usage
```
docs3.py [-h] [-f F] [-d D] [-o O]
```

### Arguments

#### Options

* `-f` (str): xyz file
  * Default: `42`

* `-d` (float): Maximum displacement
  * Default: `0.1`

* `-o` (str): Output file
  * Default: `42`


---

<!-- Source: maxime.py -->
## CLI Arguments: maxime

Extract atoms from an xyz file with a center point and a radius.

### Usage
```
docs3.py [-h] -f F [-o O] [-r R] [-p P]
```

### Arguments

#### Options

* `-f` (str) [required]: Original XYZ file.
  * Default: `42`

* `-o` (str): Output xyz filename.
  * Default: `None`

* `-r` (float): Sphere radius.
  * Default: `1.0`

* `-p` (str): Central sphere point.
  * Default: `0 0 0`


---

<!-- Source: rodrigo.py -->
## CLI Arguments: rodrigo

Generates the TDDFT absorption spectre from an ORCA TDDFT calculation.         It generates spectre .jpg image.
         It can generate a peaks.dat and a spectra.dat that contains the peaks and spectre.

### Usage
```
docs3.py [-h] -F F [-o O] [-u U] [--gaussian] [-gauss_disp GAUSS_DISP] [-exp EXP] [-exp_units EXP_UNITS] [-exp_shift EXP_SHIFT] [-semi SEMI] [--interactive]
```

### Arguments

#### Options

* `-F` (str) [required]: Output ORCA file
  * Default: `42`

* `-o` (str): Output image filename. Default is output_file_spectre.jpg

* `-u` (str): Units, nm or eV. Default is eV.
  * Default: `eV`

* `--gaussian`: Include gaussian wrap for peaks.
  * Default: `False`

* `-gauss_disp` (float): Gaussian dispersion in the current energy unit.

* `-exp` (str): Experimental data for comparison.

* `-exp_units` (str): Experimental data energy units. Default is eV.
  * Default: `eV`

* `-exp_shift` (float): Experimental spectrum shift in current energy units. Default is 0.
  * Default: `0.0`

* `-semi` (str): Semiclassical nx log file.

* `--interactive`: Interactive mode, will not save a .jpg file.
  * Default: `False`


---

<!-- Source: stef.py -->
## CLI Arguments: stef

Takes the ground state geometry of a molecule and compares it to optimized TDDFT geometries run in ORCA.        
It compares the energy, RMSD and final root, plotting the results in an image. 

### Usage
```
docs3.py [-h] -f F [-en EN] [-u U] [-o O] [-dat DAT] [--landscape] [--spline] [--barrier]
```

### Arguments

#### Options

* `-f` (str) [required]: Trajectory .xyz file
  * Default: `42`

* `-en` (float): Total absolute energy in Hartree at Franc-Condon
  * Default: `0`

* `-u` (str): Energy units for the plot (default="eV")
  * Default: `eV`

* `-o` (str): Output image name.
  * Default: `default.jpg`

* `-dat` (str): Save data to .dat file.
  * Default: `default.dat`

* `--landscape`: Use landscape layout
  * Default: `False`

* `--spline`: Use spline
  * Default: `False`

* `--barrier`: Plot absolute barrier
  * Default: `False`


---

# Classes

<!-- Source: alberto.py -->
### Alberto

*Full path: `alberto.Alberto`*

Takes an excited state optimization in ORCA and plots the energy of each state at each
step of the optimization and the actual root.

This class parses ORCA output files from TD-DFT/TDA calculations to track and
visualize the energy evolution of electronic states during geometry optimization.

Attributes
----------
filename : str
    Name of the ORCA output file.
reference_energy : float
    Value of the reference energy (typically ground state energy) in Hartree.
units : str
    Units for the energy values in plots. Default is "eV".
content_list : List[str]
    Content of the ORCA output file as a list of lines.
state_indexes : List[int]
    List of state indexes found in the file.
number_of_states : int
    Total number of excited states found.
number_of_steps : int
    Number of optimization steps.
cis_array : np.ndarray
    Array containing SCF energies at each step.
energy_array : np.ndarray
    2D array containing energy values for each state at each step.
bulk_energy_array : np.ndarray
    Copy of energy_array before reference shift.
curr_energy : np.ndarray
    Energy of the current root at each optimization step.
_starting_root : int
    Initial root used for the optimization.
_final_root : int
    Final root after optimization.

Methods
-------
_initialize() -> None
    Initialize all data structures and parse the output file.
_get_file_content() -> None
    Open and read the output file.
_get_state_indexes() -> None
    Extract state indexes from the output file.
_get_cis_energies() -> None
    Extract SCF energies from the output file.
_get_abs_energy_array() -> None
    Calculate absolute energies for all states.
_get_relative_energy() -> None
    Calculate energies relative to the reference energy.
_get_current_energy() -> None
    Extract the current energy surface at each step.
set_units(units: str) -> None
    Change the energy units used for plotting.
_rescale_units(prev_units: str) -> None
    Rescale energy values according to the selected units.
set_imagename(image_name: str) -> None
    Set the output image name.
_prepare_plot() -> Tuple[plt.Figure, plt.Axes]
    Prepare the energy plot.
plot() -> None
    Display the energy plot.
generate_image(output_image_name: str) -> None
    Save the energy plot to an image file.

### Properties

#### `final_root`

Get the final root after optimization.

Returns
-------
int
    Index of the final root.

#### `starting_root`

Get the initial root used for the optimization.

Returns
-------
int
    Index of the starting root.

### Methods

#### generate_image

```python
def generate_image(self, output_image_name: str) -> None:
```

Save the energy plot to an image file.

Parameters
----------
output_image_name : str
    Name of the output image file.

#### plot

```python
def plot(self) -> None:
```

Display the energy plot showing all states and the current surface.

Creates and displays a plot showing the energy evolution of all states
during the optimization, highlighting the current surface being followed.

#### set_imagename

```python
def set_imagename(self, image_name: str) -> None:
```

Sets the output image name for saving plots.

Parameters
----------
image_name : str
    Image name to set.

#### set_units

```python
def set_units(self, units: str) -> None:
```

Change the energy units used for plotting.

Parameters
----------
units : str
    Energy units to use ("eV" or "Hartree").
    
Raises
------
ValueError
    If an unsupported unit is specified.

---

<!-- Source: aurora.py -->
### Aurora

*Full path: `aurora.Aurora`*

A class for extracting and analyzing potential energy surface (PES) paths.

This class processes xyz files containing energy data to extract PES branch information,
calculate relative energies, and visualize energy paths along different electronic states.

Parameters
----------
fc_filename : str
    Path to the Franck-Condon xyz file with energy as comment
last_branch_filename : str
    Path to the xyz file of the last point in the branch

### Methods

#### get_relative_energies

```python
def get_relative_energies(self, xyz_file: str) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]:
```

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

#### pes_DataFrame

```python
def pes_DataFrame(self) -> pandas.core.frame.DataFrame:
```

Get a copy of the PES dataframe.

Returns
-------
pandas.DataFrame
    A copy of the dataframe containing all PES data

#### plot_branch

```python
def plot_branch(self) -> None:
```

Plot the PES branch and display the figure.

This method calls _preplot_branch to prepare the plot and then displays it.

Returns
-------
None

#### save_branch_plot

```python
def save_branch_plot(self, name: Optional[str] = None) -> None:
```

Save the PES branch plot to a file.

Parameters
----------
name : str, optional
    The output filename, if None the last branch filename with jpg extension will be used
    
Returns
-------
None

---

<!-- Source: aurora.py -->
### AuroraConfiguration

*Full path: `aurora.AuroraConfiguration`*

Configuration class for Aurora with default parameters.

Parameters
----------
fc_file : str, optional
    Path to the Franck-Condon xyz file with energy as comment, default is None
last_branch_file : str, optional
    Path to the xyz file of the last point in the branch, default is None

---

<!-- Source: dani.py -->
### Dani

*Full path: `dani.Dani`*

A class for extracting and processing data from Gaussian output files.

This class extracts final geometry coordinates and vibrational frequencies 
from Gaussian output files, and can generate XYZ format files.

Attributes
----------
gaussian_output_filename : str
    Path to the Gaussian output file
atomic_dict : Dict[int, str]
    Dictionary mapping atomic numbers to element symbols

### Methods

#### extract_vibrational_modes

```python
def extract_vibrational_modes(self, filename: Optional[str] = None) -> None:
```

Extract vibrational frequencies and IR intensities from the Gaussian output.

Parameters
----------
filename : str, optional
    Output filename, if None the input filename with _frequencies.dat extension will be used
    
Returns
-------
None

#### generate_xyzfile

```python
def generate_xyzfile(self, filename: Optional[str] = None) -> None:
```

Generate an XYZ format file from the extracted coordinates.

Parameters
----------
filename : str, optional
    Output filename, if None the input filename with .xyz extension will be used
    
Returns
-------
None

---

<!-- Source: emma.py -->
### Emma

*Full path: `emma.Emma`*

_No class documentation provided_

### Methods

#### internal_deviation

```python
def internal_deviation(self):
```

_No documentation provided_

#### rmsd_report_to_md

```python
def rmsd_report_to_md(self, md_file: str = ''):
```

_No documentation provided_

---

<!-- Source: fonsi.py -->
### Fonsi

*Full path: `fonsi.Fonsi`*

_No class documentation provided_

### Methods

#### plot

```python
def plot(self):
```

_No documentation provided_

#### save_fig

```python
def save_fig(self, name: str = ''):
```

_No documentation provided_

---

<!-- Source: javi.py -->
### Javi

*Full path: `javi.Javi`*

_No class documentation provided_

### Properties

#### `beta`

_No documentation provided_

#### `ci_energy`

_No documentation provided_

#### `is_rotation_needed`

_No documentation provided_

### Methods

#### E_A

```python
def E_A(self, x: float, y: float) -> float:
```

_No documentation provided_

#### E_B

```python
def E_B(self, x: float, y: float) -> float:
```

_No documentation provided_

#### average_energy

```python
def average_energy(self, x: float, y: float) -> float:
```

_No documentation provided_

#### energy_difference

```python
def energy_difference(self, x: float, y: float) -> float:
```

_No documentation provided_

#### generate_force_file

```python
def generate_force_file(self, xyz_file: str):
```

_No documentation provided_

#### plot_CI

```python
def plot_CI(self, max_grid: float = 1):
```

_No documentation provided_

---

<!-- Source: jeremy.py -->
### Jeremy

*Full path: `jeremy.Jeremy`*

_No class documentation provided_

### Properties

#### `atom_labels`

_No documentation provided_

#### `internal_list`

Return copy of the internal coordinate list.

Returns
-------
list

#### `internal_values`

Return copy of the list containing internal coordinate values.

Returns
-------
list

#### `position_array`

Return a copy of the xyz cooridnates array.

Returns
-------
ArrayLike

### Methods

#### clear_cache

```python
def clear_cache(self) -> None:
```

Clear cached properties.

Returns
-------
None

#### compare_internals

```python
def compare_internals(self, reference: 'Jeremy') -> numpy.ndarray:
```

_No documentation provided_

#### override_connectivity_matrix

```python
def override_connectivity_matrix(self, cm: Union[numpy._typing._array_like._SupportsArray[numpy.dtype[Any]], numpy._typing._nested_sequence._NestedSequence[numpy._typing._array_like._SupportsArray[numpy.dtype[Any]]], bool, int, float, complex, str, bytes, numpy._typing._nested_sequence._NestedSequence[bool | int | float | complex | str | bytes]]) -> None:
```

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

#### rmsdiff

```python
def rmsdiff(self, reference):
```

_No documentation provided_

---

<!-- Source: jeremy.py -->
### MoleculeHandler

*Full path: `jeremy.MoleculeHandler`*

_No class documentation provided_

### Properties

#### `atom_labels`

_No documentation provided_

#### `internal_list`

Return copy of the internal coordinate list.

Returns
-------
list

#### `internal_values`

Return copy of the list containing internal coordinate values.

Returns
-------
list

#### `position_array`

Return a copy of the xyz cooridnates array.

Returns
-------
ArrayLike

### Methods

#### clear_cache

```python
def clear_cache(self) -> None:
```

Clear cached properties.

Returns
-------
None

#### compare_internals

```python
def compare_internals(self, reference: 'Jeremy') -> numpy.ndarray:
```

_No documentation provided_

#### override_connectivity_matrix

```python
def override_connectivity_matrix(self, cm: Union[numpy._typing._array_like._SupportsArray[numpy.dtype[Any]], numpy._typing._nested_sequence._NestedSequence[numpy._typing._array_like._SupportsArray[numpy.dtype[Any]]], bool, int, float, complex, str, bytes, numpy._typing._nested_sequence._NestedSequence[bool | int | float | complex | str | bytes]]) -> None:
```

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

#### rmsdiff

```python
def rmsdiff(self, reference):
```

_No documentation provided_

---

<!-- Source: laura.py -->
### Laura

*Full path: `laura.Laura`*

_No class documentation provided_

### Methods

#### generate_xyzfile

```python
def generate_xyzfile(self, filename=None, print_energies=True):
```

_No documentation provided_

#### get_coordinates_dataframe

```python
def get_coordinates_dataframe(self):
```

_No documentation provided_

---

<!-- Source: maxime.py -->
### Jeremy

*Full path: `maxime.Jeremy`*

_No class documentation provided_

### Properties

#### `atom_labels`

_No documentation provided_

#### `internal_list`

Return copy of the internal coordinate list.

Returns
-------
list

#### `internal_values`

Return copy of the list containing internal coordinate values.

Returns
-------
list

#### `position_array`

Return a copy of the xyz cooridnates array.

Returns
-------
ArrayLike

### Methods

#### clear_cache

```python
def clear_cache(self) -> None:
```

Clear cached properties.

Returns
-------
None

#### compare_internals

```python
def compare_internals(self, reference: 'Jeremy') -> numpy.ndarray:
```

_No documentation provided_

#### override_connectivity_matrix

```python
def override_connectivity_matrix(self, cm: Union[numpy._typing._array_like._SupportsArray[numpy.dtype[Any]], numpy._typing._nested_sequence._NestedSequence[numpy._typing._array_like._SupportsArray[numpy.dtype[Any]]], bool, int, float, complex, str, bytes, numpy._typing._nested_sequence._NestedSequence[bool | int | float | complex | str | bytes]]) -> None:
```

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

#### rmsdiff

```python
def rmsdiff(self, reference):
```

_No documentation provided_

---

<!-- Source: maxime.py -->
### Maxime

*Full path: `maxime.Maxime`*

_No class documentation provided_

**Inherits from:** Jeremy

### Properties

#### `atom_labels`

_No documentation provided_

#### `internal_list`

Return copy of the internal coordinate list.

Returns
-------
list

#### `internal_values`

Return copy of the list containing internal coordinate values.

Returns
-------
list

#### `position_array`

Return a copy of the xyz cooridnates array.

Returns
-------
ArrayLike

### Methods

#### clear_cache

```python
def clear_cache(self) -> None:
```

Clear cached properties.

Returns
-------
None

#### compare_internals

```python
def compare_internals(self, reference: 'Jeremy') -> numpy.ndarray:
```

_No documentation provided_

#### complete_molecule

```python
def complete_molecule(self, index):
```

_No documentation provided_

#### complete_molecules

```python
def complete_molecules(self, frontier_df: pandas.core.frame.DataFrame) -> pandas.core.frame.DataFrame:
```

_No documentation provided_

#### extract_radius_atoms

```python
def extract_radius_atoms(self, point: list[float] = [0, 0, 0], r: float = 1, new_filename='None') -> list:
```

_No documentation provided_

#### override_connectivity_matrix

```python
def override_connectivity_matrix(self, cm: Union[numpy._typing._array_like._SupportsArray[numpy.dtype[Any]], numpy._typing._nested_sequence._NestedSequence[numpy._typing._array_like._SupportsArray[numpy.dtype[Any]]], bool, int, float, complex, str, bytes, numpy._typing._nested_sequence._NestedSequence[bool | int | float | complex | str | bytes]]) -> None:
```

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

#### rmsdiff

```python
def rmsdiff(self, reference):
```

_No documentation provided_

---

<!-- Source: misc/NXSpecReader.py -->
### NXSpecReader

*Full path: `NXSpecReader.NXSpecReader`*

Reads and processes spectral data files to generate convoluted spectra.

Parameters
----------
filename : str
    Path to the input spectral data file
npoints : int, optional
    Number of points for spectrum generation (default: 1000)

### Methods

#### generate_spectrum

```python
def generate_spectrum(self) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]:
```

Generate full spectrum across all states.

Returns
-------
NDArray[np.float64]
    Complete convolved spectrum

#### return_semiclassical_spectrum

```python
def return_semiclassical_spectrum(self):
```

_No documentation provided_

#### save_to_file

```python
def save_to_file(self, filename: str) -> None:
```

Save generated spectrum to file.

Parameters
----------
filename : str
    Output file path

---

<!-- Source: rodrigo.py -->
### Rodrigo

*Full path: `rodrigo.Rodrigo`*

_No class documentation provided_

### Methods

#### export_gaussian_eV

```python
def export_gaussian_eV(self, dispersion: float = 0.2):
```

_No documentation provided_

#### get_peaks_eV

```python
def get_peaks_eV(self):
```

_No documentation provided_

#### get_peaks_nm

```python
def get_peaks_nm(self):
```

_No documentation provided_

#### load_experimental

```python
def load_experimental(self, experimental_filename, units='eV'):
```

_No documentation provided_

#### plot

```python
def plot(self):
```

_No documentation provided_

#### plot_additional_spectra

```python
def plot_additional_spectra(self, data, units='eV', label='Imported Spectrum', shift=0.0):
```

_No documentation provided_

#### plot_gaussian

```python
def plot_gaussian(self, dispersion=None):
```

_No documentation provided_

#### plot_vertical

```python
def plot_vertical(self):
```

_No documentation provided_

#### save_image

```python
def save_image(self, imagename):
```

_No documentation provided_

---

<!-- Source: stef.py -->
### NebExtractor

*Full path: `stef.NebExtractor`*

_No class documentation provided_

### Methods

#### export_data

```python
def export_data(self, name: str = 'default.dat'):
```

_No documentation provided_

#### savefig

```python
def savefig(self, name: str = 'default.jpg', dpi: int = 600, landscape: bool = False, spline=False, annotate_barrier: bool = False):
```

_No documentation provided_

---
