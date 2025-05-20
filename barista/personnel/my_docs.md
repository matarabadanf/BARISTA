# 📘 Python Class Documentation (Alphabetical Order)

## 📚 Table of Contents

- [Alberto](#alberto)
- [Aurora](#aurora)
- [AuroraConfiguration](#auroraconfiguration)
- [Dani](#dani)
- [Emma](#emma)
- [Fonsi](#fonsi)
- [Javi](#javi)
- [Jeremy](#jeremy)
- [Jeremy](#jeremy)
- [Laura](#laura)
- [Maxime](#maxime)
- [NebExtractor](#nebextractor)
- [NXSpecReader](#nxspecreader)
- [Rodrigo](#rodrigo)

---

## 🧱 Class: `Alberto` <a name="alberto"></a>
_Defined in: `alberto.py`_

Takes an excited state optimization in ORCA and plots the energy of each state at each
step of the optimization and the actual root.

This class parses ORCA output files from TD-DFT/TDA calculations to track and
visualize the energy evolution of electronic states during geometry optimization.

Attributes
 ---
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

### Methods:
- **`__init__`**: Initialize the Alberto class with the specified parameters.

**Parameters**
 ---
filename : str
    Path to the ORCA output file to analyze.
reference_energy : float, optional
    Reference energy value in Hartree, by default 0.0.
    Used to calculate relative energies.
units : str, optional
    Energy units for plotting ("eV" or "Hartree"), by default "eV".
- **`_initialize`**: Initialize all data structures and run the parsing methods in sequence.

This method coordinates the full initialization process by calling all
the necessary parsing and processing methods in the correct order.
- **`_get_file_content`**: Open the ORCA output file and store its content as a list of strings.

Raises
------
ValueError
    If there is no filename specified.
FileNotFoundError
    If the file does not exist.
- **`_get_state_indexes`**: Extract state indexes from the ORCA output file.

Finds all states mentioned in the TD-DFT/TDA section and stores their indices.
Sets the number_of_states attribute based on discovered states.

Raises
------
ValueError
    If no result file was found.
- **`_get_cis_energies`**: Extract SCF energies from the ORCA output file.

Collects all SCF energies from the file and stores them in cis_array.
If no reference energy was provided, the first SCF energy is used.
Sets the number_of_steps attribute based on collected energies.
- **`_get_abs_energy_array`**: Parse the ORCA output file to obtain the absolute energies in a.u.

The procedure performed is the following:
 - Extract the energy increase from the current s0 to the excited states
 - Set the first row of the energy array to the SCF energy
 - Add the energy of the excited states to the SCF in each state

This returns the total absolute energy of each of the states
(typically in the test energies between -562.25 and -561.85 a.u.).
- **`_get_relative_energy`**: Calculate energies relative to the reference energy.

Subtracts the reference energy from all energy values to obtain
relative energies that are more suitable for visualization.
- **`_get_current_energy`**: Extract the current energy surface at each optimization step.

Identifies which electronic state is being followed during the optimization
at each step and extracts its energy. Also determines the starting and final 
roots used in the optimization.
- **`set_units`**: Change the energy units used for plotting.

**Parameters**
 ---
units : str
    Energy units to use ("eV" or "Hartree").
    
Raises
------
ValueError
    If an unsupported unit is specified.
- **`_rescale_units`**: Rescale energy values according to the selected units.

Converts energy values between eV and Hartree using the conversion
factor 27.2114 eV/Hartree.

**Parameters**
 ---
prev_units : str
    Previous units used for energy values.
- **`set_imagename`**: Sets the output image name for saving plots.

**Parameters**
 ---
image_name : str
    Image name to set.
- **`_prepare_plot`**: Prepare the energy plot showing all states and the current surface.

Creates a figure with plots for all electronic states and marks the
current surface being followed during optimization.

**Returns**
 
Tuple[plt.Figure, plt.Axes]
    The created figure and axes objects.
- **`plot`**: Display the energy plot showing all states and the current surface.

Creates and displays a plot showing the energy evolution of all states
during the optimization, highlighting the current surface being followed.
- **`generate_image`**: Save the energy plot to an image file.

**Parameters**
 ---
output_image_name : str
    Name of the output image file.
- **`starting_root`**: Get the initial root used for the optimization.

**Returns**
 
int
    Index of the starting root.
- **`final_root`**: Get the final root after optimization.

**Returns**
 
int
    Index of the final root.

---

## 🧱 Class: `Aurora` <a name="aurora"></a>
_Defined in: `aurora.py`_

A class for extracting and analyzing potential energy surface (PES) paths.

This class processes xyz files containing energy data to extract PES branch information,
calculate relative energies, and visualize energy paths along different electronic states.

**Parameters**
 ---
fc_filename : str
    Path to the Franck-Condon xyz file with energy as comment
last_branch_filename : str
    Path to the xyz file of the last point in the branch

### Methods:
- **`__init__`**: _No docstring_
- **`from_files`**: Factory method to create and initialize an Aurora instance.

**Parameters**
 ---
fc_file : str
    Path to the Franck-Condon xyz file
last_branch_file : str
    Path to the xyz file of the last point in the branch
    
**Returns**
 
Aurora
    An initialized Aurora instance
- **`get_relative_energies`**: Obtain energy values of an xyz file that has energy values as comment string.

This method reads the xyz file, extracts the energy values from the comment line,
and calculates relative energies with respect to the reference energy in eV.

**Parameters**
 ---
xyz_file : str
    Path to the xyz file containing energy data
    
**Returns**
 
numpy.ndarray
    Array of relative energies in eV
    
Raises
------
ValueError
    If the reference energy has not been determined
- **`pes_DataFrame`**: Get a copy of the PES dataframe.

**Returns**
 
pandas.DataFrame
    A copy of the dataframe containing all PES data
- **`plot_branch`**: Plot the PES branch and display the figure.

This method calls _preplot_branch to prepare the plot and then displays it.

**Returns**
 
None
- **`save_branch_plot`**: Save the PES branch plot to a file.

**Parameters**
 ---
name : str, optional
    The output filename, if None the last branch filename with jpg extension will be used
    
**Returns**
 
None
- **`_initialize`**: Initialize all required data before calculations.

This method calls several private methods to set up the reference energy,
get name and identifiers, determine filenames, and fill the dataframe.

**Returns**
 
None
- **`_get_name_and_identifiers`**: Extract the prefix name and identifiers from the last branch filename.

This method parses the last branch filename to extract the pre-identifier name
and the step identifiers that define the path through the PES.

**Returns**
 
tuple
    Tuple containing the pre-identifier name and step identifiers
- **`_get_filenames`**: Generate filenames for all points in the branch.

This method creates a list of filenames for all points in the branch
based on the pre-identifier name and step identifiers.

**Returns**
 
list
    List of filenames for all points in the branch
    
Raises
------
ValueError
    If the general name or path step identifiers have not been determined
- **`_set_reference_energy`**: Set the reference energy from the Franck-Condon state.

This method reads the FC file, extracts the energies, sets the reference energy
to the ground state energy, and initializes the PES dataframe with FC energies.

**Returns**
 
float
    The reference energy
- **`_fill_dataframe`**: Fill the dataframe with all the energies and steps of the branch.

This method iterates through all filenames, extracts the energies,
and populates the PES dataframe with the relative energies for each state and step.

**Returns**
 
None
- **`_preplot_branch`**: Prepare the plot for the PES branch.

This method sets up the plot with appropriate labels, plots all states,
and highlights the path through the PES.

**Returns**
 
matplotlib.pyplot
    The prepared plot
- **`_get_path_energies`**: Get the energies along the highlighted path.

This method extracts the energy values for each step along the highlighted path
through the PES.

**Returns**
 
tuple
    Tuple containing x-coordinates and corresponding energy values
- **`_get_plot_labels`**: Generate x-axis labels for the plot.

This method creates appropriate labels for each step in the PES path.

**Returns**
 
tuple
    Tuple containing x-coordinates and corresponding tick labels
- **`_get_highlighted_path`**: Get the highlighted path through the PES.

This method determines the sequence of electronic states that form
the path through the PES.

**Returns**
 
list
    List of state indices forming the highlighted path

---

## 🧱 Class: `AuroraConfiguration` <a name="auroraconfiguration"></a>
_Defined in: `aurora.py`_

Configuration class for Aurora with default parameters.

**Parameters**
 ---
fc_file : str, optional
    Path to the Franck-Condon xyz file with energy as comment, default is None
last_branch_file : str, optional
    Path to the xyz file of the last point in the branch, default is None

_No methods found._

---

## 🧱 Class: `Dani` <a name="dani"></a>
_Defined in: `dani.py`_

A class for extracting and processing data from Gaussian output files.

This class extracts final geometry coordinates and vibrational frequencies 
from Gaussian output files, and can generate XYZ format files.

Attributes
 ---
gaussian_output_filename : str
    Path to the Gaussian output file
atomic_dict : Dict[int, str]
    Dictionary mapping atomic numbers to element symbols

### Methods:
- **`__init__`**: Initialize the Dani class with a Gaussian output file.

**Parameters**
 ---
gaussian_output_filename : str
    Path to the Gaussian output file
- **`from_file`**: Factory method to create and initialize a Dani instance from a file.

**Parameters**
 ---
gaussian_output_filename : str
    Path to the Gaussian output file
    
**Returns**
 
Dani
    An initialized Dani instance
- **`generate_xyzfile`**: Generate an XYZ format file from the extracted coordinates.

**Parameters**
 ---
filename : str, optional
    Output filename, if None the input filename with .xyz extension will be used
    
**Returns**
 
None
- **`extract_vibrational_modes`**: Extract vibrational frequencies and IR intensities from the Gaussian output.

**Parameters**
 ---
filename : str, optional
    Output filename, if None the input filename with _frequencies.dat extension will be used
    
**Returns**
 
None
- **`_initialize`**: Initialize the object by reading the file and extracting coordinates.

**Returns**
 
None
- **`_get_file_content`**: Read the content of the Gaussian output file.

**Returns**
 
None

Raises
------
FileNotFoundError
    If the specified file cannot be found
- **`_get_coordinates_section`**: Extract the section containing the final coordinates from the Gaussian output file.

This method searches backwards through the file to find the last occurrence of
the coordinates section, which should contain the final optimized geometry.

**Returns**
 
tuple
    A tuple containing the start and end indices of the coordinates section

---

## 🧱 Class: `Emma` <a name="emma"></a>
_Defined in: `emma.py`_

_No docstring provided._

### Methods:
- **`__init__`**: _No docstring_
- **`_determine_targets`**: _No docstring_
- **`_get_gs_energy`**: _No docstring_
- **`_generate_geometries`**: _No docstring_
- **`_generate_geometry_objects`**: _No docstring_
- **`internal_deviation`**: _No docstring_
- **`_generate_dataframe`**: _No docstring_
- **`rmsd_report_to_md`**: _No docstring_

---

## 🧱 Class: `Fonsi` <a name="fonsi"></a>
_Defined in: `fonsi.py`_

_No docstring provided._

### Methods:
- **`__init__`**: _No docstring_
- **`_initialize`**: _No docstring_
- **`_load_energies`**: _No docstring_
- **`_preplot`**: _No docstring_
- **`save_fig`**: _No docstring_
- **`plot`**: _No docstring_

---

## 🧱 Class: `Javi` <a name="javi"></a>
_Defined in: `javi.py`_

_No docstring provided._

### Methods:
- **`__init__`**: _No docstring_
- **`_init`**: _No docstring_
- **`from_xy`**: _No docstring_
- **`ci_energy`**: _No docstring_
- **`_load_vectors`**: _No docstring_
- **`g_ab`**: _No docstring_
- **`s_ab`**: _No docstring_
- **`h_ab`**: _No docstring_
- **`_pre_beta`**: _No docstring_
- **`_rotate_for_beta`**: _No docstring_
- **`beta`**: _No docstring_
- **`is_rotation_needed`**: _No docstring_
- **`_g_tilde`**: _No docstring_
- **`_h_tilde`**: _No docstring_
- **`x`**: _No docstring_
- **`y`**: _No docstring_
- **`pitch`**: Pitch \delta_gh. 
- **`asymmetry`**: Asymmetry \Delta_gh. 
- **`energy_difference`**: _No docstring_
- **`average_energy`**: _No docstring_
- **`E_A`**: _No docstring_
- **`E_B`**: _No docstring_
- **`theta_s`**: _No docstring_
- **`sigma`**: _No docstring_
- **`p`**: _No docstring_
- **`b`**: _No docstring_
- **`plot_CI`**: _No docstring_
- **`generate_force_file`**: _No docstring_

---

## 🧱 Class: `Jeremy` <a name="jeremy"></a>
_Defined in: `jeremy.py`_

_No docstring provided._

### Methods:
- **`__init__`**: _No docstring_
- **`override_connectivity_matrix`**: Override connectivity matrix to generate internals from a reference.

Redefines the internal coordinates based in the reference connectivity.
Useful when comparing two structures, being able to determine their
'likeness' by obtaining the same internals of both and their values.

**Parameters**
 ---
connectivity_matrix : ArrayLike
    Connectivity matrix of the reference system.

Raises
------
ValueError
    DESCRIPTION.

**Returns**
 
None.
- **`clear_cache`**: Clear cached properties.

**Returns**
 
None
- **`rmsdiff`**: _No docstring_
- **`compare_internals`**: _No docstring_
- **`position_array`**: Return a copy of the xyz cooridnates array.

**Returns**
 
ArrayLike
- **`atom_labels`**: _No docstring_
- **`internal_list`**: Return copy of the internal coordinate list.

**Returns**
 
list
- **`internal_values`**: Return copy of the list containing internal coordinate values.

**Returns**
 
list
- **`connectivity_matrix`**: Build connectivity matrix.

If below the threshold, will fill with 1. If no connection below the
threshold, it will search for the nearest atom and assign the value 10
in order to leave no atom unconnected.


**Parameters**
 ---
bond_thresh : float, optional
    Bond threshold. The default is 1.6.

**Returns**
 
connectivity_matrix : ArrayLike
    Connectivity matrix of the molecule.
- **`xyz_df`**: Return copy of the cartesian dataframe.

**Returns**
 
pd.DataFrame
- **`distance_matrix`**: _No docstring_
- **`r_vector_matrix`**: _No docstring_
- **`bond_list`**: _No docstring_
- **`angle_list`**: Build the angle list.

Done by adding an atom to the end of the bonds forwards and backwards.
If the angle is new it is stored. If not, discarded.

**Returns**
 
list
    Angle list with atomic indices.
- **`dihedral_list`**: Build the dihedral list.

It is done by adding an atom to the end of the angles forwards and
backwards. If the dihedral is new it is stored. If not, discarded
- **`bond_values`**: _No docstring_
- **`angle_values`**: Calculate the angles.

Uses theta = arccos((r_ba · r_bc)/(||r_ba||*||r_bc||)) * 180/pi.
- **`dihedral_values`**: Calculate the dihedral angles.

Done with the atan(sin_phi, cos_phi) formula.
- **`_read_xyzfile`**: _No docstring_
- **`_extract_atoms`**: _No docstring_
- **`_extract_positions`**: _No docstring_
- **`_extend_labels`**: _No docstring_

---

## 🧱 Class: `Jeremy` <a name="jeremy"></a>
_Defined in: `maxime.py`_

_No docstring provided._

### Methods:
- **`__init__`**: _No docstring_
- **`override_connectivity_matrix`**: Override connectivity matrix to generate internals from a reference.

Redefines the internal coordinates based in the reference connectivity.
Useful when comparing two structures, being able to determine their
'likeness' by obtaining the same internals of both and their values.

**Parameters**
 ---
connectivity_matrix : ArrayLike
    Connectivity matrix of the reference system.

Raises
------
ValueError
    DESCRIPTION.

**Returns**
 
None.
- **`clear_cache`**: Clear cached properties.

**Returns**
 
None
- **`rmsdiff`**: _No docstring_
- **`compare_internals`**: _No docstring_
- **`position_array`**: Return a copy of the xyz cooridnates array.

**Returns**
 
ArrayLike
- **`atom_labels`**: _No docstring_
- **`internal_list`**: Return copy of the internal coordinate list.

**Returns**
 
list
- **`internal_values`**: Return copy of the list containing internal coordinate values.

**Returns**
 
list
- **`connectivity_matrix`**: Build connectivity matrix.

If below the threshold, will fill with 1. If no connection below the
threshold, it will search for the nearest atom and assign the value 10
in order to leave no atom unconnected.


**Parameters**
 ---
bond_thresh : float, optional
    Bond threshold. The default is 1.6.

**Returns**
 
connectivity_matrix : ArrayLike
    Connectivity matrix of the molecule.
- **`xyz_df`**: Return copy of the cartesian dataframe.

**Returns**
 
pd.DataFrame
- **`distance_matrix`**: _No docstring_
- **`r_vector_matrix`**: _No docstring_
- **`bond_list`**: _No docstring_
- **`angle_list`**: Build the angle list.

Done by adding an atom to the end of the bonds forwards and backwards.
If the angle is new it is stored. If not, discarded.

**Returns**
 
list
    Angle list with atomic indices.
- **`dihedral_list`**: Build the dihedral list.

It is done by adding an atom to the end of the angles forwards and
backwards. If the dihedral is new it is stored. If not, discarded
- **`bond_values`**: _No docstring_
- **`angle_values`**: Calculate the angles.

Uses theta = arccos((r_ba · r_bc)/(||r_ba||*||r_bc||)) * 180/pi.
- **`dihedral_values`**: Calculate the dihedral angles.

Done with the atan(sin_phi, cos_phi) formula.
- **`_read_xyzfile`**: _No docstring_
- **`_extract_atoms`**: _No docstring_
- **`_extract_positions`**: _No docstring_
- **`_extend_labels`**: _No docstring_

---

## 🧱 Class: `Laura` <a name="laura"></a>
_Defined in: `laura.py`_

_No docstring provided._

### Methods:
- **`__init__`**: _No docstring_
- **`from_file`**: _No docstring_
- **`_initialize`**: _No docstring_
- **`converged`**: _No docstring_
- **`_get_file_content`**: _No docstring_
- **`_get_nroots`**: _No docstring_
- **`_is_gradient_calculation`**: _No docstring_
- **`_get_gs_energy`**: _No docstring_
- **`_get_excitation_energies`**: _No docstring_
- **`_get_coordinates_section`**: _No docstring_
- **`_extract_coordinates_to_dataframe`**: _No docstring_
- **`get_coordinates_dataframe`**: _No docstring_
- **`generate_xyzfile`**: _No docstring_

---

## 🧱 Class: `Maxime` <a name="maxime"></a>
_Defined in: `maxime.py`_

_No docstring provided._

### Methods:
- **`extract_radius_atoms`**: _No docstring_
- **`complete_molecules`**: _No docstring_
- **`complete_molecule`**: _No docstring_

---

## 🧱 Class: `NebExtractor` <a name="nebextractor"></a>
_Defined in: `stef.py`_

_No docstring provided._

### Methods:
- **`__init__`**: _No docstring_
- **`file_content`**: _No docstring_
- **`energy_array`**: _No docstring_
- **`forward_barrier`**: _No docstring_
- **`reverse_barrier`**: _No docstring_
- **`absolute_barrier`**: _No docstring_
- **`_preplot`**: _No docstring_
- **`savefig`**: _No docstring_
- **`export_data`**: _No docstring_

---

## 🧱 Class: `NXSpecReader` <a name="nxspecreader"></a>
_Defined in: `misc/NXSpecReader.py`_

Reads and processes spectral data files to generate convoluted spectra.

**Parameters**
 ---
filename : str
    Path to the input spectral data file
npoints : int, optional
    Number of points for spectrum generation (default: 1000)

### Methods:
- **`__init__`**: _No docstring_
- **`_initialize`**: Initialize all necessary data structures and perform initial calculations.
- **`_read_file`**: Read the input file and store its contents.
- **`_get_max_state`**: Determine the highest excited state number in the file.
- **`_generate_dictionaries`**: Initialize dictionaries to store energy and intensity data for each state.
- **`_extract_states`**: Extract energy and oscillator strength values for each excited state.
- **`_transform_dictionaries_to_df`**: Convert energy and intensity dictionaries to a pandas DataFrame.

**Returns**
 
pd.DataFrame
    DataFrame containing state, energy, and oscillator strength data
- **`_gaussian`**: Generate a Gaussian function.

**Parameters**
 ---
X : NDArray[np.float64]
    X values for Gaussian calculation
height : float, optional
    Peak height (default: 1)
center : float, optional
    Peak center (default: 15)
spread : float, optional
    Peak width (default: 0.05)

**Returns**
 
NDArray[np.float64]
    Gaussian function values
- **`_generate_x_grid`**: Generate energy grid for spectrum calculation.
- **`_generate_xgrid_nm`**: Generate wavelength grid in nanometers.
- **`_convolute_spectrum`**: Generate convolved spectrum for given states.

**Parameters**
 ---
filtered_df : pd.DataFrame
    DataFrame containing energy and oscillator strength values

**Returns**
 
NDArray[np.float64]
    Convolved spectrum
- **`generate_spectrum`**: Generate full spectrum across all states.

**Returns**
 
NDArray[np.float64]
    Complete convolved spectrum
- **`return_semiclassical_spectrum`**: _No docstring_
- **`save_to_file`**: Save generated spectrum to file.

**Parameters**
 ---
filename : str
    Output file path

---

## 🧱 Class: `Rodrigo` <a name="rodrigo"></a>
_Defined in: `rodrigo.py`_

_No docstring provided._

### Methods:
- **`__init__`**: _No docstring_
- **`_initialize`**: _No docstring_
- **`_load_file`**: _No docstring_
- **`_locate_peaks_in_file`**: _No docstring_
- **`_extract_peaks`**: _No docstring_
- **`_normalize_peaks`**: _No docstring_
- **`get_peaks_nm`**: _No docstring_
- **`get_peaks_eV`**: _No docstring_
- **`_gaussian`**: Generate a Gaussian function.

**Parameters**
 ---
X : NDArray[np.float64]
    X values for Gaussian calculation
height : float, optional
    Peak height (default: 1)
center : float, optional
    Peak center (default: 15)
spread : float, optional
    Peak width (default: 0.05)

**Returns**
 
NDArray[np.float64]
    Gaussian function values
- **`_plot_gaussian_nm`**: _No docstring_
- **`export_gaussian_eV`**: _No docstring_
- **`_plot_gaussian_eV`**: _No docstring_
- **`plot_gaussian`**: _No docstring_
- **`load_experimental`**: _No docstring_
- **`_plot_vertical_nm`**: _No docstring_
- **`_plot_vertical_eV`**: _No docstring_
- **`plot_vertical`**: _No docstring_
- **`plot_additional_spectra`**: _No docstring_
- **`plot`**: _No docstring_
- **`save_image`**: _No docstring_

---

