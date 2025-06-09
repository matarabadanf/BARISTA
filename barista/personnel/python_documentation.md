# Python Code Documentation

This document provides comprehensive documentation for Python scripts.

## Table of Contents

### Command Line Interfaces

* [alberto](#alberto)
* [aurora](#aurora)
* [cris](#cris)
* [dani](#dani)
* [emma](#emma)
* [fonsi](#fonsi)
* [javi](#javi)
* [laura](#laura)
* [lorena](#lorena)
* [maxime](#maxime)
* [rodrigo](#rodrigo)
* [stef](#stef)

# Command Line Interfaces

<!-- Source: alberto.py -->
## CLI Arguments: alberto

Takes a excited state optimization in ORCA and plots the energy of each state at each step of the opitmizaiton and the actual root.

### Usage
```
doc2.py [-h] -f F [-o O] [-en EN] [-u U]
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
doc2.py [-h] -f F -b B [-o O]
```

### Arguments

#### Options

* `-f` (str) [required]: FC xyz file with energy as comment.
  * Default: `42`

* `-b` (str) [required]: Last point in the branch.

* `-o` (str): Output image name.


---

<!-- Source: cris.py -->
## CLI Arguments: cris

Takes the lower and upper state gradients in a CI with the CDV vector and characterizes the type oc CI.        
Results can be plotted interactively using plotly.

### Usage
```
doc2.py [-h] [-molcas_file MOLCAS_FILE] [-g0 G0] [-g1 G1] [-cdv CDV] [-ref_xyz REF_XYZ] [--generate_forces] [--generate_displacements] [-disp_mod DISP_MOD] [--plot_polar] [--plot_lower] [--plot_upper]
               [--interactive]
```

### Arguments

#### 
Core Data

* `-molcas_file` (str): Molcas .out file with gradients and cdv.
  * Default: `0.0`

* `-g0` (str): Lower state gradient file.

* `-g1` (str): Higher state gradient file.

* `-cdv` (str): cdv vector file.


#### 
Vector Visualization And Geometry Displacement

* `-ref_xyz` (str): Reference xyz file for forces and displacement.
  * Default: `None`

* `--generate_forces`: Generate xyz file with  x, y, and -s_ab vectors for visualization.
  * Default: `False`

* `--generate_displacements`: Generate displaced geometries in x, -x, y, -y -s_ab and mirrored -s_ab. A folder will be generated with these geometries.
  * Default: `False`

* `-disp_mod` (float): Displacement module of the vectors to generate displaced geometries. Default is 0.1.
  * Default: `0.1`


#### 
Pes Visualization

* `--plot_polar`: Plot surfaces in polar coordinates in 2D.
  * Default: `False`

* `--plot_lower`: Plot lower state pes in 2D.
  * Default: `False`

* `--plot_upper`: Plot upper state pes in 2D.
  * Default: `False`

* `--interactive`: Open interactive plot.
  * Default: `False`


---

<!-- Source: dani.py -->
## CLI Arguments: dani

Extract the final geometry of an optimization run in Gaussian. Extract vibrational frequencies if available

### Usage
```
doc2.py [-h] -f F [-o O] [--frequencies]
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
doc2.py [-h] -fc FC -ex EX [EX ...] [-O O] [-o O] [-r R] [-md MD] [-i I]
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
doc2.py [-h] -f F [-o O] [-en EN] [-u U] [--interactive]
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
doc2.py [-h] -g0 G0 -g1 G1 -nac NAC [--interactive]
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
doc2.py [-h] -f F [-o O] [--no_energies]
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
doc2.py [-h] [-f F] [-d D] [-o O]
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
doc2.py [-h] -f F [-o O] [-r R] [-p P]
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
doc2.py [-h] -F F [-o O] [-u U] [--gaussian] [-gauss_disp GAUSS_DISP] [-exp EXP] [-exp_units EXP_UNITS] [-exp_shift EXP_SHIFT] [-semi SEMI] [--interactive]
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
doc2.py [-h] -f F [-en EN] [-u U] [-o O] [-dat DAT] [--landscape] [--spline] [--barrier]
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
