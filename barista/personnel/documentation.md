# Command Line Tools Documentation

This document provides documentation for command line tools.

## Table of Contents

* [dani](#dani)
* [fonsi](#fonsi)
* [aurora](#aurora)
* [laura](#laura)
* [javi](#javi)
* [rodrigo](#rodrigo)
* [lorena](#lorena)
* [alberto](#alberto)
* [emma](#emma)
* [stef](#stef)
* [maxime](#maxime)

---

<!-- Source: dani.py -->
## dani

Extract the final geometry of an optimization run in Gaussian. Extract vibrational frequencies if available

### Usage
```
generate_docs.py [-h] -f F [-o O] [--frequencies]
```

### Arguments

#### Options

* `-f` (str) [required]: Input file name
  * Default: `42`

* `-o` (str): Output file name. Default is Input_filename.xyz

* `--frequencies`: Extract frequencies to Input_filename_freq.dat
  * Default: `False`


---

<!-- Source: fonsi.py -->
## fonsi

Takes a brewer optimizaiton energies.dat and plots the energy at each step.

### Usage
```
generate_docs.py [-h] -f F [-o O] [-en EN] [-u U] [--interactive]
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

<!-- Source: aurora.py -->
## aurora

Exctract the potential energy surface path and energies of of a PES branch.
    Takes as input the ground state equilibrium geometry xyz file  and the xyz file of the last point of the PES branch.
    
    

This script determines the path, highlights, and generates a pandas dataframe that can be later plotted and saved as an image.

### Usage
```
generate_docs.py [-h] -f F -b B [-o O]
```

### Arguments

#### Options

* `-f` (str) [required]: FC xyz file with energy as comment.
  * Default: `42`

* `-b` (str) [required]: Last point in the branch.

* `-o` (str): Output image name.


---

<!-- Source: laura.py -->
## laura

Exctract the final geometry of an optimization run in ORCA.

### Usage
```
generate_docs.py [-h] -f F [-o O] [--no_energies]
```

### Arguments

#### Options

* `-f` (str) [required]: Input file name
  * Default: `42`

* `-o` (str): Output file name

* `--no_energies`: Dont include energies in the xyzfile as comment.
  * Default: `True`


---

<!-- Source: javi.py -->
## javi

Takes the lower and upper state gradients in a CI with the NAC vector and characterizes the type oc CI.        
Results can be plotted interactively using plotly.

### Usage
```
generate_docs.py [-h] -g0 G0 -g1 G1 -nac NAC [--interactive]
```

### Arguments

#### Options

* `-g0` (str) [required]: Lower state gradient file.

* `-g1` (str) [required]: Higher state gradient file.

* `-nac` (str) [required]: NAC vector file.

* `--interactive`: Open interactive plot.
  * Default: `False`


---

<!-- Source: rodrigo.py -->
## rodrigo

Generates the TDDFT absorption spectre from an ORCA TDDFT calculation.         It generates spectre .jpg image.
         It can generate a peaks.dat and a spectra.dat that contains the peaks and spectre.

### Usage
```
generate_docs.py [-h] -F F [-o O] [-u U] [--gaussian] [-gauss_disp GAUSS_DISP] [-exp EXP]
                        [-exp_units EXP_UNITS] [-exp_shift EXP_SHIFT] [-semi SEMI] [--interactive]
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

<!-- Source: lorena.py -->
## lorena

Takes the an xyz file and shakes the geometry

### Usage
```
generate_docs.py [-h] [-f F] [-d D] [-o O]
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

<!-- Source: alberto.py -->
## alberto

Takes a excited state optimization in ORCA and plots the energy of each state at each step of the opitmizaiton and the actual root.

### Usage
```
generate_docs.py [-h] -f F [-o O] [-en EN] [-u U]
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

<!-- Source: emma.py -->
## emma

Takes the ground state geometry of a molecule and compares it to optimized TDDFT geometries run in ORCA.        
It compares the energy, RMSD and final root, plotting the results in an image. 

### Usage
```
generate_docs.py [-h] -fc FC -ex EX [EX ...] [-O O] [-o O] [-r R] [-md MD] [-i I]
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

<!-- Source: stef.py -->
## stef

Takes the ground state geometry of a molecule and compares it to optimized TDDFT geometries run in ORCA.        
It compares the energy, RMSD and final root, plotting the results in an image. 

### Usage
```
generate_docs.py [-h] -f F [-en EN] [-u U] [-o O] [-dat DAT] [--landscape] [--spline]
                        [--barrier]
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

<!-- Source: maxime.py -->
## maxime

Extract atoms from an xyz file with a center point and a radius.

### Usage
```
generate_docs.py [-h] -f F [-o O] [-r R] [-p P]
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
