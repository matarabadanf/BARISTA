# BARISTA 

A collection of scripts and small programs to assist in PES exploration in excited states. BARISTA is divided in two main funcitonalities. BREWER is a connical intersection optimizer based in the gradient projection and updated branching plane methods that is modular and can be interfaced with electronic structure programs. Personnel contains a collection of modules that can be executed as command-line scripts or imported for python scripting with various functionalities.  

## Installation
BARISTA can be installed with 

```bash
pip install -e . 
```

To ensure that the scripts can be executed, use:
```Bash
chmod 700 config/setup_barista.py
setup_barista.sh
```
Which will set the scripts as executables. In order to use them anywhere, add the `/personnel` and `/brewer` folders to the path. 

To configure the electronic structure launchers, use:
```bash
config/setup_launchers.py
```
This sets up the electronic structure program launchers in a batch system such as `run_orca5.0.3` for the BREWER program.  

# Brewer
BREWER is a small homebrew inteface with electronic structure programs (currently only ORCA) based in the python module `ase.io`. Via `ase.io` it is possible to interface electronic structure calculations with some highly customizable algorithms. 

Currently supported calcuations are:
- MECI search with TDDFT in ORCA.

## MECI calculations

MECI calculations can be performed with BREWER using the following algorithms: 
 - Penalty algorithm.
 - Projection method (no CDV).
 - Updated branching plane (Estimated CDV).

The first is implemented as in [J. Phys. Chem. B 2008, 112, 2, 405–413](https://pubs.acs.org/doi/10.1021/jp0761618). The second and third are implemented as in [J. Chem. Theory Comput. 2010, 6, 1538–1545](https://pubs.acs.org/doi/abs/10.1021/ct1000268).

### Configuring Brewer
Brewer needs to be configured: it needs the path to a suitable ORCA (currently only orca, could be expanded) launcher. To set up the launcher scripts, use `setup_launchers.py`. The `setup_launchers.py` script can be run as such, or as:
```Bash
setup_launchers.py < configuration.dat
```

Where the configuration file should look like:
```
~/bin/run_orca.sh ! orca launcher
```

**IMPORTANT**: 
- ORCA launchers should be able to copy `.xyz` and `.gbw` files to the `$WorkDir`. We include an example launcher for ORCA in the `data/` folder. This is in order to avoid recoputing the SCF optimization in both states, eventually saving some computation time. 

## Using Brewer
Brewer requires an input file named `brewer.TEMPLATE` input file. The `barista/brewer/brewer_setup.py` script generates a directory with the `brewer.TEMPLATE` input file and copies the required files for the calculation with:

```Bash
brewer_setup.py -f geometry_to_optimize.xyz 
```

Along the files that are in the directory there is: 
- `opt.py`: Main function, manages reading the input and performs the optimization
- `flavs.py`: Contains the calculator called by ASE. ASE requires the use of calculators to obtain properties such as forces and energies.
- `brewerjobmanager.py`: manages the input and optimization. 
- `run_opt.sh`: Main script for the queue system.
- `javi.py`: manages MECI characterization after optimization. 

To run a calculation with BREWER in a SLURM queue system, it can be done as:
```Bash
sbatch -args run_opt.sh
```

Consider that the calculation will call an electronic structure program in the following way:

```Bash
/path/to/launcher/laucher.sh input.in
```

Launchers should have been already configured with `setup_launchers.py`. 

After the calculation has finished, a folder with auxiliary files `AUX_FILES/` such as electronc structure program outputs will generated.


## MECI calculations
To perform meci calculations, a `brewer.TEMPLATE` file containing the following fields would suffice:
```
name = molecule.xyz
mode = ci
```
Where the `molecule.xyz` is the starting geometry with and `mode` is the type of job. There are some optional arguments for the CI optimization:
#### Keywords 
| Keyword    | Information|
|------------|------------|
| `profile`  | |
| `iroot`    |           |
| `jroot`    |         |
| `n_roots`  |      |
| `alpha`    |        |
| `sigma`    |          |

| Keyword    | Information|
|------------|------------|
| `profile`       |  MECI search algoritm. Options are `penalty`, `gp` (Gradient projection) and `ubp` (Updated Branching plane).                    |
| `iroot`         |  Lower state potential energy surface. Default is $0$.                     |
| `jroot`         |  Higher state potential energy surface. Default is `iroot`$+1$.                       |
| `n_roots`       |  Number of excited states calculated by the electronic structure program. Default is $10$.                      |
| `basis`         |  Basis for calculation in electronic structure program. Default is cc-pVDZ                    |
| `charge`        |  Charge of the system.                    |
| `mult`          |  Multiplicity of the system.                    |

Choice of electronic structure program.

| Keyword    | Information|
|------------|------------|
| `program`       |  Choice of electronic structure program. Default is ORCA.                    |

In TDDFT:

| Keyword    | Information|
|------------|------------|
| `functional`    |  Coice of cuntional. Default is CAM-B3LYP.                    |

Penalty algorithm:
| Keyword    | Information|
|------------|------------|
| `alpha`         |  Alpha value for penalty method (see reference). Default is $3.5$.                        |
| `sigma`         |  Sigma value for penalty method (see reference). Default is $0.02$.                     |

UBP and GP algorithms
| Keyword    | Information|
|------------|------------|
| `slow_start`    | Add a weight ion the first steps to rescale the forces in order to ensure a smooth start. Default is False.                    |
| `slow_increase` | Increase the weight periodically to ensure smooth transition after first steps.  Default is False.                       |
| `convergence`   | Type of convergence. `Natural` converges when the force is strictly under the threshold. `Forced` consideres converged after 10 steps with a energy deviation smaller than  $0.00001$ Hartree.                   |

MOLCAS interface (work in progress):
| Keyword    | Information|
|------------|------------|
| `rassccf_params`|                      |
| `caspt2_params` |                      |
| `imag`          |                      |

The default `brewer.TEMPLATE` configuration is:

```
geom= cytisine_min_1_1_1.xyz

mode = ci
profile = ubp
calc_nacme = False

functional = CAM-B3LYP
basis = cc-pvdz
iroot = 0
jroot = 1
nroots = 10

charge = 0
mult = 1
program = 'ORCA'
label = None
convergence = 'natural'

slow_start = True
```

## Results
Once a MECI calculation has finished, The resulting files are:
- `brewer_opt.log`: contains the optimizer steps, energy difference in Hartree and max_force at each step. Includes conical intersection characterization at that point. 
- `energies.dat`: shows the energy of the requested states at each step and their energy difference in eV.
- `molecule_ci_search_traj.xyz`: shows the geometrical trajectory of the optimization.
- `AUX_FILES/` folder containing the last step gradient calculations, geometries, output files and CDV and DGV if it was a GP or UBP calculation.

The energy profile of the states of interest can be plotted with `barista/personnel/fonsi.py -f energies.dat`. 

# Conical intersectioni characterization 


## Personnel
This is a collection of scripts to make short tasks. These modules serve as an interface between this program and electronic structure codes. [See the documentation](https://github.com/matarabadanf/BARISTA/blob/main/docs/Personnel.md).
