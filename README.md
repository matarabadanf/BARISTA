# BARISTA 

A collection of scripts and small programs to assist in PES exploration in excited states. BARISTA is divided in two main funcitonalities:
- **BREWER** is a connical intersection optimizer based in the gradient projection and updated branching plane methods that is modular and can be interfaced with electronic structure programs:
  - Optimizes Conical intersections with: Penalty algorithm, Projection method (no CDV) or Updated branching plane (Estimated CDV).
  - Characterizes the resulting conical intersections.
  - Plots and has some tools to visualize branching plane directions, minimum directions... Can modify geometries for further optimizations. 
- **Personnel** contains a collection of modules that can be executed as command-line scripts or imported for python scripting with various functionalities.  

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
- The run_opt.sh script should include all module imports required for the calculation. This must include a python version that has ASE installed. 
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
Input keywords:

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
| `slow_start`    | Add a weight on the first steps to rescale the forces in order to ensure a smooth start. Default is False.                    |
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
Conical intersections can be characterized with the scripts `javi.py` and `cris.py`. The difference lies in the extended functionality of Cris. After a MECI optimization, these modules are called and a characterization takes place using Fdez. Galvan et. al. algorithm ([J. Chem. Theory Comput. 2016, 12, 8, 3636–3653](https://pubs.acs.org/doi/10.1021/acs.jctc.6b00384)). Results are displayed as: 

```
Beta has been rotated 2 * pi/2. Beta =   3.2719 (187.4639 deg)

pitch     =   0.0949
asymmetry =   0.5321
tilt      =   0.9552
theta     =   0.0000 (  0.0000 deg)

With the rotation, the conditions 0 < theta < pi/2 and 0 < asymmetry are fulfiled: True

P and B values are:

  0.5955, Peaked
  1.0727, Single Path
```

Cris.py can be used indipendently to work in CLI using a molcas output file (even though characterization is performed by molcas) for plotting reasons, displacement of geometries or visualization of directions. 

# Personnel
This is a collection of scripts to make short tasks. These modules serve as an interface between this program and electronic structure codes. 

## Command Line Interfaces

* [alberto](docs/Personnel.md#alberto): Takes a excited state optimization in ORCA and plots the energy of each state at each step of the opitmizaiton and the actual root.
* [aurora](docs/Personnel.md#aurora): Exctract the potential energy surface path and energies of of a PES branch.
* [cris](docs/Personnel.md#cris): Performs conical intersectioni characterization and includes extra functionality. Generates plots, force files to visualize the x and y vectors. Generated displaced geometries to optimize after the conical intersection. 
* [dani](docs/Personnel.md#dani): Extracts frequencies from gaussian calculation.
* [emma](docs/Personnel.md#emma): Takes the ground state geometry of a molecule and compares it to optimized TDDFT geometries run in ORCA. It compares the energy, RMSD and final root, plotting the results in an image. 
* [fonsi](docs/Personnel.md#fonsi): Plots energy of states of interest in BREWER optimization. 
* [javi](docs/Personnel.md#javi): Performs MECI characterization in a minimal way to work in all systems with less requirements than Cris. 
* [laura](docs/Personnel.md#laura): Extractes geometry and energies of a single point or final point in an optimization and adds energies as comment line to xyz file. This was useful when working with multiple geometries to use dataframes.
* [lorena](docs/Personnel.md#lorena): Randommly shakes a molecule. 
* [maxime](docs/Personnel.md#maxime): Extracts a spehere in a molecular crystal with the complete involved molecules.
* [rodrigo](docs/Personnel.md#rodrigo): Spectra ploting, modification and various utilites (works better as a module)
* [stef](docs/Personnel.md#stef): Extracts energies from ORCA NEB calculation. 
