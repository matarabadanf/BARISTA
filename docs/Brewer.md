# Brewer
BREWER is a small homebrew inteface with electronic structure programs (currently only ORCA) based in the python module `ase.io`. Via `ase.io` it is possible to interface electronic structure calculations with some highly customizable algorithms in order to perform complex tasks such as MECI calculations, NEBs... 


## Implemented Brewer functionality 
Currently supported calcuations are:
- MECI search with TDDFT in ORCA.

### MECI calculations

MECI calculations can be performed with BREWER using the following algorithms: 
 - Penalty algorithm.
 - Updated branching plane (no NACME).
 - Updated branching plane (Estimated NACME).

The first is implemented as in [J. Phys. Chem. B 2008, 112, 2, 405–413](https://pubs.acs.org/doi/10.1021/jp0761618). The second is implemented as in [J. Chem. Theory Comput. 2010, 6, 1538–1545](https://pubs.acs.org/doi/abs/10.1021/ct1000268).

#### Brief theoretical background of the MECI seatch methods



## Configuring Brewer
Brewer needs to be configured: it needs the path to a suitable ORCA (currently only orca, could be expanded) launcher. To set up the launcher scripts, use `setup_launchers.py`. The `setup_launchers.py` script can be run as such, or as:
```Bash
setup_launchers.py < configuration.dat
```

Where the configuration file should look like:
```
~/bin/run_orca.sh ! orca launcher
```

**IMPORTANT**: ORCA launchers should be able to copy `.xyz` and `.gbw` files to the `$WorkDir`. We include an example launcher for ORCA in the `data/` folder. This is in order to avoid recoputing the SCF optimization in both states, eventually saving some computation time. 

Due to the way we handled input generation, it is not necessary that the ORCA launcher manages the parrallelization `\%pal n_cores` section in orca inputs. 

## Using Brewer
Brewer requires a `brewer.TEMPLATE` input file. The `barista/brewer/brewer_setup.py ` script generates a directory with the `brewer.TEMPLATE` input file and copies the required files for the calculation with:
```Bash
brewer_setup.py -f geometry_to_optimize.xyz 
```
Along the files that are in the directory there is: 
- `opt.py`: Main function, manages reading the input and performs the optimization
- `flavs.py`: Contains the calculator called by ASE. ASE requires the use of calculators to obtain properties such as forces and energies.
- `run_opt.sh`: Main script for the queue system.

To run a calculation with BREWER in a SLURM queue system, it can be done as:
```Bash
sbatch -args run_opt.sh
```

Consider that the calculation will call an electronic structure program in the following way:

```Bash
/path/to/launcher/laucher.sh input.in
```

Launchers should have been already configured with `setup_launchers.py`. 

After the calculation has finished, a folder with results will be generated (hierarchy in progress). There results can be extracted. 


### MECI calculations
To perform meci calculations, a `brewer.TEMPLATE` file containing the following fields would suffice:
```
name = molecule.xyz
mode = ci
```
Where the `molecule.xyz` is the starting geometry with and `mode` is the type of job. There are some optional arguments for the CI optimization:
#### Keywords 
| Keyword    | Information|
|------------|------------|
| `profile`  | MECI search algoritm. Options are `penalty`, `gp` (Gradient projection) and `ubp` (Updated Branching plane).|
| `iroot`    | Lower state potential energy surface. Default is $0$.           |
| `jroot`    | Higher state potential energy surface. Default is `iroot`$+1$.           |
| `n_roots` | Number of excited states calculated by the electronic structure program. Default is $10$.       |
| `alpha`    | Alpha value for penalty method (see reference). Default is $3.5$.           |
| `sigma`    | Sigma value for penalty method (see reference). Default is $0.02$.          |

The default `brewer.TEMPLATE` configuration is:

```
name = molecule.xyz
mode = ci
profile = ubp 
iroot = 0
jroot = 1
n_roots = 10
```

#### Results
Once a MECI calculation has finished, The resulting files are:
- `optimization.log`: contains the optimizer steps, energies and max_force at each step. 
- `energies.dat`: shows the energy of the requested states at each step and their energy difference in eV.
- `molecule_ci_search_traj.xyz`: shows the geometrical trajectory of the optimization.
- `AUX_FILES/` folder containing the last step gradient calculations, geometries, output files and CDV and DGV if it was a GP or UBP calculation.  
