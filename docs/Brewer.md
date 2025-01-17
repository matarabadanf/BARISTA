# Brewer
BREWER is a small homebrew inteface with electronic structure programs (currently only ORCA) based in the python module `ase.io`. Via `ase.io` it is possible to interface electronic structure calculations with some highly customizable algorithms in order to perform complex tasks such as MECI calculations, NEBs (in development)... 


## Implemented Brewer functionality 

### MECI calculations

MECI calculations can be performed with BREWER using the following algorithms: 
 - Penalty algorithm
 - Updated branching plane (no NACME)
 - Updated branching plane (Estimated NACME)

The first is implemented as in [](). The second is implemented as in ![]().

#### Brief theoretical background of the MECI seatch methods



## Configuring Brewer
Brewer needs to be configured: it needs the path to a suitable ORCA (currently only orca, could be expanded) launcher. To set up the launcher scripts, use `setup_launchers.py`. Launchers should be able to copy `.xyz` and `.gbw` files to the `\$WorkDir`. We include an example launcher for ORCA in the `data/` folder. 

Due to the way we handled input generation, it is not necessary that the ORCA launcher manages the parrallelization `\%pal n_cores' section in orca inputs. 

## Using Brewer
Brewer requires a `brewer.TEMPLATE` input file. The `barista/brewer/setup_brewer_calculation.py` script generates the brewer.TEMPLATE input file and copies the required files for the calculation. These can include:
- `opt.py`: Main function, manages reading the input and performs the optimization

- `flavs.py`: Contains the calculator called by ASE. ASE requires the use of calculators to obtain properties such as forces and energies. A calculator should be able to: 1. write an input file, 2. perform the calculation of the input file with an external command (such as a call to ORCA) and 3. parse the results of the file and obtain the forces to feed to the optimization algorithm. 

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
| `n\_roots` | Number of excited states calculated by the electronic structure program. Default is $10$.       |
| `alpha`    | Alpha value for penalty method (see reference). Default is $3.5$.           |
| `sigma`    | Sigma value for penalty method (see reference). Default is $0.02$.          |

The default 'brewer.TEMPLATE' configuration is:

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
- `molecule_ci_search.traj`: shows the geometrical trajectory of the optimization.
- `AUX_FILES/` folder containing the last step gradient calculations, geometries, output files and CDV and DGV if it was a GP or UBP calculation.  
