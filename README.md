# BARISTA 

A collection of scripts and programs for computational chemistry packages.
## Brewer
This is a coffee machine in development that would be able to perform tasks that require some interfacing between programs. It is based in the python module ase.io. The idea is that via ase.io it is possible to interface electronic structure calculations with some highly customizable algorithms in order to perform complex tasks such as search of MECIs, perform NEBs (in development) and many more. As of today it should be possible to:
 - Obtain MECIs using the penalty and update branch methods. In theory, but the coffe machine was not used in a while and needs to be tested again. 

### Configuring Brewer
Brewer needs to be configured: it needs the path to a suitable ORCA launcher. #todo

### Using Brewer
Brewer requires a `brewer.TEMPLATE` file that contains at least two fields:
```
name = molecule.xyz
mode = ci
```
Where the molecule.xyz is the molecule to be worked with and the mode is the type of optimization (for now only CI). There are some optional arguments for the CI optimization:
```
profile = penalty
iroot = 0
jroot = 1
n_roots = 10
alpha  = 0.02 
sigma = 3.5
``` 
These are default optional values. The profile selects the CI optimization method. Penalty is based on the penalty algorithm [citation needed]. UPB is based on the updated branching plane algorithm [citation needed] (Not yet implemented).

Once a brewer.TEMPLATE file is generated, by running (or adding to a queue) the run_opt.sh script the optimization takes place. The resulting files are:
- `optimization.log`: contains the optimizer steps, energies and max_force
- `energies.dat`: shows the energy of the requested states at each step and their energy difference in eV
- `molecule_ci_search.traj`: shows the geometrical trajectory of the optimization 

### How Brewer works 
Brewer is based in the ASE.io python package. It utilizes its classes and methods to perform the required calculations. For this, there are two python files:
- `opt.py`: Main function, manages reading the input and performs the optimization
- `flavs.py`: Contains the calculator called by ASE. ASE requires the use of calculators to obtain properties such as forces and energies. A calculator should be able to: 1. write an input file, 2. perform the calculation of the input file with an external command (such as a call to ORCA) and 3. parse the results of the file and obtain the forces to feed to the optimization algorithm. 

The difference in `modes` can be set up in the calculator, as the forces to optimized are defined differently in each CI algorithm. 

With the calculator defined, the optimization algorithm workflow is the general one: 1. Calculate forces, 2. Displace geometry, repeat. 

## Personel
These are small scripts to make short tasks more bearable (like coffee). Each person is specialized in a certain flavour:
 - Emma.py: In the context of excited state optimizations, to distinguish if different states fall to the same minima, the difference from the ground state structure can be seen by checking the RMSD of the resulting excited optimized structure. This script takes the ground state .xyz file and the list of optimized structures and plots the energy and RMSD of the excited state minima found. 
 - Rodrigo.py: Script to plot the UV absorption spectra from a TDDFT ORCA calculation with options such as cuttoff energies, peaks and units.  
 - Lorena: Script to shake the geometry of an xyz file.
