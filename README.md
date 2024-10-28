# BARISTA 

A collection of scripts and programs for computational chemistry packages.
## Brewer
This is a coffee machine in development that would be able to perform tasks that require some interfacing between programs. It is based in the python module ase.io. The idea is that via ase.io it is possible to interface electronic structure calculations with some highly customizable algorithms in order to perform complex tasks such as search of MECIs, perform NEBs (in development) and many more. As of today it should be possible to:
 - Obtain MECIs using the penalty and update branch methods. In theory, but the coffe machine was not used in a while and needs to be tested again. 

## Personel
These are small scripts to make short tasks more bearable (like coffee). Each person is specialized in a certain flavour:
 - Emma.py: In the context of excited state optimizations, to distinguish if different states fall to the same minima, the difference from the ground state structure can be seen by checking the RMSD of the resulting excited optimized structure. This script takes the ground state .xyz file and the list of optimized structures and plots the energy and RMSD of the excited state minima found. 
 - Rodrigo.py: Script to plot the UV absorption spectra from a TDDFT ORCA calculation with options such as cuttoff energies, peaks and units.  
 - Lorena: Script to shake the geometry of an xyz file.
