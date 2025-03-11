# Personnel
All the scripts used are written as the instance of different Classes, allowing the possibilityu to interface said classes for more complex functionality. All scripts can be made executable and are found in `barista/personnel/`.  
## Geometry and PES scripts  
 - Emma.py: In the context of excited state optimizations, to distinguish if different states fall to the same minima, the difference from the ground state structure can be seen by checking the RMSD of the resulting excited optimized structure. This script takes the ground state .xyz file and the list of optimized structures and plots the energy and RMSD of the excited state minima found. Also now it interfaces with the jeremy class, allowing finer differenciation of conformations. 
 - Rodrigo.py: Script to plot the UV absorption spectra from a TDDFT ORCA calculation with options such as cuttoff energies, peaks and units.  
 - Lorena: Script to shake the geometry of an xyz file.
 - Laura.py: Extracts the las geometry of an orca optimization and sets the comment string to the root energies.
 - Jeremy.py: using a geometry, calculates the internal coordinates. Useful to detemine the internal coordinate displacement between geometries.
 - Dani.py: Extracts geomtry and frequencies from orca calculation.
 - Stef.py: extracts NEB energies from orca calculation.
 - Alberto.py: extracts excited optimization path of orca calculation.
 - Aurora.py: Plots the branch diagram of a PES. It will be nterfaced with Stef for whole path characterization. 

## CI scripts 
- fonsi.py: shows the path enmergies in a CI optimization
- javi.py: characterizes and plots a conical intersection and the branching plane. Characterization is made based in the P and B values as implemented in [J. Chem. Theory Comput. 2016, 12, 8, 3636–3653](https://pubs.acs.org/doi/10.1021/acs.jctc.6b00384)
