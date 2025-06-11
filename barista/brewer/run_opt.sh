#!/bin/bash

#SBATCH -n 4 -c 1 -N 1 
#SBATCH -t 3-00:00:00
 module load gcc/11.2.0
 module load openmpi/gnu/4.1.1_gnu_11.2.0
 module load orca/5.0.3
module load orca/5.0.3
module load hdf5 && module load python/3.10.0 # MN module requirement for ase
# module load python/3.10.13 # CCC module requirement for ase

python3 opt.py

