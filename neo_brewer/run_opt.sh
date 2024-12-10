#!/bin/bash

#SBATCH -n 4 -c 1 -N 1 
#SBATCH -t 3-00:00:00

module load orca/5.0.3
module load hdf5 && module load python/3.12.1

python3 opt.py

