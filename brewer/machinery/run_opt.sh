#!/bin/bash

#SBATCH -n 4 -c 1 -N 1 
#SBATCH -t 3-00:00:00
#SBATCH --qos=gp_resa
#SBATCH -A uam77

module load python/3.10.13

python3 opt.py

