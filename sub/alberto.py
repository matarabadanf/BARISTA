import argparse
import os
import numpy as np 
import re
import matplotlib.pyplot as plt 
import sys

# Parser is used to input via terminal the required arguments for Emma
parser=argparse.ArgumentParser(
    description='''Takes the ground state geometry of a molecule and compares it to optimized TDDFT geometries run in ORCA.\
        \nIt compares the energy, RMSD and final root, plotting the results in an image. ''',
    epilog="""It should run appropriately with ORCA 5.0.X""")

parser.add_argument('-f', type=str, required=True, default=42, help='Optimized ground state .xyz file')
parser.add_argument('-o', type=str, default=42, help='Output image filename')

# display help message if there is no arguments

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)


args=parser.parse_args()

# get the file and content

file = open(args.f, 'r')
cont = file.readlines()

# grep the lines of interest 

for line in cont:
    if 'STATE' in line and 'TD-DFT/TDA EXCITED STATES (SINGLETS)' not in line and 'EXCITED STATE GRADIENT DONE' not in line:
        print(line)

