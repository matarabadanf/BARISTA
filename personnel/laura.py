#!/usr/bin/env python3
import argparse
import os
import numpy as np 
import matplotlib.pyplot as plt 
import sys

# Parser is used to input via terminal the required arguments for Emma
parser=argparse.ArgumentParser(
    description='''Exctract the final geometry of an optimization run in ORCA.''',
    epilog="""It should run appropriately with ORCA 5.0.X""")

parser.add_argument('-f', type=str, default=42, help='Input file name')
parser.add_argument('-o', type=str, default=True, help='Output file name')


def laura(filename: str, default_name:str = True) -> None:    
    # The directory content is listed and the gs geometry is searched. 
    
    f = open(filename, 'r')
    cont = f.readlines()
    
    for i in range(len(cont)):
        if 'CARTESIAN COORDINATES (ANGSTROEM)' in cont[i]:
            starting_index = i + 2
        elif 'CARTESIAN COORDINATES (A.U.)' in cont[i]:
            end_index = i -2 
        elif 'FINAL SINGLE POINT ENERGY' in cont[i]:
            ener = cont[i].strip().split()[-1]
    
    print(len(cont[starting_index:end_index]))
    print(ener)
    for i in cont[starting_index:end_index]:
        print(i.strip())
    
    if default_name == True:
        output_name = args.f.replace('.in', '').replace('.out', '') + '.xyz'
    else:
        output_name = default_name
    out_file = open(output_name, 'w')
    out_file.write(str(len(cont[starting_index:end_index])) + '\n')
    out_file.write(str(ener) + '\n')
    for i in cont[starting_index:end_index]:
        out_file.write(i.strip()+'\n')


if __name__ == '__main__':

    # display help message if there is no arguments
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args=parser.parse_args()
    
    laura(args.f, args.o)
