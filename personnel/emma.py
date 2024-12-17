#!/usr/bin/env python3
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

parser.add_argument('-fc', type=str, required=True, default=42, help='Optimized ground state .xyz file')
parser.add_argument('-ex', type=str, required=True, nargs='+', default=42, help='Excited state optimization .out "filenames" (it is important that the name ends with *_N.in.out)')
parser.add_argument('-O', type=str, default=42, help='Output file extension (default=".in.out")')
parser.add_argument('-o', type=str, default=42, help='Output image filename')
parser.add_argument('-r', type=str, default=42, help='Generate a report table in output file')
parser.add_argument('-md', type=str, default=True, help='generate report in MD format (default is True)')

# display help message if there is no arguments

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)


args=parser.parse_args()



# The directory content is listed and the gs geometry is searched. 
dir_cont = os.listdir(".")

if args.fc not in dir_cont:
    print('No ground state geometry file %s was found in this directory' % args.fc)
    exit()


# Check if there is a user defined output extension
if isinstance(args.O, str):
    excited_strings = [i.replace('N%s' % args.O, '') for i in args.ex]
    output_extension = args.O
else:
    excited_strings = [i.replace('N.in.out', '') for i in args.ex]
    output_extension = '.in.out'


# Search for the excited optimization files of interest
excited_files = []
starting = []
follow = []

for string in excited_strings:
    for directory_file in dir_cont:
        if string == directory_file[0:len(string)] and output_extension in directory_file:
            try:
                int(directory_file[len(string)])
                excited_files.append(directory_file)
            except:
                pass
if excited_files == []:
    for i in excited_strings:
        print('No output files with the format %sN.in.out were found, please check' % ''.join(i))
    exit

# Ground state file is loaded

gs_file = open(args.fc, 'r')
gs_contents = gs_file.readlines()
gs_energy = float(gs_contents[1].strip().split()[-1])
gs_coordinates = []

for atom in gs_contents[2:]: 
    to_append = [float(a) for a in atom.strip().split()[1:]]
    gs_coordinates.append(np.array(to_append))

gs_file.close()

# Parsing the orca Output_files 

excited_optimization_coordinates = []
excited_optimization_states = []
excited_optimization_final_roots = []
excited_optimization_final_energies = []
non_converged = []

for file in excited_files:
    ex_file = open(file, 'r')
    if re.findall('did not converge', ex_file.read()):
        non_converged.append(file)
        ex_file.close()
        continue
    
    ex_file = open(file, 'r')
    if re.findall('ORCA TERMINATED NORMALLY', ex_file.read()):
        ex_file.seek(0)
        if 'follow' in file:
            follow.append(True)
        else:
            follow.append(False)
        starting.append(file.replace('.in.out', '').split('_')[-1])
        ex_cont = ex_file.readlines()
        energies = []
        geom_start = []
        geom_end = []
        iroots = []
        coordinates = []
        print('Reading data from file: ', file)
        for line in ex_cont:
            if 'FINAL SINGLE POINT ENERGY' in line:
                energies.append(float(line.strip().split()[-1]))
            elif '(Root' in line:
                iroots.append(int(line.replace(')', '').strip().split()[-1]))
            elif 'FINAL ENERGY EVALUATION AT THE STATIONARY POINT' in line:
                geom_start = ex_cont.index(line) + 6
                geom_end = geom_start + len(gs_contents) - 2
                for atom in [l.strip().split() for l in ex_cont[geom_start:geom_end]]:
                    coordinates.append(np.array([float(atom[1]), float(atom[2]), float(atom[3])])) 
        
        if iroots == [] or energies == []:
            print(file, ' did not contain extractable information')

        else:
            excited_optimization_states.append(file)   
            excited_optimization_final_roots.append(iroots[-1])
            excited_optimization_final_energies.append(energies[-1])
            excited_optimization_coordinates.append(coordinates)
 
# Data has been parsed, now RMSD (or something alike is computed)

RMSDs = []

for molecule in range(0, len(excited_optimization_states)):
    RMSD = 0
    natoms = len(excited_optimization_coordinates[molecule])
    for atom in range(natoms):
        RMSD = np.sqrt(1/natoms) * np.sqrt(np.linalg.norm(excited_optimization_coordinates[molecule][atom] - gs_coordinates[atom])**2)
    RMSDs.append(RMSD)

# If report is requested, perform it 
longest = 0
for name in excited_optimization_states:
    if len(name) > longest:
        longest = len(name)

whitespaces = ''.join(['0' for i in range(longest)])

for name in range(len(excited_optimization_states)):
    excited_optimization_states[name] += ''.join([' ' for i in range(longest-len(excited_optimization_states[name]))])

final_array = []

for i in range(len(RMSDs)):
    final_array.append(np.array([int(starting[i]), str(follow[i]), float(excited_optimization_final_energies[i]), float(excited_optimization_final_energies[i] - gs_energy), float((excited_optimization_final_energies[i] - gs_energy)*27.211), int(excited_optimization_final_roots[i]), float(RMSDs[i])]))

final_array = sorted(final_array,key=lambda x: x[3])

if isinstance(args.r, str):
    file = open(args.r, 'w')
    if args.md == True:
        file.write('|STARTING ROOT | FOLLOWIROOT | FINAL ENERGY | $\Delta E$ / Hartree |  $\Delta E$ / eV | Final State | RMSD|\n')
        file.write('|--------------|-------------|--------------|----------------------|------------------|-------------|-----|\n')
    else:
        file.write('STARTING ROOT, FOLLOWIROOT, FINAL ENERGY, DELTA_E (Hartree), Delta_E (eV), Final State, RMDS\n')
    for datat in final_array:
        data = [int(datat[0]), str(datat[1]), float(datat[2]), float(datat[3]), float(datat[4]), int(datat[5]), float(datat[6])]
        if args.md == True:
            file.write(f'|${data[0]:<3}$|{data[1]:<6}|${data[2]:^20.6f}$|${data[3]:^20.6f}$|${data[4]:^20.2f}$|${data[5]:^20}$|${data[6]:^20.4f}$|\n'.replace(' ', ''))
        else:
            file.write(f'{data[0]:<3}{data[1]:<6}{data[2]:^20.6f}{data[3]:^20.6f}{data[4]:^20.2f}{data[5]:^20}{data[6]:^20.4f}\n')
    for fi in non_converged:
        file.write('%s did not converge in geometry optimization\n' % fi)

# Obtain and save the Figure:

number_of_final_roots = len(set(excited_optimization_final_roots))

cmap = plt.get_cmap("tab10")

for geom in range(len(RMSDs)):
    plt.scatter(RMSDs[geom], (excited_optimization_final_energies[geom]-gs_energy)*27.214, color=cmap(excited_optimization_final_roots[geom]), label='Final root: %i' % excited_optimization_final_roots[geom])

   

plt.ylabel('$\Delta$E')
plt.xlabel('RMSD')
plt.legend()

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1, 1))


if isinstance(args.o, str):
    plt.savefig(args.o, dpi=600, bbox_inches='tight')
else:
    plt.show()

