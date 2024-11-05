import argparse
import os
import numpy as np 
import re
import matplotlib.pyplot as plt 
import sys

# Parser is used to input via terminal the required arguments for Emma
parser=argparse.ArgumentParser(description='''Takes the ground state geometry of a molecule and compares it to optimized TDDFT geometries run in ORCA.\ \nIt compares the energy, RMSD and final root, plotting the results in an image. ''', epilog="""It should run appropriately with ORCA 5.0.X""")

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
counter = 0 
n_os = []

for line in cont:
    if 'TD-DFT/TDA EXCITED STATES (SINGLETS)' in line:
        counter+=1
        print('aaaa')
        if counter >= 2:
            break
    elif 'STATE' in line and 'TD-DFT/TDA EXCITED STATES (SINGLETS)' not in line and 'EXCITED STATE GRADIENT DONE' not in line:
        print(line)
        n_os.append(int(line.strip().split()[1].replace(':', '')))


total_list = [[] for i in range(0, max(n_os))]
cis_energies = []

for line in cont:
    if 'STATE' in line and 'TD-DFT/TDA EXCITED STATES (SINGLETS)' not in line and 'EXCITED STATE GRADIENT DONE' not in line:
        index = int(line.strip().split()[1].replace(':', '')) - 1
        energy = float(line.strip().split()[3].replace(':', ''))
        total_list[index].append(energy)
    elif 'E(SCF)' in line:
        cis_energies.append(float(line.strip().split()[2]))

# print(total_list)

cis_array = np.array(cis_energies)
total_arrays = [np.array(l) for l in total_list]

x = np.arange(1, len(total_list[0])+1, 1)

print(cis_array)

plt.plot(x, cis_array, label='Ground state')

# obtain the state of the actual root
curr_energy = []
for line in cont:
    if 'DE(CIS) =' in line:
        print(line)
        curr_energy.append(float(line.strip().split()[-4].replace(')', '')))

print(total_arrays)

print(curr_energy)

for i in range(len(total_arrays)):
    plt.plot(x, cis_array + total_arrays[i], label='Root %i'%(i+1))
plt.scatter(x, cis_array + curr_energy, label='Active root', marker='x', c='rebeccapurple')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('Step')
plt.ylabel('Energy / Hartree')
plt.savefig('asdf', dpi=600, bbox_inches='tight')
plt.show()
