#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Exctract the final geometry of an optimization run in ORCA.""",
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="FC xyz file with energy as comment"
)

parser.add_argument(
    "-b", type=str, required=True, default=None, help="Last point in the branch"
)

parser.add_argument("-o", type=str, default=True, help="Output file name")

parser.add_argument(
    "-g", type=str, default="True", help="Generate an xyz file"
)



def aurora(
    fc_file: str, 
    last_branch_file: str,
) -> None:
    # The directory content is listed and the gs geometry is searched.

    with open(fc_file, 'r') as f:
        cont = f.readlines()
    
    reference_energy = float(cont[1].strip().split()[0])

    fc_energies = (np.array([float(energy) for energy in cont[1].strip().split()])-reference_energy)*27.211

    print(fc_energies)
    print(reference_energy)

    split_name = last_branch_file.replace('.xyz', '').replace('ci', '-1_-1_-1').split('_')

    for i, item in enumerate(split_name):
        try:
            int(item)
            index = i
            break
        except:
            print(item)

    content_to_weave = split_name[index:]

    general_name = split_name[:index]
    print(general_name)

    print(content_to_weave)

    content_to_weave = [int(i) for i in content_to_weave if i != '']

    identifiers = np.array(content_to_weave).reshape([-1,3])

    print(identifiers)

    a, b, c = identifiers[0]

    string_0 = '_'.join(general_name)+'_%i_%i_%i.xyz' % (a,b,c)

    filenames = [string_0]

    for identifier in identifiers[1:]:
        print(identifier)
        if max(identifier) == -1:
            filenames.append(filenames[-1].replace('.xyz', '') + '_ci.xyz')
        else:
            a, b, c = identifier
            filenames.append(filenames[-1].replace('.xyz', '') + '_%i_%i_%i.xyz' % (a,b,c))
    
    pes_energies = [fc_energies]

    print(filenames)

    for name in filenames:
        print(name)
        with open(name, 'r') as f:
            cont = f.readlines()
        pes_energies.append((np.array([float(energy) for energy in cont[1].strip().split()])-reference_energy)*27.2114)
    print(pes_energies)

    pes_energies = np.array(pes_energies).T

    print(pes_energies)

    x = np.arange(len(pes_energies[0]))

    for i, energy in enumerate(pes_energies):
        plt.plot(x, energy, alpha=.5, linestyle='--')
    plt.show()

    plt.savefig(filenames[-1].replace('xyz', 'jpg'))

if __name__ == "__main__":

    # display help message if there is no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    print(args.g)
    a = aurora(args.f, args.b)
