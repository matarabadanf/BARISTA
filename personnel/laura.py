#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
import numpy as np

np.set_printoptions(linewidth=np.inf)

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Exctract the final geometry of an optimization run in ORCA.""",
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="Input file name"
)
parser.add_argument("-o", type=str, default=True, help="Output file name")
parser.add_argument(
    "-g", type=str, default="True", help="Generate an xyz file"
)


def laura(
    filename: str, default_name: str = True, generate_file: bool = False
) -> [float, pd.DataFrame]:
    # The directory content is listed and the gs geometry is searched.

    f = open(filename, "r")
    cont = f.readlines()

    for line in cont:
        if "Number of roots to be determined" in line:
            n_roots = int(line.strip().split()[-1])
        elif 'ORCA SCF GRADIENT CALCULATION' in line:
            engrad = True

    gs_energy = 0

    engrad = False

    if engrad == True:
        for line in cont:
            if 'Total Energy ' in line:
                gs_energy = float(line.strip().split()[-4])

    for line in cont:
         if "E(SCF)  =" in line:
             gs_energy = float(line.strip().split()[-2])
             print('ground state energy is', gs_energy)

    energies = np.zeros(n_roots + 1)
    energies[0] = gs_energy

    counter = 0 

    for line in cont:
        if (
            "STATE" in line
            and "TD-DFT/TDA EXCITED STATES (SINGLETS)" not in line
            and "EXCITED STATE GRADIENT DONE" not in line
            and "NATURAL" not in line
        ):
            index = int(line.replace(':',' ').strip().split()[1])
            energies[index] = float(line.strip().split()[3])

    energies[1:] += gs_energy

    print(energies)

    for i in range(len(cont)):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in cont[i]:
            starting_index = i + 2
        elif "CARTESIAN COORDINATES (A.U.)" in cont[i]:
            end_index = i - 2
        elif "FINAL SINGLE POINT ENERGY" in cont[i]:
            ener = cont[i].strip().split()[-1]

    print(len(cont[starting_index:end_index]))
    print(ener)
    for i in cont[starting_index:end_index]:
        print(i.strip())

    if generate_file == "True":
        if default_name:
            output_name = (
                filename.replace(".in", "").replace(".out", "") + ".xyz"
            )
        else:
            output_name = default_name
        out_file = open(output_name, "w")
        out_file.write(str(len(cont[starting_index:end_index])) + "\n")
        out_file.write(str(energies).replace(']', '').replace('[', '') + "\n")
        for i in cont[starting_index:end_index]:
            out_file.write(i.strip() + "\n")

    positions = np.array(
        [
            np.array(vector.strip().split()[1:]).astype(float)
            for vector in cont[starting_index:end_index]
        ]
    )
    labels = [
        np.array(vector.strip().split()[0])
        for vector in cont[starting_index:end_index]
    ]

    df = pd.DataFrame(
        {
            "Atom labels": labels,
            "x": positions[:, 0],
            "y": positions[:, 1],
            "z": positions[:, 2],
        }
    )

    return ener, df


if __name__ == "__main__":

    # display help message if there is no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    print(args.g)
    a = laura(args.f, args.o, args.g)
    print(a[1])
