#!/usr/bin/env python3
import argparse
import random
import numpy as np
import sys

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Takes the an xyz file and shakes the geometry""",
    epilog="""It should run appropriately with any .xyz file""",
)

parser.add_argument("-f", type=str, default=42, help="xyz file")
parser.add_argument("-d", type=float, default=0.1, help="Maximum displacement")
parser.add_argument("-o", type=str, default=42, help="Output file")

# display help message if there is no arguments


def lorena(filename: str, output_file: str, displacement: float = 0.1) -> None:
    f = open(filename, "r")
    cont = f.readlines()

    symbols = []
    coordinates = []

    for atom in cont[2:]:
        symbols.append(atom.strip().split()[0])
        coordinates.append(
            np.array([float(coord) for coord in atom.strip().split()[1:]])
        )

    print(symbols)
    print(coordinates)

    d = displacement

    pos_neg = [1, -1]
    for i in range(len(symbols)):
        random_vector = np.array(
            [
                np.random.rand() * random.choice(pos_neg) * d,
                np.random.rand() * random.choice(pos_neg) * d,
                np.random.rand() * random.choice(pos_neg) * d,
            ]
        )
        print(random_vector)
        coordinates[i] += random_vector

    if args.o == 42:
        pass
    else:
        resultfile = open(output_file, "w")
        resultfile.write(cont[0])
        resultfile.write(cont[1])
        for i in range(len(coordinates)):
            data = [
                symbols[i],
                coordinates[i][0],
                coordinates[i][1],
                coordinates[i][2],
            ]
            resultfile.write(
                f"{data[0]:<3}{data[1]:^20.6f}{data[2]:^20.6f}{data[3]:^20.2f}\n"
            )


if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    a = lorena(args.f, args.o, args.d)
