#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os 

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Takes a excited state optimization in ORCA and plots the energy of each state at each step of the opitmizaiton and the actual root.""",
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="ORCA optimization file"
)
parser.add_argument("-o", type=str, default=True, help="Output image filename")
parser.add_argument(
    "-en",
    type=float,
    default=0,
    help="Reference energy of the ground state optimized energy",
)
parser.add_argument("-u", type=str, default="eV", help="Energy units")

# display help message if there is no arguments

class Fonsi:

    def __init__(self, fc_energy, energies_filename):
        self.fc_energy = fc_energy
        self.energies_filename = energies_filename


    def _load_energies:
        try:
            self.energies_array = np.loadtxt(self.energies_filename) - self.reference_energy
        except:
            pass
    def _test(self):
        pass 

def fonsi(
    filename: str,
    reference_energy: float = 0,
    output_image: str = True,
    units: str = "eV",
) -> np.array:
    # get the file and content

    energies = np.loadtxt(filename) - reference_energy

    if units == 'eV':
        energies *= 27.2114

    x = np.arange(1, len(energies) + 1, 1)

    plt.plot(x, energies[:,0])
    plt.plot(x, energies[:,1])


    if '/' in filename:
        title = filename.strip().split('/')[-2]
    else:
        title = os.getcwd().strip().split('/')[-1]


    if output_image is True:
        output_image = title
       
        if '/' in filename:
            output_image = '/'.join(filename.strip().split('/')[0:-1])+'/'+filename.strip().split('/')[-2]
        else:
            output_image = os.getcwd().strip().split('/')[-1]

    plt.title(title)
    plt.xlabel("Step")
    plt.ylabel("Energy / %s" % units)
    print(output_image)
    plt.savefig(
        output_image.replace(".dat", ".jpg"), dpi=250, bbox_inches="tight"
    )
    
    plt.show()
    
    print(energies[:,1][-1])
    
    return None


if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    a = fonsi(args.f, args.en, args.o, args.u)
