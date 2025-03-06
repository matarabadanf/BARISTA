#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os 

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Takes a brewer optimizaiton energies.dat and plots the energy at each step.""",
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="energies.dat file"
)

parser.add_argument("-o", type=str, default=None, help="Output image filename")

parser.add_argument(
    "-en",
    type=float,
    default=0,
    help="Reference energy of the ground state optimized energy",
)

parser.add_argument("-u", type=str, default="eV", help="Energy units")

parser.add_argument(
    "--interactive",
    action='store_true',
    required=False,
    help="Open interactive plot.",
)
# display help message if there is no arguments

class Fonsi:

    def __init__(self,
            reference_energy:float, 
            energies_filename:str, 
            units:str='eV'
        ):
        self._reference_energy = reference_energy
        self._energies_filename = energies_filename
        self._units = units 

        self._initialize()

    def _initialize(self):
        self._load_energies()


    def _load_energies(self):
        self._energies_array = np.loadtxt(self._energies_filename) - self._reference_energy
        if self._units == 'eV':
            self._energies_array *= 27.2114 

    def _preplot(self):
        x = np.arange(1, len(self._energies_array.T[0]) + 1, 1)

        plt.plot(x, self._energies_array[:,0])
        plt.plot(x, self._energies_array[:,1])

        plt.title('CI Convergence')
        plt.xlabel("Step")
        plt.ylabel(f"Energy / {self._units}")

        return plt 

    def save_fig(self, name:str=''):
        self._preplot()
        if name == '':
            name = 'CI_convergence.jpg'
        if '.' not in name:
            name = name + '.jpg'

        plt = self._preplot()
        plt.savefig(name, dpi=600)

    def plot(self):
        plt = self._preplot()
        plt.show()

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

    a = Fonsi(args.en, args.f, args.u)
    print(a._energies_array)
    if args.o is not None:
        a.save_fig(args.o)
    if args.interactive:
        a.plot()
