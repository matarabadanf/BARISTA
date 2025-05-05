#!/usr/bin/env python3
import matplotlib.pyplot as plt 
import numpy as np 
import argparse
import sys
from functools import cached_property
from scipy.interpolate import make_interp_spline

parser = argparse.ArgumentParser(
    description="""Takes the ground state geometry of a molecule and compares it to optimized TDDFT geometries run in ORCA.\
        \nIt compares the energy, RMSD and final root, plotting the results in an image. """,
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-f",
    type=str,
    required=True,
    default=42,
    help="Trajectory .xyz file",
)
parser.add_argument(
    "-en",
    type=float,
    default=0,
    help='Total absolute energy in Hartree at Franc-Condon',
)

parser.add_argument(
    "-u",
    type=str,
    default='eV',
    help='Energy units for the plot (default="eV")',
)

parser.add_argument(
    "-o",
    type=str,
    default='default.jpg',
    help='Output image name.',
)

parser.add_argument(
    '--landscape', 
    action='store_true', 
    help='Use landscape layout'
)

parser.add_argument(
    '--spline', 
    action='store_true', 
    help='Use spline'
)

parser.add_argument(
    '--barrier', 
    action='store_true', 
    help='Plot absolute barrier'
)

class NebExtractor:

    def __init__(
            self,
            traj_file:str, 
            reference_energy:float=0,
            units:str = 'eV'
        ):
        
        self.traj_file = traj_file
        self.reference_energy = reference_energy 
        self.units = units 
        
    @cached_property 
    def file_content(self) -> list:
        with open(self.traj_file, 'r') as f:
            file_content = f.readlines()

        return file_content 

    @cached_property
    def energy_array(self) -> np.ndarray:
        energies = []

        for line in self.file_content:
            if 'E' in line:
                energies.append(float(line.strip().split()[-1]))

        if self.units == 'eV': 
            converter = 27.2114 
        else:
            converter = 1.0

        return (np.array(energies) - self.reference_energy) * converter

    @cached_property
    def forward_barrier(self) -> float:
        return abs(max(self.energy_array) - self.energy_array[0])

    @cached_property
    def reverse_barrier(self) -> float:
        return abs(max(self.energy_array) - self.energy_array[-1])

    @cached_property
    def absolute_barrier(self) -> float:
        return abs(max(self.energy_array) - min(self.energy_array))

    def _preplot(
            self, 
            spline:bool=False,
            annotate_barrier:bool=True, 
            title:str='None', 
            landscape:bool=False
        ):

        x = np.arange(0,len(self.energy_array))
        y = self.energy_array

        xi_yi_max = max(zip(x, y), key=lambda pair: pair[1])
        print(f"Max energy: {xi_yi_max[1]} at index {xi_yi_max[0]}")

        xi_yi_min = min(zip(x, y), key=lambda pair: pair[1])
        print(f"Min energy: {xi_yi_min[1]} at index {xi_yi_min[0]}")

        if landscape:
            plt.figure(figsize=(12, 6))
            ax = plt.gca()
            for spine in ax.spines.values():
                spine.set_visible(False)

        X_Y_Spline = make_interp_spline(x, y)
        X_ = np.linspace(x.min(), x.max(), 500)
        Y_ = X_Y_Spline(X_)
        if spline == True:
            plt.plot(X_, Y_)
        else:
            plt.plot(x, y)

        if annotate_barrier:
            ymin = min(y)
            ymax = max(y)
            barrier = ymax - ymin
            idx_ymin = xi_yi_min[0]
            idx_ymax = xi_yi_max[0]

            # Draw the double arrow using annotate twice
            plt.annotate(
                '', xy=(idx_ymax, ymax), xytext=(idx_ymax, ymin),
                arrowprops=dict(arrowstyle='<->', color='red', lw=1.5)
            )

            # Add energy barrier label
            plt.text(
                idx_ymax + 0.1, (ymax + ymin) / 2,
                f'{barrier:.2f} eV',
                color='red', va='center', fontsize=10
            )
            
        if title != 'None':
            plt.title(title)

        plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False
        )

    def savefig(self, name:str='default.jpg', dpi:int= 600, landscape:bool=False, spline=False, annotate_barrier:bool=False):
        self._preplot(landscape=landscape, spline=spline, annotate_barrier=annotate_barrier)

        plt.savefig(name, dpi=dpi)

        plt.show()


if __name__ == '__main__':
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    s = NebExtractor(args.f, args.en, args.u)

    if args.o != 'default.jpg':
        s.savefig(args.o, landscape=args.landscape, spline=args.spline, annotate_barrier=args.barrier)
        

    #s = NebExtractor('asdf.xyz', reference_energy=-533.24663318)
    #print(s.energy_array)

    #print(f'Forward barrier: {s.forward_barrier:6.4f}')
    #print(f'Reverse barrier: {s.reverse_barrier:6.4f}')
    #print(f'Absolute barrier: {s.absolute_barrier:6.4f}')

    #s._preplot()
   

