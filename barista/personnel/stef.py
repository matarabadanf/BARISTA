#!/usr/bin/env python3
import matplotlib.pyplot as plt 
import numpy as np 
import argparse
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
    def file_content(self):
        with open(self.traj_file, 'r') as f:
            file_content = f.readlines()

        return file_content 

    @cached_property
    def energy_array(self):
        energies = []

        for line in self.file_content:
            if 'E' in line:
                energies.append(float(line.strip().split()[-1]))

        if self.units == 'eV': 
            converter = 27.2114 

        return (np.array(energies) - self.reference_energy) * converter

    @cached_property
    def forward_barrier(self):
        return abs(max(self.energy_array) - self.energy_array[0])

    @cached_property
    def reverse_barrier(self):
        return abs(max(self.energy_array) - self.energy_array[-1])

    @cached_property
    def absolute_barrier(self):
        return abs(max(self.energy_array) - min(self.energy_array))

    def _preplot(self, spline:bool=False):
        x = np.arange(0,len(self.energy_array))
        y = self.energy_array
        X_Y_Spline = make_interp_spline(x, y)
        X_ = np.linspace(x.min(), x.max(), 50)
        Y_ = X_Y_Spline(X_)
        if spline == True:
            plt.plot(X_, Y_)
        else:
            plt.plot(x, y)
        plt.show()


if __name__ == '__main__':
    pass 
    

    s = NebExtractor('asdf.xyz', reference_energy=-533.24663318)
    print(s.energy_array)

    print(f'Forward barrier: {s.forward_barrier:6.4f}')
    print(f'Reverse barrier: {s.reverse_barrier:6.4f}')
    print(f'Absolute barrier: {s.absolute_barrier:6.4f}')

    s._preplot()
   

