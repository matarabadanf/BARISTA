import os

import numpy as np
import ase.io.xyz
from ase.units import Hartree, Bohr
from ase.calculators.calculator import FileIOCalculator, ReadError


class penalty(FileIOCalculator):
    # requires an engrad0 and engrad1 files, where iroots must be iroots i and j therefore two calculations are performed
    implemented_properties = ['energy', 'forces']
    command = 'bash run_template.sh PREFIX 1 && bash run_template.sh PREFIX 2'

    def __init__(self, restart=None, ignore_bad_restart_file=False, label='geom', atoms=None, sigma=3.5, alpha=0.02, **kwargs):
        FileIOCalculator.__init__(
            self, restart, ignore_bad_restart_file, label, atoms, **kwargs)
        self.pcpot = None
        self.sigma = sigma
        self.alpha = alpha

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def write_input(self, atoms, properties=None, system_changes=None):
        FileIOCalculator.write_input(self, atoms, properties, system_changes)

        cmt = "energy"
        if "forces" in properties:
            cmt = "gradient"

        with open(self.label+".xyz", "w") as fd:
            ase.io.write(fd, self.atoms, comment=cmt, format="xyz")

    def read(self, label):
        FileIOCalculator.read(self, label)
        if not os.path.isfile("run_template.sh"):
            raise ReadError

        self.read_results()

    def read_results(self):
        with open('run_template.sh', 'r') as otempl:
            cont = otempl.readlines()
            for line in cont:
                if 'run_orca' in line:
                    orca_templatename = line.strip().split()[1][:-2]

        if os.path.isfile(orca_templatename + '1.engrad') and os.path.isfile(orca_templatename + '2.engrad'):
            for root in ["1", "2"]:
                with open(orca_templatename + root + ".engrad", "r") as enfile:
                    cont = enfile.readlines()
                    energies = []
                    grads = []
                    for i in range(len(cont)):
                        print(cont[i].strip())
                        if cont[i].strip() == '# The current total energy in Eh' and cont[i+1].strip() == '#':
                            ener_start = i+2
                        elif cont[i].strip() == '#' and cont[i+1].strip() == '# The current gradient in Eh/bohr':
                            ener_end = i-1
                        elif cont[i].strip() == '# The current gradient in Eh/bohr' and cont[i+1].strip() == '#':
                            grad_start = i+2
                        elif cont[i].strip() == '#' and cont[i+1].strip() == '# The atomic numbers and current coordinates in Bohr':
                            grad_end = i-1

                print('In the parsing ', ener_start,
                      ener_end, grad_start, grad_end)
                energies = np.array([float(j)
                                    for j in cont[ener_start:ener_end+1]])
                grads = np.array([float(j)
                                 for j in cont[grad_start:grad_end+1]])
                print(grads, energies)

                with open(self.label + ".energy" + root, "w") as enerfile:
                    for energy in energies:
                        enerfile.write(str(energy) + '\n')

                with open(self.label + ".grad" + root, "w") as gradfile:
                    for grad in grads:
                        gradfile.write(str(grad) + '\n')

        if not os.path.isfile(self.label + '.energy1'):
            raise ReadError
        if not os.path.isfile(self.label + '.energy2'):
            raise ReadError
        
        en1 = np.loadtxt(f'{self.label}'+".energy1")
        en2 = np.loadtxt(f'{self.label}'+".energy2")
        
        en = (en2-en1)**2/(en2-en1+self.alpha)
        en = self.sigma * en + (en2+en1)/2.
        
        self.results['energy'] = en * ase.units.Hartree
        
        if os.path.isfile(self.label+'.grad1') and os.path.isfile(self.label+'.grad2'):
            grad1 = np.loadtxt(f'{self.label}'+".grad1")
            grad2 = np.loadtxt(f'{self.label}'+".grad2")

            kk = 2*self.alpha*(en2-en1)
            kk += (en2-en1)**2
            kk /= (en2-en1+self.alpha)**2
            kk *= self.sigma

            grad_pen = (grad2-grad1)*kk
            grad_sa = (grad2+grad1)/2.
            grad = grad_pen + grad_sa

            grad = np.reshape(grad, self.atoms.positions.shape)
            self.results['forces'] = grad
            self.results['forces'] *= - Hartree / Bohr
            print(self.results['forces'])
