# -*- coding: utf-8 -*-
import numpy as np
import ase.io
from ase.calculators.calculator import FileIOCalculator


class CICalculator(FileIOCalculator):
    # Tthe calculator needs a command that is what the program will execute. 
    # Therefore it is needed to prepare the inputs. 
    command = "bash ~/bin/run_orca5.0.3 engrad_0.in  && ~/bin/run_orca5.0.3  engrad_1.in "
    
    # ASE requires to define the implemented properties of a calculator 
    implemented_properties = ["forces", "energy"]

    def __init__(
        self,
        restart=None,
        ignore_bad_restart_file=None,
        label=None,
        atoms=None,
        command=command,
        mode=None,
        profile=None,
        n_roots=10,
        iroot=0,
        jroot=1,
        functional="CAM-B3LYP",
        basis="cc-pvdz",
        alpha=0.02,
        sigma=3.5,
        n_procs=1,
        geom=None
    ):
        # The class methods have to be inherited or they will die
        super().__init__(
            restart=None,
            ignore_bad_restart_file=None,
            label=None,
            atoms=None,
            command=command,
            profile=None,
        )
        
        self.iroot = iroot
        self.jroot = jroot
        self.directory = "."
        self.prefix = label
        self.n_roots = (n_roots,)
        self.functional = (functional,)
        self.basis = basis
        self.atoms = atoms
        self.alpha = alpha
        self.sigma = sigma
        self.n_procs = n_procs
        # this prepares a current geometry
        with open(self.label + ".xyz", "w") as fd:
            ase.io.write(fd, self.atoms, format="xyz")
        # this saves the original geometry
        with open(self.label + "_original.xyz", "w") as fd:
            ase.io.write(fd, self.atoms, format="xyz")


    def write_input(self, atoms=None, properties=None, system_changes=None):
        FileIOCalculator.write_input(
            self, atoms, properties, system_changes
        )  # because the docs demmand so
        # we will only update the .xyz file as the input for calculations will remain the same throughout the whole calculation
        with open(self.label + ".xyz", "w") as fd:
            ase.io.write(fd, self.atoms, format="xyz")

        # this prepares the orca inputs
        for index, root in enumerate([self.iroot, self.jroot]):
            with open("engrad_%i.in" % index, "w") as engrad_file:
                engrad_file.write(
                    "! engrad {} {} \n! NoAutostart\n\n\n".format(self.functional[0], self.basis)
                )
                engrad_file.write("")
                engrad_file.write(
                    "%%tddft nroots %s iroot %s tda TRUE end\n"
                    % (self.n_roots[0], root)
                )
                engrad_file.write("* xyz 0 1\n")
                with open(self.label + ".xyz", "r") as xyzfile:
                    cont = xyzfile.readlines()
                for line in cont[2:]:
                    engrad_file.write(line)
                engrad_file.write("*\n\n")
                if self.n_procs != 1:
                    engrad_file.write("%%pal NPROCS %i END" % self.n_procs)


    def read_results(self):
        # read results from the ORCA calculation
        for index in range(0, 2):
            with open("engrad_%i.in.out" % index, "r") as enfile:
                cont = enfile.readlines()
                energies = np.array([0.0])

                # parse the energies and locate the gradients
                for jndex, line in enumerate(cont):
                    if "FINAL SINGLE POINT ENERGY" in line:
                        energies[0] = float(line.strip().split()[-1])
                    if "CARTESIAN GRADIENT" in line:
                        start_gradients = jndex + 3
                    if "Difference to translation invariance:" in line:
                        end_gradients = jndex - 1

                grads = np.zeros([end_gradients - start_gradients, 3])

                # parse the gradients
                for kndex, gradient in enumerate(
                    cont[start_gradients:end_gradients]
                ):
                    gradient = gradient.strip().split()
                    grads[kndex] = [
                        float(gradient[-3]),
                        float(gradient[-2]),
                        float(gradient[-1]),
                    ]

                # save the values in a new file
                with open("engrad_%i_energy.dat" % index, "w") as enerfile:
                    for energy in energies:
                        enerfile.write(
                            str(energy).replace("[", "").replace("]", "") + "\n"
                        )
                with open("engrad_%i_gradient.dat" % index, "w") as grd:
                    for grad in grads:
                        grd.write(
                            str(grad).replace("[", "").replace("]", "") + "\n"
                        )
        with open('energies.dat', 'a') as endat:
            en1 = np.loadtxt("engrad_0_energy.dat")
            en2 = np.loadtxt("engrad_1_energy.dat")
            endat.write('%.6f %.6f %.6f\n' % (en1, en2, (en2-en1)*27.211))

        self.penalty_results()
   
 
    def penalty_results(self):
        # management of the energies
        en1 = np.loadtxt("engrad_0_energy.dat")
        en2 = np.loadtxt("engrad_1_energy.dat")

        en = (en2 - en1) ** 2 / (en2 - en1 + self.alpha)
        en = self.sigma * en + (en2 + en1) / 2.0

        # management of the gradients
        grad1 = np.loadtxt("engrad_0_gradient.dat")
        grad2 = np.loadtxt("engrad_1_gradient.dat")

        # penalty algorithm
        kk = 2 * self.alpha * (en2 - en1)
        kk += (en2 - en1) ** 2
        kk /= (en2 - en1 + self.alpha) ** 2
        kk *= self.sigma

        grad_pen = (grad2 - grad1) * kk
        grad_sa = (grad2 + grad1) / 2.0
        grad = grad_pen + grad_sa
        grad = np.reshape(grad, self.atoms.positions.shape)

        # link results to class variable
        self.results["energy"] = en * ase.units.Hartree
        self.results["forces"] = grad
        self.results["forces"] *= -ase.units.Hartree / ase.units.Bohr
