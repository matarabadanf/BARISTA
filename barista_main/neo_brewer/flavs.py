# -*- coding: utf-8 -*-
import numpy as np
import ase.io
from ase.calculators.calculator import FileIOCalculator


class CICalculator(FileIOCalculator):
    # Tthe calculator needs a command that is what the program will execute. 
    # Therefore it is needed to prepare the inputs. 
    command = " ~/bin/run_orca.sh engrad_0.in; cp engrad_0.gbw tmp.gbw;  ~/bin/run_orca.sh engrad_1.in"
    
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
        profile='penalty',
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
        self.profile = profile 
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
        try:
            with open('tmp.gbw', 'r') as f:
                prev_orb = True
        except:
            prev_orb = False

        # this prepares the orca inputs
        for index, root in enumerate([self.iroot, self.jroot]):
            with open("engrad_%i.in" % index, "w") as engrad_file:
                engrad_file.write(
                    "! engrad {} {} \n! NoAutostart\n\n\n".format(self.functional[0], self.basis)
                )
                if prev_orb == True: 
                    engrad_file.write('%moinp "tmp.gbw" \n%scf guess moread end\n'
                )
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
#                if self.n_procs != 1:
#                    engrad_file.write("%%pal\n NPROCS %i\nEND\n" % self.n_procs)


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
        
        if self.profile == 'penalty':
            self.penalty_results()
        if self.profile == 'gp':
            self.gradient_projection_forces()
        if self.profile == 'ubp':
            self.ubp_forces()
 
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
    
    def gradient_projection_forces(self):
        en1 = np.loadtxt("engrad_0_energy.dat")
        en2 = np.loadtxt("engrad_1_energy.dat")

        # management of the gradients
        grad1 = np.loadtxt("engrad_0_gradient.dat")
        grad2 = np.loadtxt("engrad_1_gradient.dat")
        
        grad1 = grad1.reshape(-1)
        grad2 = grad2.reshape(-1)
        
        x = (grad1-grad2)
        x /= np.linalg.norm(x)
         
        g_diff = 2 * (en1 - en2) * x
        
        P = np.identity(len(grad1)) - np.outer(x,x)

        total_gradient = g_diff + P @ (grad1 + grad2)/2
        total_gradient = total_gradient.reshape([-1,3])
        
        self.results["energy"] = en2 - en1  
        self.results["forces"] = - total_gradient * (ase.units.Hartree / ase.units.Bohr)

    def ubp_forces(self):
        en1 = np.loadtxt("engrad_0_energy.dat")
        en2 = np.loadtxt("engrad_1_energy.dat")
    
        # management of the gradients
        grad1 = np.loadtxt("engrad_0_gradient.dat").reshape(-1)
        grad2 = np.loadtxt("engrad_1_gradient.dat").reshape(-1)
    
        try:
            x_minus_one = np.loadtxt('x_minus_one.dat')
            y_minus_one = np.loadtxt('y_minus_one.dat')
            first_iter = False
        except:
            first_iter = True


        # x is the unit vector parallel to the difference gradient vector
        x = (grad1-grad2).reshape(-1)
        x /= np.linalg.norm(x)
        
        # this is a section to reshape everything so it has the correct dimensions 
        if not first_iter:
            x_minus_one = x_minus_one.reshape(-1)
            y_minus_one = y_minus_one.reshape(-1)
        x = x.reshape(-1)

        # in the first iteration, the (normalized) mean energy gradient vector is chosen as y_0    
        if first_iter == True:
            y = (grad1 + grad2).reshape(-1)/2
            y /= np.linalg.norm(y)
        
        # in the rest of iterations, y is calculated from the previous steps as the formula states. 
        else:
            y = (np.dot(y_minus_one, x) * x_minus_one - np.dot(x_minus_one, x) * y_minus_one) / (np.dot(y_minus_one, x)**2 + np.dot(x_minus_one, x)**2)**0.5
            y /= np.linalg.norm(y)
        
        y = y.reshape(-1)

        # All vectors x, y, x_minus and y_minus should have dimension 1 
    
        g_diff = 2 * (en1 - en2) * x
    #     g_diff /= np.linalg.norm(g_diff)
        mean_grad = (grad1 + grad2)/2
    #     mean_grad /= np.linalg.norm(mean_grad)
    #    print('gdiff is:')
    #    print(g_diff)
    #    print('mean_grad is:')
    #    print(mean_grad)
    
        x = x.reshape(-1)
        y = y.reshape(-1)
    #    print(x)
    #    print(y)
        P = np.identity(len(grad1.reshape(-1))) - np.outer(x,x) - np.outer(y,y)
    #    print('P is:')
    #    print(P)
    #    print(len(P))
    
    #    print(len(g_diff), 'is the dimension x of g_diff')
    #    print(len(mean_grad), 'is the dimension x of mean_grad')
    
        total_gradient = g_diff.reshape(-1) + P @ mean_grad.reshape(-1)
    
        total_gradient = total_gradient.reshape([-1,3])
        
        x = x.reshape([-1,3])
        y = y.reshape([-1,3])
    
    #    print(total_gradient)
        self.results["energy"] = en2 - en1
        self.results["forces"] = - total_gradient * (ase.units.Hartree / ase.units.Bohr)
    
        with open('x_minus_one.dat', 'w') as x_min_file:
            for atom in x:
                x_min_file.write('%.6f %.6f %.6f\n' % (atom[0], atom[1], atom[2]))
    
        with open('y_minus_one.dat', 'w') as y_min_file:
            for atom in y:
                y_min_file.write('%.6f %.6f %.6f\n' % (atom[0], atom[1], atom[2]))
