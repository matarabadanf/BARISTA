# -*- coding: utf-8 -*-
import numpy as np
import ase.io
from ase.calculators.calculator import FileIOCalculator
import scipy


class CICalculator(FileIOCalculator):
    # Tthe calculator needs a command that is what the program will execute. 
    # Therefore it is needed to prepare the inputs. 
    command = " ~/bin/run_orca.sh orca_engrad.in"
    
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
        geom=None,
        profile:str ='PENALTY',
        n_roots:int =10,
        iroot:int=0,
        jroot:int=1,
        functional:str="CAM-B3LYP",
        basis:str="cc-pvdz",
        alpha:float=0.02,
        sigma:float=3.5,
        n_procs:int=1,
        program: str= 'ORCA',
        calc_nacme: bool = False,  # type: ignore
        charge:int = 0, # type: ignore
        mult:int = 1,
        convergence:str='forced',
        caspt2_params:str='',
        rassccf_params:str='',
        imag:float=0.3,
        shake:float=0.01,
        orbitals:str='',
        rasscf:str = ''
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
        self.n_roots = n_roots
        self.functional = functional
        self.basis = basis
        self.atoms = atoms
        self.alpha = alpha
        self.sigma = sigma
        self.n_procs = n_procs
        self.profile = profile.upper()
        self.charge = charge 
        self.convergence = convergence
        self.mult = mult
        self.imag = imag
        self.orbitals = orbitals
        self.shake = shake 
        self.rasscf = rasscf

        self.program = program.strip().upper()
        # this prepares a current geometry
        with open(self.label + ".xyz", "w") as fd:
            ase.io.write(fd, self.atoms, format="xyz")
        # this saves the original geometry
        with open(self.label + "_original.xyz", "w") as fd:
            ase.io.write(fd, self.atoms, format="xyz")

        # Select the command to run the different programs. 
        if self.program == 'ORCA':
            CICalculator.command = " ~/bin/run_orca.sh orca_engrad.in"
        if self.program == 'MOLCAS':
            CICalculator.command = " ~/bin/run_openmolcas.sh molcas_engrad.in"

        

        self.calc_nacme = calc_nacme
        if self.calc_nacme and [iroot,jroot] == [0,1] and self.program == 'ORCA': 
            print('ORCA supports NACME calculation for states 0 and 1, this will be used')
        elif self.calc_nacme and [iroot,jroot] != [0,1] and self.program == 'ORCA': 
            self.calc_nacme = False 
            print('ORCA only supports NACME calculation for states 0 and 1, nacme will not be calculated')


    def write_input(self, atoms=None, properties=None, system_changes=None):
        FileIOCalculator.write_input(
            self, atoms, properties, system_changes
        )  # because the docs demmand so
        # we will only update the .xyz file as the input for calculations will remain the same throughout the whole calculation

        if self.program == 'ORCA':
            # print('\nOrca input is being generated')
            self.generate_orca_input()
        elif self.program == 'MOLCAS':
            self.generate_molcas_input()

    def generate_orca_input(self):
        '''Generate an orca input. 

        First the geometry of the current step is exported from the ase object. 
        Then the file is read to obtain the xyz coordinates. 
        The input file is written.
        '''
        with open(self.label + ".xyz", "w") as fd:
            ase.io.write(fd, self.atoms, format="xyz")

        with open(self.label + ".xyz", "r") as xyzfile:
            cont = xyzfile.readlines()
            
        coordlist = cont[2:]

        # include nacme in orca calculation. Default is False. 
        nacme_str = ''
        if self.calc_nacme:
            nacme_str = 'NACME True'

        with open('orca_engrad.in', 'w') as inp_file:
            # this prepares the orca inputs

            inp_file.write(
                f"! engrad {self.functional} {self.basis} \n! NoAutostart\n\n\n"
            )

            # Control the TDDFT block
            inp_file.write(
                f"%tddft nroots {self.n_roots} iroot {self.iroot} tda TRUE end\n"
            )

            # Coords, charge and mult of the first calculation
            inp_file.write(f"* xyz {self.charge} {self.mult}\n")
            for atom in coordlist:
                inp_file.write(atom.strip()+'\n')
            inp_file.write("*\n")

            # generate second job
            inp_file.write('$new_job\n')
            inp_file.write(
                f"! engrad {self.functional} {self.basis} \n! NoAutostart\n\n\n"
            )

            # Control the TDDFT block
            inp_file.write(
                f"%tddft nroots {self.n_roots} iroot {self.jroot} tda TRUE {nacme_str} end\n"
            )

            # Coords, charge and mult of the first calculation
            inp_file.write(f"* xyz {self.charge} {self.mult}\n")
            for atom in coordlist:
                inp_file.write(atom.strip()+'\n')
            inp_file.write("*\n\n")

    def parse_orca(self):
        '''Parse the gradient of the lower and upper states and nacme if computed.

        Parses and converts to a 1D arraythe arrays and returns them. 
        Returns an empty array for nacme if the nacme was not requested. 
        '''
        with open('orca_engrad.in.out', 'r') as output_file:
            output_list = output_file.readlines()

        # parse the energies of the two states
        energies = np.array([float(l.strip().split()[-1]) for l in output_list if 'FINAL SINGLE' in l])
        gradients = []

        gradient_indices = [i+3 for i, line in enumerate(output_list) if 'CARTESIAN GRADIENT' in line]
        gradient_ends = [i-1 for i, line in enumerate(output_list) if 'Difference to translation invariance:' in line]

        for gradstart, gradend in zip(gradient_indices, gradient_ends):
            grad_str = [l.strip().split()[3:] for l in output_list[gradstart:gradend]]
            grad = np.array(grad_str, dtype=float)
            gradients.append(grad)

        # Parse nacme if requested 
        if self.calc_nacme:

            nacme_start = [i+4 for i, line in enumerate(output_list) if 'CARTESIAN NON-ADIABATIC COUPLINGS' in line][0]
            nacme_end = [i-1 for i, line in enumerate(output_list) if 'Difference to translation invariance:' in line][-1] 
            nacme_str = [l.strip().split()[3:] for l in output_list[nacme_start:nacme_end]]
            nacme = np.array(nacme_str, dtype=float)
            nacme /= np.linalg.norm(nacme)

        else:
            nacme = np.zeros_like(gradients[0])

        return energies, gradients[0], gradients[1], nacme 


    def generate_molcas_input(self):
        
        with open(self.label + ".xyz", "w") as fd:
            ase.io.write(fd, self.atoms, format="xyz")
    
        # This prepares the Molcas input
        with open('molcas_engrad.in', 'w') as inp_file:
    
            #GATEWAY
            # I have created a shake parameter (the default can be 0.00)
            inp_file.write(
            
            f" &GATEWAY \n Title=molecule \n coord={self.label+'.xyz'} \n basis={self.basis} \
            \n group=Nosym \n ricd \n shake \n {self.shake}" )
    
        
            # SEWARD
            inp_file.write( 
            
            f"&SEWARD \n  ")
            
    
            #  I created another parameter called orbitals in case a RasOrb orbital file exists, otherwise we compute SCF
            if self.orbitals != '':
                inp_file.write(
                f">>COPY $CurrDir/orbitals.RasOrb $WorkDir/$Project.RasOrb \n \
                >>COPY $WorkDir/$Project.RasOrb INPORB \n >>> ENDIF <<<  " )
            else:
                inp_file.write(
                f"&SCF" )
    
    
            #RASSCF specifications for i root
            inp_file.write(    
            f"&RASSCF \n {self.rasscf}")
    
            #CASPT2 specifications (create imaginary shift parameter)
            inp_file.write(f"&CASPT2 \n imag \n {self.imag} \n XMULtistate \n all \n ipea \n 0.0 ")
    
            #ALASKA
            inp_file.write(f"&ALASKA ")

    def parse_molcas(self):
        with open('molcas_engrad.in.out', 'r') as output_file:
            output_list = output_file.readlines()
    
            #Parse the energies of the two states
            energies = []
    
            energies_index = [i for i, line in enumerate(output_list) if 'Point #  42 of   42 done' in line]
            energies_ends = [i+1 for i, line in enumerate(output_list) if 'Total XMS-CASPT2 energies:' in line]
            for energies_ends, energies_index in zip(energies_ends, energies_index):
                en_str = [l.strip().split()[6] for l in output_list[energies_ends:energies_index]]
                en = np.array(en_str, dtype=float)
                energies.append(en)
    
    
    
    
            #Parse the molecular gradients
    
            gradient_indices = [i+4 for i, line in enumerate(output_list) if 'Numerical gradient, root     ' in line]
    
    
            gradient_ends = []
            for indice in gradient_indices:
                gradient_ends.append(next((i for i, line in enumerate(output_list) if '---------------------------------------------' in line and i > indice), None))
    
    
            gradients = []
            for gradstart, gradend in zip(gradient_indices, gradient_ends):
                grad_str = [l.strip().split()[1:] for l in output_list[gradstart:gradend]]
                grad = np.array(grad_str, dtype=float)
                gradients.append(grad)
    
            #no nacme in molcas 
            nacme = np.zeros_like(gradients[0])
    
            return energies, gradients[0], gradients[1], nacme 
    
    def write_results(self):
        '''Writes the results in .dat files, updates energies.dat and trajectory files'''
        if self.program == 'ORCA':
            energies, engrad_0, engrad_1, nacme = self.parse_orca()
        elif self.program == 'MOLCAS':
            energies, engrad_0, engrad_1, nacme = self.parse_molcas()
 
        np.savetxt('ener.dat', energies)
        np.savetxt('engrad0.dat', engrad_0)
        np.savetxt('engrad1.dat', engrad_1)
 
        try: 
            iteration = int(np.loadtxt('iteration.dat'))
        except:
            iteration = 0 
 
        np.savetxt('iteration.dat', np.array([iteration]) + 1)
 
        # calculates nacme if requested 
        if self.calc_nacme:
            np.savetxt('nacme.dat', nacme)
 
        # update logfiles
        with open('energies.dat', 'a') as f:
            f.write(f'{energies[0]:10.6f} {energies[1]:12.8f} {(energies[1] - energies[0])*27.2114:6.4f}\n')
        
        with open(f'{self.label}_traj.xyz', 'a') as f:
            with open(self.label + ".xyz", "r") as xyzfile:
                cont = xyzfile.readlines()
            for line in cont:
                f.write(line)


    def read_results(self):
        '''Writes and reads the .dat files and calculates the effective gradient for the optimizer.'''
        # parse and write the results to files
        self.write_results()
        # print('Calculation results were extracted to .dat files')
        print('\nAt this optimization step, the state is:\n')

        # calculate effective gradient with selected method
        if self.profile == 'PENALTY':
            self.penalty_results()
        elif self.profile == 'GP':
            self.gradient_projection_forces()
        elif self.profile == 'UBP':
            self.ubp_forces()
 
        iteration = float(np.loadtxt('iteration.dat'))
        
        # If the convergence forced, check last 10 energies and if variation is small, stop the calculation.
        if iteration > 11 and self.convergence == 'forced':
            ener = np.loadtxt('energies.dat')[-10:]
            print(f'\nener0 deviation in the last 10 cycles =  {np.std(ener[:,0]):10.8f} hartree')
            print(f'ener1 deviation in the last 10 cycles =  {np.std(ener[:,1]):10.8f} hartree')
            print(f'D_E deviation in the last 10 cycles =  {np.std(ener[:,2]):10.8f} eV\n')
            if np.std(ener[:,0]) < 0.00001 and np.std(ener[:,1]) < 0.00001 and np.std(ener[:,1]) < 0.01:
                # if np.std(ener[:,0]) < 0.1 and np.std(ener[:,1]) < 0.1 and np.std(ener[:,1]) < 0.1:
                print('Optimization converged due to lack of change in the last 10 iterations\n')
                self.results["forces"] = np.zeros_like(np.copy(self.results['forces']))


    def penalty_results(self):
        '''Calculates effective force based in the penalty algorithm'''
        # management of the energies
        en1, en2 = np.loadtxt("ener.dat")

        en = (en2 - en1) ** 2 / (en2 - en1 + self.alpha)
        en = self.sigma * en + (en2 + en1) / 2.0

        # management of the gradients
        grad1 = np.loadtxt("engrad0.dat")
        grad2 = np.loadtxt("engrad1.dat")

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
        '''Calculates effective force based in the gradient projection algorithm'''
        # management of the energies
        en1, en2 = np.loadtxt("ener.dat")

        # management of the gradients
        grad1 = np.loadtxt("engrad0.dat")
        grad2 = np.loadtxt("engrad1.dat")
        
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
        '''Calculates effective force based in the updated branching plane algorithm with or without nacme'''
        # management of the energies
        en1, en2 = np.loadtxt("ener.dat")
    
        # management of the gradients
        grad1 = np.loadtxt("engrad0.dat").reshape(-1)
        grad2 = np.loadtxt("engrad1.dat").reshape(-1)
    
        # check if previous iteration values are in memory
        try:
            x_minus_one = np.loadtxt('x_minus_one.dat')
            y_minus_one = np.loadtxt('y_minus_one.dat')
            first_iter = False
        except:
            first_iter = True

        # x is the unit vector parallel to the difference gradient vector
        x = (grad1-grad2).reshape(-1)
        x /= np.linalg.norm(x)

        # the definition of y can be: first step, updated, exact
        # in the first iteration, the (normalized) mean energy gradient vector is chosen as y_0    
        # otherwise the x and y previous steps are loaded
        if first_iter:
            y_space = scypi.linalg.null_space(x)
            y =  y_space[:, 0]

        elif not self.calc_nacme:
            x_minus_one = x_minus_one.reshape(-1)
            y_minus_one = y_minus_one.reshape(-1)
            
            # the y vector is calculated from the previous step
            y = (np.dot(y_minus_one, x) * x_minus_one - np.dot(x_minus_one, x) * y_minus_one) / (np.dot(y_minus_one, x)**2 + np.dot(x_minus_one, x)**2)**0.5
            y /= np.linalg.norm(y)

        else:
            y = np.loadtxt('nacme.dat').reshape(-1)
            y /= np.linalg.norm(y)

        # ensure correct dimensions
        x = x.reshape(-1)
        y = y.reshape(-1)

        # All vectors x, y, x_minus and y_minus should have dimension 1 
    
        # generate the effective gradients
        g_diff = 2 * (en1 - en2) * x
        mean_grad = (grad1 + grad2)/2
    
        # Generate the projector here the 0.3 weight was chosen semi-arbitrarily. 
        P = np.identity(len(grad1.reshape(-1))) - np.outer(x,x) - 0.3 * np.outer(y,y)
    
        # calculate the total effective gradient with the projector
        total_gradient = g_diff.reshape(-1) + P @ mean_grad.reshape(-1)
        total_gradient = total_gradient.reshape([-1,3])
        
        # reshape to save the data
        x = x.reshape([-1,3])
        y = y.reshape([-1,3])
    
        # feed the optimizer with energies and gradientsprint(total_gradient)
        self.results["energy"] = en2
        self.results["forces"] = - total_gradient * (ase.units.Hartree / ase.units.Bohr)
    
        # save results for next iteration
        if not self.calc_nacme:
            np.savetxt('x_minus_one.dat', x)
            np.savetxt('y_minus_one.dat', y)
