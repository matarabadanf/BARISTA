import ase.io
import os
from ase.optimize import FIRE as opti
try:
    from barista.brewer.flavs import CICalculator
    from barista.personnel.javi import Javi
except:
    from flavs import CICalculator
    from javi import Javi

class BrewerJobManager:
    def __init__(self, brewer_templatefile:str='brewer.TEMPLATE'):
        self.brewer_templatefile = brewer_templatefile
    
    def run_job(self):
        self.select_mode()
        self.run_optimization()
        # self.dump_trajectory()
        self.analyze_conical()
        self.cleanup()
        print('Auxiliary files were moved to AUX_FILES/\n\n')
        self.goodbye()

    def read_inputfile(self) -> dict:
        # read the template to determine the geometry and particularities of the calculations such as basis 
        try:
            with open('brewer.TEMPLATE', 'r') as template_file:
                content = template_file.readlines()
            content = [tuple(line.replace('=', ' ').strip().split()) for line in content if len(line.replace('=', ' ').strip().split()) == 2 ]
            input_keyword_dict = dict(content)

            return input_keyword_dict
        except FileNotFoundError:
            print('No brewer.TEMPLATE file was found.')

    def select_mode(self):
        input_keyword_dict = self.read_inputfile()

        print(f'Input keywords were read from {self.brewer_templatefile}\n')

        # set default calculation parameters to the following
        if input_keyword_dict['mode'] == 'ci':
            parameters = dict((
                ('profile', 'ubp'),
                ('n_roots', 10),
                ('iroot', 0),
                ('jroot', 1),
                ('functional', "CAM-B3LYP"),
                ('basis', "cc-pvdz"),
                ('alpha', 0.02),
                ('sigma', 3.5),
                ('calc_nacme', False),
                ('charge', 0),
                ('mult', 1),
                ('program', 'ORCA'),
                ('label', None), 
                )
            )

        # override default parameters with the input file ones 
        for key in input_keyword_dict.keys():
            parameters[key] = input_keyword_dict[key]

        if 'fmax' in  input_keyword_dict.keys():
            self.fmax = input_keyword_dict.keys['fmax']
        else:
            self.fmax = 0.05

        if parameters['label'] is None:
            parameters['label'] = parameters['geom'].replace('.xyz', '') + '_ci'


        # recast to correct type
        try:
            for key in ('n_roots','iroot','jroot', 'charge','mult'):
                parameters[key] = int(parameters[key])
                
            for key in ('alpha', 'sigma'):
                parameters[key] = float(parameters[key])
        except:
            raise TypeError('A parameter float or int could not be casted to the appropriate type')
        
        # recast boolean 
        if str(parameters['calc_nacme']).lower() == 'true':
            parameters['calc_nacme'] = True
        else:
            parameters['calc_nacme'] = False

        #automatically select nroots in case it is not defined in the input 
        if 'n_roots' not in input_keyword_dict.keys():
            parameters['n_roots'] = parameters['jroot'] + 2

        self.parameters = parameters

        # print parameters
        print('Final input parameters: \n')
        for key in self.parameters:
            print(f'    {key} = {self.parameters[key]}')
 
    @property 
    def job_cores(self): 
        try:
            return int(os.popen('echo $SLURM_NTASKS').read())
        except:
            return 1 
    
    def generate_molecule_and_calc(self) -> 'molecule': 
        molecule = ase.io.read(self.parameters['geom'])
        print(f'\nGeometry was read from {self.parameters['geom']}.\n')
        # molecule.calc = CICalculator(atoms=molecule, n_procs=1, **self.parameters)
        molecule.calc = CICalculator(atoms=molecule, n_procs=self.job_cores, **self.parameters)
        print(f'The number of processors used is: {self.job_cores}\n')
        print(f'\nOptimization convergence threshold is: Max effective force < {self.fmax}\n')
        return molecule

    def run_optimization(self):
        molecule = self.generate_molecule_and_calc()
        opt = opti(molecule, trajectory=f"{self.parameters['label']}_ci_search.traj")
        opt.run(fmax=self.fmax)
        # ase.io.write(f"{self.parameters['label']}_ci_search_opt_geom.xyz", molecule)
        print('\n')
        print('**************************')
        print('\n')
        print('OPTIMIZATION HAS CONVERGED')
        print('\n')
        print('**************************\n\n')

    def dump_trajectory(self):
        a = ase.io.trajectory.TrajectoryReader(f"{self.parameters['label']}_ci_search.traj")
        for geom in a:
            ase.io.write(
                filename=f"{self.parameters['label']}_ci_search_opt_geom_traj.xyz", 
                images=geom, format='xyz', 
                append=True
            )
    def analyze_conical(self):
        if self.parameters['calc_nacme']:
            conical = Javi(
                'engrad0.dat',
                'engrad1.dat',
                'nacme.dat'
            )
        else:  
            conical = Javi(
                'engrad0.dat',
                'engrad1.dat',
                'y_minus_one.dat'
            )

        print('\n\nConical intersection characterization after convergence:')
        print(f"P = {conical.p[0]:5.3f} -> {conical.p[1]} conical.\nB = {conical.b[0]:5.3f} -> {conical.b[1]} conical.")

    def cleanup(self):
        os.system('mkdir AUX_FILES')
        os.system('mv *engrad* AUX_FILES')
        os.system('mv *.ener* AUX_FILES')
        os.system('rm -rf __pycache__/')
        os.system('rm *.traj')
        os.system(f'rm {self.parameters["label"]}.xyz')

    def goodbye(self):
        print('*****************************')
        print('          Meow meow          ')
        print('*****************************')
