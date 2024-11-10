import ase.io
import os
from contextlib import redirect_stdout
from ase.optimize import FIRE as opti
from flavs import Penalty 

# determine the number of cores for the task  
n_procs = int(os.popen('echo $SLURM_NTASKS').read())

# read the template to determine the geometry and particularities of the calculations such as basis 
try:
    with open('brewer.TEMPLATE', 'r') as template_file:
        content = template_file.readlines()

    content = [tuple(line.replace('=', ' ').strip().split()) for line in content if len(line.replace('=', ' ').strip().split()) == 2 ]
    cont = dict(content)

except:
    print('There is no brewer.TEMPLATE file in this directory. "geom.xyz" is the default geometry filename to run the calculations (this will likely fail).')

# setting up the parameters
if 'geom' not in cont.keys():
    cont['geom'] = 'geom.xyz'

if 'label' not in cont.keys():
    cont['label'] = cont['geom'].replace('.xyz', '')

print(cont)

# generate atoms object, associate calculator and get forces 
molecule = ase.io.read(cont['geom'])
molecule.calc = Penalty(label=cont['label'], atoms=molecule, n_procs=n_procs)
molecule.get_forces()


with open('optimization.log', "w") as f:
    with redirect_stdout(f):
        opt = opti(molecule, trajectory="%s_ci_search.traj" % cont['label'])
        opt.run()
        ase.io.write("%s_ci_search_opt_geom.xyz" % cont['label'], molecule)

a = ase.io.trajectory.TrajectoryReader("%s_ci_search.traj" % cont['label'])

for geom in a:
    ase.io.write(filename="%s_ci_search_opt_geom_traj.xyz" % cont['label'], images=geom, format='xyz', append=True)

os.system('mkdir AUX_FILES')
os.system('mv *engrad* AUX_FILES')
os.system('mv *.energy* AUX_FILES')
os.system('rm -rf __pycache__/')
