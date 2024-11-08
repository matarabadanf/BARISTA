import ase.io
import os
from contextlib import redirect_stdout
from ase.optimize import FIRE as opti
from flavs import Penalty 

molecule = ase.io.read("water.xyz")
molecule.calc = Penalty(label='water', atoms=molecule)

molecule.get_forces()


with open('optimization.log', "w") as f:
    with redirect_stdout(f):
        opt = opti(molecule, trajectory="water_ci_search.traj")
        opt.run()
        ase.io.write("water_ci_search_opt_geom.xyz", molecule)

a = ase.io.trajectory.TrajectoryReader("water_ci_search.traj")

for geom in a:
    ase.io.write(filename="water_ci_search_opt_geom_traj.xyz", images=geom, format='xyz', append=True)

os.system('mkdir AUX_FILES')

os.system('mv *engrad* AUX_FILES')
os.system('mv *.energy* AUX_FILES')
os.system('rm -rf __pycache__/')
