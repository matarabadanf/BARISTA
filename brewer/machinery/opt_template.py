import ase.io
from ase.optimize import FIRE as opti
from flavours import $flavour as calc

geom=ase.io.read("geom.xyz")
cal=calc(label="$name")
geom.calc = cal

geom.get_forces()

opt=opti(geom, trajectory="$name.traj")
opt.run()
ase.io.write("$name_opt_geom.xyz", geom)

a = ase.io.trajectory.TrajectoryReader("$name.traj")

for geom in a:
    ase.io.write(filename="$name_opt_geom_traj.xyz", images=geom, format='xyz', append=True)
