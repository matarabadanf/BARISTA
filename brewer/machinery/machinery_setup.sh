#!/bin/bash

flavour=$1
name=$2 
MACHINERY_PATH=`dirname "$0"`

echo import ase.io > opt.py
echo from ase.optimize import FIRE as opt >> opt.py
echo from flavours import $flavour as calc >> opt.py
echo  >> opt.py
echo 'geom=ase.io.read("geom.xyz")'  >> opt.py
echo 'cal=calc(label="'${name}'")'  >> opt.py
echo geom.calc = cal  >> opt.py
echo   >> opt.py
echo 'geom.get_forces()'  >> opt.py
echo   >> opt.py
echo 'opt=opti(geom, trajectory="'$name'.traj")'  >> opt.py
echo 'opt.run()'  >> opt.py
echo 'ase.io.write("'$name'_opt_geom.xyz",geom)'  >> opt.py

cp $MACHINERY_PATH/flavours.py .

cat $MACHINERY_PATH/run_opt.sh > run_opt.sh

