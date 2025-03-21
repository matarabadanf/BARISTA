#!/usr/bin/bash


chmod 700 setup_launchers.py

chmod 700 ../barista/personnel/*.py 
chmod 700 ../barista/brewer/*.py
brewer_path="$(realpath "$(pwd)/../barista/brewer/")"
sed -i "s|BREWER_PATH = '.*'|BREWER_PATH = '${brewer_path}/'|" ${brewer_path}/brewer_setup.py
