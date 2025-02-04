#!/usr/bin/env python3
import os

print('Path to a suitable orca launcher:')
orcapath = input()

print(orcapath)

with open('../barista/brewer/flavs.py', 'r') as f:
    flavs_content = f.readlines()

with open('configuration.dat', 'r') as f:
    current_config = f.readlines()

with open('temporary_flavs.py', 'w') as f:
    for line in flavs_content:
        f.write(line.replace(current_config[0].split()[0], orcapath.split()[0]))

with open('configuration.dat', 'w') as f:
    for line in current_config:
        f.write(line.replace(current_config[0].split()[0], orcapath.split()[0]))

os.rename('temporary_flavs.py', '../barista/brewer/flavs.py')

os.system('setup_barista.py')
