#!/usr/bin/env python3
import argparse
import sys
import os

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Exctract the final geometry of an optimization run in ORCA.""",
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-f", type=str, required=True, default=42, help="Geometry xyz file to start the optimization"
)

if len(sys.argv) == 1:
   parser.print_help(sys.stderr)
   sys.exit(1)

args = parser.parse_args()


BREWER_PATH = '/gpfs/home/uam/uam121718/bin/BARISTA/barista/brewer/' 

if not BREWER_PATH:
    print('Configure BREWER_PATH with setup_brewer.py')
    exit

if isinstance(args.f, str):
    workdir = args.f.replace('.xyz', '_CI')
    os.system('mkdir %s' % workdir)
    os.system('cp %sbrewer.TEMPLATE %s' % (BREWER_PATH, workdir))
    os.system('mv %s %s' % (args.f, workdir))
    os.system('cp %sopt.py %s' % (BREWER_PATH, workdir))
    os.system('cp %srun_opt.sh %s' % (BREWER_PATH, workdir))
    os.system('cp %sflavs.py %s' % (BREWER_PATH, workdir))
    os.system('cp %s/../personnel/javi.py %s' % (BREWER_PATH, workdir))
    os.system('cp %sBrewerJobManager.py %s' % (BREWER_PATH, workdir))

    with open('%s/brewer.TEMPLATE' % workdir, 'r') as f:
        cont = f.readlines()

    cont = [line.replace('water.xyz', args.f) for line in cont]

     
    with open('%s/brewer.TEMPLATE' % workdir, 'w') as f:
        for line in cont:
            f.write(line)

