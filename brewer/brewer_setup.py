import subprocess 
import os 


BREWER_PATH = os.path.dirname(os.path.abspath(__file__))


print(BREWER_PATH)

BREWER_PATH = os.path.join(BREWER_PATH)

print(BREWER_PATH)

basis, functional, geometry, name, tp, state0, state1, nroots, flavour= [None for i in range(9)]

try:
    with open('brewer_TEMPLATE', 'r') as inp:
        cont = inp.readlines()
except:
    os.system('cp ' + str(BREWER_PATH+'/brewer_TEMPLATE_empty') + ' brewer_TEMPLATE')
for line in cont:
    print(line.strip().upper())
    if 'BASIS' in line.upper():
        basis = line.strip().upper().replace('BASIS', '').replace('=', '')
    elif 'FUNCTIONAL' in line.upper():
        functional = line.strip().upper().replace('FUNCTIONAL', '').replace('=', '')
    elif 'GEOMETRY' in line.upper():
        geometry = line.strip()[9:]
    elif 'NAME' in line.upper():
        name = line[5:].strip()
    elif 'STATE0' in line.upper():
        state0 = line.strip().upper().replace('STATE0', '').replace('=', '')
    elif 'STATE1' in line.upper():
        state1 = line.strip().upper().replace('STATE1', '').replace('=', '')
    elif 'TYPE' in line.upper():
        tp = line.strip().upper().replace('TYPE', '').replace('=', '').split()[0]
        flavour = line.strip().upper().replace('TYPE', '').replace('=', '').split()[1].lower()

print(basis, functional, geometry, name, tp, state0, state1, flavour, nroots)

os.system('cp ' + geometry + ' geom.xyz')
os.system('cp ' + str(BREWER_PATH+'/machinery/run_orca5.0.3') + ' .')

with open(str(BREWER_PATH+'/machinery/opt_template.py'), 'r') as opttemplate:
    cont = opttemplate.readlines()

with open('opt.py', 'w') as optfile:
    for line in cont:
        optfile.write(line.replace('$name', name).replace('$flavour', flavour))


if tp == 'CI':
    if not nroots:
        nroots = 10

    os.system(str(BREWER_PATH+'/flavours/CI/engrad_mixer.sh') + ' '  + basis + ' ' + functional + ' ' + str(state0) + ' ' + str(1))
    os.system(str(BREWER_PATH+'/flavours/CI/engrad_mixer.sh') + ' '  + basis + ' ' + functional + ' ' + str(state1) + ' ' + str(2))

#    os.system(str(BREWER_PATH+'/machinery/machinery_setup.sh') + ' '  + flavour + ' ' + name)
    os.system('cp ' + str(BREWER_PATH+'/flavours/CI/run_template.sh') + ' .')
    os.system('cp ' + str(BREWER_PATH+'/machinery/run_opt.sh') + ' .')
    os.system('cp ' + str(BREWER_PATH+'/machinery/flavours.py') + ' .')
