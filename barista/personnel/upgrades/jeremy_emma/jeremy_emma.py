from barista.personnel.jeremy import Jeremy
import os 

script_dir = os.path.dirname(os.path.abspath(__file__))

orca_ci = Jeremy(os.path.join(script_dir, 'opt_cont.xyz'))
brew_ci = Jeremy(os.path.join(script_dir, 'beta-carboline_min_1_1_1_ci.xyz'))


orca_ci.compare_internals(brew_ci)

print(orca_ci.atom_labels[1], orca_ci.atom_labels[2], orca_ci.atom_labels[3], orca_ci.atom_labels[4])

print(orca_ci.internal_list[57])
print(orca_ci.internal_values[57])
print(brew_ci.internal_values[57])

print(brew_ci.internal_values[57] - abs(orca_ci.internal_values[57]))
