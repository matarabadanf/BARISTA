import numpy as np
import sys

filename = sys.argv[1]

print(filename)

def gaussian(array, center, height, width=10):
    return array * height * np.exp(-(array-center)**2/(2*width**2))

fi = open(filename, 'r')
cont = fi.readlines()

for line in cont:
    if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
        init_index = cont.index(line)
    elif 'ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS' in line:
        end_index  = cont.index(line)

spectra_peaks = np.array([np.array([float(line.strip().split()[2]), float(line.strip().split()[3])]) for line in cont[init_index +5:end_index-2]])

default_length = 2000

startx = 100

endx = 450

spectra_array = np.zeros(default_length)
energy_array = np.linspace(startx, endx, default_length)

print(spectra_array)
print(energy_array)


for peak in spectra_peaks:
    spectra_array += gaussian(energy_array, peak[0], peak[1])

#spectra_array /= max(spectra_array)
print(spectra_array)

spectra_filename = filename + '.spectra'

spectra_file = open(spectra_filename, 'w')

for i in range(default_length):
    spectra_file.write('%f %f\n' % (energy_array[i], spectra_array[i]))
