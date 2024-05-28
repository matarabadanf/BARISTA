import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]

def gaussian(X, a, b, c):
    return a * np.exp(- (X-b)**2/(2*c**2))

fi = open(filename, 'r')
cont = fi.readlines()

for line in cont:
    if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
        init_index = cont.index(line)
    elif 'ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS' in line:
        end_index  = cont.index(line)

spectra_peaks = np.array([np.array([float(line.strip().split()[2]), float(line.strip().split()[3])]) for line in cont[init_index +5:end_index-2]])
spectra_peaks[:,1] = spectra_peaks[:,1]/max(spectra_peaks[:,1])


X = np.linspace(100, 450, 1000)
final_spectra = np.zeros(1000)


for peak in spectra_peaks:
    print(peak)
    final_spectra += gaussian(X, peak[1], peak[0], 10)
    plt.vlines(x=peak[0], ymin=0, ymax=peak[1], colors='plum')

#plt.scatter(spectra_peaks[:,0], spectra_peaks[:,1]/max(spectra_peaks[:,1]), marker='^', c='violet')


plt.plot(X, final_spectra/max(final_spectra), c='yellowgreen')

plt.title('Absorption Spectra\n$\omega$B97X-D3\nTDA')
plt.xlabel('Wavelength / nm')
plt.ylabel('Oscillator strenght / ')
plt.legend()
plt.savefig(filename[:-7]+'_spectra.jpg', dpi=800)
