from barista.personnel.javi import Javi
import numpy as np 
import os 

# I have to use CI derivative coupling vector from molcas. The rest follows:

script_dir = os.path.dirname(os.path.abspath(__file__))

ethylene_path = os.path.join(script_dir, 'ethylene_b')
grad0_path = os.path.join(ethylene_path, 'engrad0.dat')
grad1_path = os.path.join(ethylene_path, 'engrad1.dat')
nacme_path = os.path.join(ethylene_path, 'nacme2.dat')

j =  Javi(
    grad0_path,
    grad1_path,
    nacme_path
)

ga = np.loadtxt(grad0_path).flatten()
gb = np.loadtxt(grad1_path).flatten()
hab = np.loadtxt(nacme_path).flatten()
g_ab = ( gb - ga )/2
print(f'g_ab and j.g_ab are equal: {np.all(np.isclose(j.g_ab, g_ab))}')

tan_2_beta = 2 * np.dot(g_ab, hab) / (np.dot(g_ab, g_ab) - np.dot(hab, hab))
beta_2 = np.atan(tan_2_beta) 
beta = beta_2 / 2 
print(f'beta and j.beta are equal: {np.all(np.isclose(j.beta, beta))}')
print((beta + np.pi/2 % np.pi) /np.pi*360, (j.beta % np.pi) /np.pi*360)

ref_pitch = 0.1026
ref_asymm = 0.4886
ref_sigma = 0.5668
ref_theta = 0 

print(f'Reference Pitch and j.pitch differ in          : {abs((j.pitch - ref_pitch)):9.5f}')
print(f'Reference asymmetry and j.asymmetry differ in  : {abs((j.asymmetry - ref_asymm)):9.5f}')
print(f'Reference sigma and j.sigma differ in          : {abs((j.sigma - ref_sigma)):9.5f}')
print(f'Reference theta_s and j.theta_s differ in      : {abs((j.theta_s % np.pi - ref_theta % np.pi)):9.5f}')
print('Conical reference characterization')
print('(a) 0.25, 1.02 peaked single-path ')
print('Conical calculated characterization')
print(f'(a) {j.p[0]:3.2}, {j.b[0]:5.3}, {j.p[1]} {j.b[1]}')
