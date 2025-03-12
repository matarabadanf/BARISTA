from barista.personnel.javi import Javi 
import numpy as np 


j = Javi('engrad0.dat', 'engrad1.dat', 'nac.dat')


x = np.loadtxt('x_vec.dat').reshape(-1)
y = np.loadtxt('y_vec.dat').reshape(-1)

print(f'Angle is {x @ y}')
beta = j.beta

print(f'Beta value is {beta}')


for del_b in range(0,4):
    b_prime = (beta+ np.pi/4*del_b) % np.pi
    print(f'Beta value is {b_prime:8.5f} (originally {beta:8.5f})')
    j.beta = b_prime
    print(j.x)
    print(f'Angle is {j.x @ j.y}')
    print(j.pitch)
    print(j.asymmetry)
# 
# print(j.p)
# print(j.b)
