#!/usr/bin/env python3

from barista.personnel import rodrigo
from barista.personnel.rodrigo import Rodrigo
import matplotlib.pyplot as plt 
import numpy as np 

## RELOCATE

vert_spec = Rodrigo(f'geom001.in.out')
x_grid,y_grid = vert_spec.export_gaussian_eV()
print(y_grid)
y_grid = np.zeros_like(x_grid)

for i in range(1,200):
    try: 
        print(f'geom{i:03d}.in.out')
        vert_spec = Rodrigo(f'geom{i:03d}.in.out')
        print(vert_spec.export_gaussian_eV(dispersion=0.2)[1])
        y_grid += vert_spec.export_gaussian_eV(dispersion=0.2)[1]
    except:
        print(i)

y_grid /= max(y_grid)
plt.plot(x_grid, y_grid)
plt.xlim([4, 8])
plt.show()
