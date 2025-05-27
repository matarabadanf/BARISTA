#!/usr/bin/env python3
from functools import cached_property
import argparse
import os
from typing import Tuple
import numpy as np
import plotly.graph_objects as go
import sys
import matplotlib.pyplot as plt 

# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Takes the lower and upper state gradients in a CI with the NAC vector and characterizes the type oc CI.\
        \nResults can be plotted interactively using plotly.""",
    epilog='It requires g0, g1 and cdv OR molcas output file',
)


core_args = parser.add_argument_group('\nCore data')

core_args.add_argument(
    "-molcas_file",
    type=str,
    default='0.0',
    help="Molcas .out file with gradients and cdv.",
)

core_args.add_argument(
    "-g0",
    type=str,
    help="Lower state gradient file.",
)

core_args.add_argument(
    "-g1",
    type=str,
    help="Higher state gradient file.",
)

core_args.add_argument(
    "-cdv",
    type=str,
    help="cdv vector file.",
)

#############################################
#          Geometries and forces            #
#############################################

forces_args = parser.add_argument_group('\nVector visualization and geometry displacement')

forces_args.add_argument(
    "-ref_xyz",
    default='None',
    type=str,
    help="Reference xyz file for forces and displacement.",
)

forces_args.add_argument(
    "--generate_forces",
    default='None',
    type=str,
    help="Generate xyz file with  x, y, and -s_ab vectors for visualization."
)

forces_args.add_argument(
    "--generate_displacements",
    action='store_true',
    required=False,
    help="Generate displaced geometries in x, -x, y, -y -s_ab and mirrored s_a. A folder will be generated with these geometries.",
)

forces_args.add_argument(
    "-disp_mod",
    default=0.1,
    type=float,
    help="Displacement module of the vectors to generate displaced geometries. Default is 0.1.",
)

#############################################
#                Visualization              #
#############################################

visualization_args = parser.add_argument_group('\nPES visualization')

visualization_args.add_argument(
    "--plot_lower",
    action='store_true',
    required=False,
    help="Plot lower state pes in 2D.",
)

visualization_args.add_argument(
    "--plot_upper",
    action='store_true',
    required=False,
    help="Plot upper state pes in 2D.",
)

visualization_args.add_argument(
    "--interactive",
    action='store_true',
    required=False,
    help="Open interactive plot.",
)

class Javi:

    def __init__(
        self,
        ga_filename: str,
        gb_filename: str,
        hab_filename: str,
        ci_energy:float = 0.0
    ) -> None:
        self._ga_filename = ga_filename
        self._gb_filename = gb_filename
        self._hab_filename = hab_filename
        self._ci_energy = ci_energy
        self.rotation_counter = 0

        self._init()

    def _init(self) -> None:
        self._load_vectors()

        for i in range(0,4):
            if self.is_rotation_needed(self.asymmetry, self.theta_s):
                self._rotate_for_beta()
            else:
                break
        if self.rotation_counter == 4:
            print('#############################################')
            print('#   WARNING: THERE IS NO ROTATION OF Beta   #')
            print('#     SUCH AS PROPERTIES ARE FULFILLED      #')
            print('#            HANDLE WITH CARE               #')
            print('#                                           #')
            print('#     THIS CONICAL MIGHT NOT BE A MECI      #')
            print('#############################################')

        print(f'\nBeta has been rotated {self.rotation_counter % 4} * pi/2. Beta = {self.beta:8.4f} ({self.beta/2/np.pi*360:8.4f} deg)')
        # print(f'beta      = {self.beta%np.pi:8.4f}')
        print(f'\npitch     = {self.pitch:8.4f}')
        print(f'asymmetry = {self.asymmetry:8.4f}')
        print(f'tilt      = {self.sigma:8.4f}')
        print(f'theta     = {self.theta_s:8.4f} ({self.theta_s/2/np.pi*360:8.4f} deg)')

        print(f'\nWith the rotation, the conditions 0 < theta < pi/2 and 0 < asymmetry are fulfiled: {not self.is_rotation_needed(self.asymmetry, self.theta_s)}')

        print('\nP and B values are:\n')
        print(f'{self.p[0]:8.4f}, {self.p[1]}')
        print(f'{self.b[0]:8.4f}, {self.b[1]}\n')

        print('\nX and Y Vectors are:')
        print('\nx vector:')
        for row in self.x.reshape(-1, 3):
            print(" ".join(f"{val:12.8f}" for val in row))

        print('\ny vector:')
        for row in self.y.reshape(-1, 3):
            print(" ".join(f"{val:12.8f}" for val in row))
   
    @classmethod
    def from_molcas(cls, output_filename: str):
            with open(f'{output_filename}', 'r') as output_file:
                output_list = output_file.readlines()
    
            #Parse the molecular gradients
    
            output_list =[line.strip() for line in output_list]
            
            gradient_pre_indices = [i for i, line in enumerate(output_list) if 'Molecular gradients' in line]
            pt_gradient_pre_indices = [i for i, line in enumerate(output_list) if 'Numerical gradient     ' in line]
    
#            print(f'Molecular gradient indices: {gradient_pre_indices}')
            if len(gradient_pre_indices) > 2:
                print('More than 2 gradients were calculated, assuming the last two are the highest theory level ones.')
#            print(f'Perturbative gradient indices: {pt_gradient_pre_indices}')

            gradient_indices = []
            gradient_ends = []

            for index in gradient_pre_indices[-2:]:
                # print(f'Molecular gradient found in line {index}\n')
                lines = []
                for linenumber, line in enumerate(output_list[index:index +20]):
                    # print(line)
                    if '-----' in line:
                        # print(f'Line was detected in line {linenumber + index}')
                        lines.append(int(linenumber + index))
                        # print(f'lines was appended {int(linenumber + index)}')
                        # print(f'\nlines = {lines}')
                # print(f'Lines are {lines}')
                gradient_indices.append(lines[-2]+1)
                gradient_ends.append(lines[-1])

            cdv_pre_index = [i for i, line in enumerate(output_list) if 'CI derivative coupling' in line][-1]
            lines = []
            for linenumber, line in enumerate(output_list[cdv_pre_index:cdv_pre_index + 20]):
                # print(line)
                if '-----' in line:
                    # print(f'Line was detected in line {linenumber + cdv_pre_index}')
                    lines.append(int(linenumber + cdv_pre_index))
                    # print(f'lines was appended {int(linenumber + cdv_pre_index)}')
                    # print(f'\nlines = {lines}')
            # print(f'Lines are {lines}')
            cdv_start, cdv_end = lines[-2]+1 , lines[-1]

            # print(gradient_indices) 
            # print(gradient_ends) 


            # print('Here stuff must be ok')

            # for i in range(0,2):
            #     print(f'start, end = {gradient_indices[i]}, {gradient_ends[i]}')
            #     print(output_list[gradient_indices[i]:gradient_ends[i]])

            # print(output_list[cdv_start:cdv_end])

            

            grad_0 = np.array([
                line.strip().split()[1:] for line in output_list[gradient_indices[0]:gradient_ends[0]]
                ], dtype=float)

            grad_1 = np.array([
                line.strip().split()[1:] for line in output_list[gradient_indices[1]:gradient_ends[1]]
                ], dtype=float)

            cdv = np.array([
                line.strip().split()[1:] for line in output_list[cdv_start:cdv_end]
                ], dtype=float)


            # gradients = []
            # for gradstart, gradend in zip(gradient_indices, gradient_ends):
            #     grad_str = [l.strip().split()[1:] for l in output_list[gradstart:gradend]]
            #     grad = np.array(grad_str, dtype=float)
            #     gradients.append(grad)
    
            # #no nacme in molcas 
            # nacme = np.zeros_like(gradients[0])

            # print(gradients)

            np.savetxt('engrad_0.dat', grad_0)
            np.savetxt('engrad_1.dat', grad_1)
            np.savetxt('cdv.dat', cdv)
    
            return cls('engrad_0.dat', 'engrad_1.dat', 'cdv.dat') 
 


    @property
    def ci_energy(self) -> float:
        return float(self._ci_energy)

    def _load_vectors(self) -> None:
        self._ga = np.loadtxt(self._ga_filename).reshape(-1)
        # self._ga /= np.linalg.norm(self._ga)

        self._gb = np.loadtxt(self._gb_filename).reshape(-1)
        # self._gb /= np.linalg.norm(self._gb)

        self._h_ab = np.loadtxt(self._hab_filename).reshape(-1)
        # self._h_ab /= np.linalg.norm(self._h_ab)

    @cached_property
    def g_ab(self) -> np.ndarray:
        return 0.5 * np.copy(self._gb - self._ga)


    @cached_property
    def s_ab(self) -> np.ndarray:
        return 0.5 * np.copy(self._gb + self._ga)
    
    @cached_property
    def h_ab(self) -> np.ndarray:
        return np.copy(self._h_ab)

    @cached_property
    def _pre_beta(self) -> np.array:
        # tan_2_beta = 2 * (np.dot(self.g_ab, self.h_ab)) / (np.linalg.norm(self.g_ab) - np.linalg.norm(self.h_ab))
        # return np.arctan(tan_2_beta)
        tan_2_beta = 2 * (np.dot(self.g_ab, self.h_ab)) / (np.dot(self.g_ab, self.g_ab) - np.dot(self.h_ab, self.h_ab))
        beta_2 = np.arctan(tan_2_beta)
        beta = beta_2 / 2
        # return beta  
        beta = 0.5 * np.arctan2(2 * np.dot(self.g_ab, self.h_ab), (np.dot(self.g_ab, self.g_ab) - np.dot(self.h_ab, self.h_ab)))

        # if beta < 0: 
        #     beta += np.pi/4


        # print(f'Beta, pi/4 = {beta} {np.pi/4}')

        # print(f'Beta mod pi/4 = {beta % np.pi/4}')
        # print(f'Beta before adjustment = {beta}')
        # while beta > np.pi/4:
        #     beta -= np.pi/2
        # beta =  beta + np.pi/4
        # print(f'Beta after adjustment = {beta}')


        beta_2 = np.arctan2(
            2 * np.dot(self.g_ab, self.h_ab),
            np.dot(self.g_ab, self.g_ab) - np.dot(self.h_ab, self.h_ab)
        )    # returns 

        beta = beta_2/2

        return [beta + i*np.pi/2 for i in range(0, 4)]
        # return 0.5 * np.arctan2(2 * np.dot(self.g_ab, self.h_ab), (np.dot(self.g_ab, self.g_ab) - np.dot(self.h_ab, self.h_ab)))

    def _rotate_for_beta(self):
        # print(f'self.rotation_counter is {self.rotation_counter}')
        self.rotation_counter += 1 

    @property
    def beta(self) -> float:
        # print(f'Rotation of {self.rotation_counter} pi halves')
        return self._pre_beta[self.rotation_counter % 4] # 2 works for x and y unit vector sign in ethylene a and b  

    def is_rotation_needed(self, asymmetry:float, theta_s:float) -> bool:
        zero_condition = 0 < self.theta_s and abs(0 - self.theta_s) < 10**-12
        
        # print(f'self.asymmetry > 0     = {self.asymmetry > 0}, {self.asymmetry}')
        # print(f'0 < self.theta_s       = {self.theta_s > -10**-12}, {self.theta_s}')
        # print(f'self.theta_s < np.pi/2 = {self.theta_s < np.pi/2}')

        return not (asymmetry > 0 and theta_s > -10**-12  and theta_s < np.pi/2)

        return not (self.asymmetry > 0 and self.theta_s > -10**-12  and self.theta_s < np.pi/2)
        

    @property
    def _g_tilde(self) -> np.ndarray:
        return self.g_ab * np.cos(self.beta) + self.h_ab * np.sin(self.beta)

    @property
    def _h_tilde(self) -> np.ndarray:
        return self.h_ab * np.cos(self.beta) - self.g_ab * np.sin(self.beta) # original

    @property
    def x(self) -> np.ndarray:
        # return np.copy(self._g_tilde / np.dot(self._g_tilde, self._g_tilde)**0.5)
        return np.copy(self._g_tilde / np.linalg.norm(self._g_tilde))

    @property
    def y(self) -> np.ndarray:
        # return np.copy(self._h_tilde / np.dot(self._h_tilde, self._h_tilde)**0.5)
        return np.copy(self._h_tilde / np.linalg.norm(self._h_tilde))

    @property
    def pitch(self) -> float:
        ''' Pitch \\delta_gh. '''

        return np.sqrt( 1/2 * (np.dot(self._g_tilde, self._g_tilde) + np.dot(self._h_tilde, self._h_tilde)))

    @property
    def asymmetry(self) -> float:
        ''' Asymmetry \\Delta_gh. '''

        asym = (np.dot(self._g_tilde, self._g_tilde) - np.dot(self._h_tilde, self._h_tilde)) / (np.dot(self._g_tilde, self._g_tilde) + np.dot(self._h_tilde, self._h_tilde))

        # return abs(float(asym))
        return float(asym) # original

    def energy_difference(self, x:float, y:float) -> float:
        return 2 * self.pitch * np.sqrt((x**2 + y**2) + self.asymmetry * (x**2 - y**2))
 
    def average_energy(self, x:float, y:float) -> float:
        return self.ci_energy + x * np.dot(self.s_ab, self.x) + y * np.dot(self.s_ab, self.y)

    def E_A(self, x:float, y:float) -> float:
        s_x = np.dot(self.s_ab, self.x) / self.pitch 
        s_y = np.dot(self.s_ab, self.y) / self.pitch
        # return x + y  
        return self.ci_energy + self.pitch * (x * s_x + y * s_y - ((x**2 + y **2) + self.asymmetry * (x**2 -y**2))**0.5)


        diff = self.energy_difference(x, y)
        aver = self.average_energy(x,y)
        return aver - diff

    def E_B(self, x:float, y:float) -> float:
        diff = self.energy_difference(x, y)
        aver = self.average_energy(x,y)
        return aver + diff

    @property
    def theta_s(self, n_points:int = 1800, radius:float=1):

        angles = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
        x_points = radius * np.cos(angles)
        y_points = radius * np.sin(angles)

        pairs = [[x_points[i], y_points[i]] for i in range(n_points)]

        samples = []

        for pair in pairs:
            x, y = pair
            sample = [x, y, self.average_energy(x, y)]
            # print(f'\n\n\n The xyz values of tilt is: {x:5.2f}, {y:5.2f}, {self.average_energy(x,y):5.2f}')
            samples.append(sample)

        x,y,z = max(samples, key=lambda item: item[2])
        # print(f'\n\n\n The xyz values of the max tilt is: {x:5.2f}, {y:5.2f}, {z:5.2f}')

        vec = np.array([x,y])
     
        theta = np.arctan2(y,x)

        return theta

        print(vec)
        cos_theta = np.dot(vec, np.array([1,0])) / np.linalg.norm(vec)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)

        theta = np.arccos(cos_theta)
 
        # theta = -theta if vec[0] < 0 else theta

        # print(f'The angle theta is {theta} or {(theta+np.pi/2)%np.pi}')
        

        # print(f'The angle of the maximum tilt is {theta:5.3} Radians or {theta*360/np.pi:3.5} degrees')
        
        return theta


    @property
    def sigma(self):
        s_x = np.dot(self.s_ab, self.x) / self.pitch
        s_y = np.dot(self.s_ab, self.y) / self.pitch
        
        return np.sqrt(s_x**2 + s_y**2)

    @property
    def p(self) -> Tuple[float, str]:

        p = self.sigma**2 / (1 - self.asymmetry**2) * (1 - self.asymmetry * np.cos(2*self.theta_s))

        p_type = 'Peaked' if p < 1 else 'Sloped'

        return (p, p_type)

    @property
    def b(self) -> Tuple[float, str]:

        b = (self.sigma**2/(4*self.asymmetry**2)) **(1/3) * (((1+self.asymmetry)*np.cos(self.theta_s)**2)**(1/3) + ((1-self.asymmetry)*np.sin(self.theta_s)**2)**(1/3))

        b_type = 'Bifurcating' if b < 1 else 'Single Path'

        return (b, b_type)

    def plot_CI(self, max_grid:float=1, filename:str = 'NONE.html'):
    
        x = np.linspace(-max_grid, max_grid, 501)
        y = np.linspace(-max_grid, max_grid, 501)
        
        X, Y = np.meshgrid(x, y)
        e_a = self.E_A(X, Y)
        e_b = self.E_B(X, Y)
        b_p = self.ci_energy + np.zeros_like(X)
        e_mean = self.average_energy(X, Y)
    
    
        threshold = 0.00001 * max_grid
        intersection_mask = np.abs(b_p - e_mean) < threshold
        
        # Extract the x, y coordinates where the intersection occurs
        x_intersect = X[intersection_mask]
        y_intersect = Y[intersection_mask]
        z_intersect = b_p[intersection_mask] 
    
        # Floor for contour 
        z_min = min(e_a.flatten()) - (max(e_a.flatten()) - min(e_a.flatten())) * 0.05  # Slightly below the minimum

        # Create a constant z-plane at the bottom to display the contour
        z_floor = np.ones(e_a.shape) * z_min
        
        max_x = max_grid * np.cos(self.theta_s)
        max_y = max_grid * np.sin(self.theta_s)
        max_value = self.average_energy(max_x,max_y)

        tilt_direction = np.array([-max_x, -max_y, -max_value]) * self.sigma

        max_x, max_y, max_value = tilt_direction
    
        # rgba(100, 120, 140, 0.7)
        # Create a 3D surface plot
        fig = go.Figure(data=[
            go.Surface(z=e_a, x=X, y=Y, showscale=False, name='Surface 1', opacity=0.8, colorscale=[[0, 'rgba(100, 120, 140, 0.7)'], [1, 'rgba(100, 120, 140, 0.7)']],
                 contours = {
                     "z": {
                     "show": True,
                     "start": np.min(e_a),
                     "end": np.max(e_a),
                     "size": (np.max(e_a) - np.min(e_a)) / 10,
                     "color":"black",
                     "width": 2,
                     "usecolormap": False,
                     }
                 }),
            go.Surface(z=e_b, x=X, y=Y, showscale=False, name='Surface 2', opacity=0.8, colorscale=[[0, 'rgba(173, 216, 230, 0.7)'], [1, 'rgba(173, 216, 230, 0.7)']],
                 contours = {
                     "z": {
                     "show": True,
                     "start": np.min(e_b),
                     "end": np.max(e_b),
                     "size": (np.max(e_b) - np.min(e_b)) / 10,
                     "color":"black",
                     "width": 2,
                     "usecolormap": False,
                     }
                 }),
            #go.Surface(z=e_mean, x=X, y=Y, showscale=False, name='Branching plane', opacity=0.5, colorscale=[[0, 'rgba(255, 0, 0, 0.7)'],[1, 'rgba(255, 0, 0, 0.7)']]),
            #go.Surface(z=b_p,    x=X, y=Y, showscale=False, name='Mean energy plane', opacity=0.5, colorscale=[[0, 'rgba(0, 255, 0, 0.7)'],[1, 'rgba(0, 255, 0, 0.7)']]),

            # Vectors
            # go.Scatter3d(x=x_intersect, y=y_intersect, z=z_intersect, mode='lines+text', line=dict(color='rgba(0, 0, 0, 0.5)', width=4), name='BP and ME intersection'),
            #go.Scatter3d(x=[0, max_x], y=[0, max_y], z=[0, max_value], mode='lines+text',line=dict(color='black', width=4), name = r'Max tilt direction', text=['',  'Max tilt direction'], textfont=dict(color='black')),
            #go.Scatter3d(x=[0, 1*max_grid], y=[0, 0], z=[0, 0], mode='lines+text',line=dict(color='black', width=4), name = r'$\hat{\mathbf{x}}$', text=['', 'x'], textfont=dict(color='black')),
            #go.Scatter3d(x=[0, 0], y=[0, 1*max_grid], z=[0, 0], mode='lines+text',line=dict(color='black', width=4), name = r'$\hat{\mathbf{y}}$', text=['', 'y'], textfont=dict(color='black')),
    
            # Dummy scatter traces for the legend to reflect the surface color
            # go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(100, 120, 140, 0.7)', size=10), name='Surface 1'),
            # go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(173, 216, 230, 0.7)', size=10), name='Surface 2'),
            # go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(255, 0, 0, 0.7)', size=10), name='Branching plane'),
            # go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(0, 255, 0, 0.7)', size=10), name='Mean energy plane'),
    
    
        ])
    
        fig.update_layout(
            # title='Conical intersection representation',
            # title_font=dict(size=28),
            width=1200,  # Make figure larger
            height=900,  # Make figure larger
            scene=dict(
            xaxis_title=None,
            yaxis_title=None,
            zaxis_title=None,
            aspectmode='cube',
            xaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                zeroline=False,
                title=None
            ),
            yaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                zeroline=False,
                title=None
            ),
            zaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                zeroline=False,
                title=None
            ),
            aspectratio=dict(x=1, y=1, z=1)
            ),
            showlegend=True,
            legend=dict(
            font=dict(size=22),
            traceorder='normal',
            )
        )
        
        fig.write_html(filename)

        fig.show()    

    def generate_force_file(self, xyz_file:str):
        with open(xyz_file, 'r') as f:
            cont = f.readlines()

        header = cont[0:1]

        coordinates = np.array([line.strip().split()[1:] for line in cont[2:]], dtype=float)
        symbols = np.array([line.strip().split()[0] for line in cont[2:]], dtype=str)

        x_force = self.x.reshape([-1,3])
        y_force = self.y.reshape([-1,3])

        dx = self.sigma * np.cos(self.theta_s + np.pi)
        dy = self.sigma * np.sin(self.theta_s + np.pi)
        
        min_tilt_force = dx * x_force + dy * y_force

        min_tilt_force /= np.linalg.norm(min_tilt_force)

        # print(y_force)

        with open('vectors.xyz', 'w') as vecfile:
            vecfile.write(header[0])
            vecfile.write('not altered\n')
            for index, coordinate in enumerate(coordinates):
                symbol = symbols[index]
                coord = coordinate
                forces = x_force[index]

                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)}\n')
                # vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces)} 3\n')

        with open('vectors.xyz', 'a') as vecfile:
            vecfile.write(header[0])
            vecfile.write('x vector \n')
            for index, coordinate in enumerate(coordinates):
                symbol = symbols[index]
                coord = coordinates[index]
                forces = x_force[index]
                # normalized_force =  x_force[index] / np.linalg.norm(x_force[index])
                # forces = normalized_force

                # vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord + forces)}\n')
                # vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces/np.linalg.norm(forces))} 3\n')
                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces)}\n')

        with open('vectors.xyz', 'a') as vecfile:
             vecfile.write(header[0])
             vecfile.write('y vector \n')
             # print(coordinates)
             for index, coordinate in enumerate(coordinates):
                 symbol = symbols[index]
                 coord = coordinate
                 forces = y_force[index]
 
                 # vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces/np.linalg.norm(forces))} 3\n')
                 vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces)} \n')
                 # vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord + forces)}\n')

        with open('vectors.xyz', 'a') as vecfile:
             vecfile.write(header[0])
             vecfile.write('min tilt force vector (normalized) \n')
             # print(coordinates)
             for index, coordinate in enumerate(coordinates):
                 symbol = symbols[index]
                 coord = coordinate
                 forces = min_tilt_force[index]
 
                 # vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces/np.linalg.norm(forces))} 3\n')
                 vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces)} \n')
                 # vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord + forces)}\n')

        other_quadrant = - dx * x_force + dy * y_force
        

        other_quadrant /= np.linalg.norm(other_quadrant)

        with open('vectors.xyz', 'a') as vecfile:
             vecfile.write(header[0])
             vecfile.write('Other quadrant for bifurcating purposes force vector (normalized) \n')
             # print(coordinates)
             for index, coordinate in enumerate(coordinates):
                 symbol = symbols[index]
                 coord = coordinate
                 forces = other_quadrant[index]
 
                 # vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces/np.linalg.norm(forces))} 3\n')
                 vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces)} \n')
                 # vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord + forces)}\n')

        # print(np.linalg.norm(x_force).reshape(-1))

    def _pre_plot_2d(self, max_grid:float = 1, surf:str='a'):
        x = np.linspace(-max_grid, max_grid, 501)
        y = np.linspace(-max_grid, max_grid, 501)
        
        X, Y = np.meshgrid(x, y)
        if surf == 'a':
            ener = self.E_A(X, Y)
        else:
            ener = self.E_B(X, Y)

        # plt.figure(figsize=(6, 5))

        # contour lines
        contour_levels = 15  
        contours = plt.contour(
            X, Y, ener,
            levels=contour_levels,
            colors='black',
            linewidths=0.7
        )
        
        # label contour lines
        plt.clabel(contours, inline=True, fontsize=8, fmt="%.2f")
        
        plt.imshow(
            ener,
            extent=(-max_grid, max_grid, -max_grid, max_grid),
            origin='lower',
            cmap='viridis',
            aspect='auto'
        )
        if surf =='a':
            dx = self.sigma * np.cos(self.theta_s + np.pi)
            dy = self.sigma * np.sin(self.theta_s + np.pi)
        else:
            dx = self.sigma * np.cos(self.theta_s)
            dy = self.sigma * np.sin(self.theta_s)

        
        plt.arrow(0, 0, dx, dy, color='red', head_width=0.05, head_length=0.1, length_includes_head=False)
        plt.arrow(0, 0, 1, 0, color='black', head_width=0.05, head_length=0.1, length_includes_head=True)
        plt.arrow(0, 0, 0, 1, color='black', head_width=0.05, head_length=0.1, length_includes_head=True)

        plt.colorbar(label="$\Delta E_{lower}$")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.xticks([-1, 0, 1])
        plt.yticks([-1, 0, 1])
        plt.tight_layout()
    
    def plot_2d(self,surf:str='a'):
        self._pre_plot_2d(surf=surf)
        plt.show()
    
    def saveplot_2d(self, filename, dpi:int=800, surf='a'):
        self._pre_plot_2d(surf=surf)
        plt.savefig(filename, dpi=dpi)

    def displace_geoms(self, directory:str='', xyz_file:str='', rescaling:float=0.1):
        with open(xyz_file, 'r') as f:
            cont = f.readlines()

        header = cont[0:1]

        ref_coordinates = np.array([line.strip().split()[1:] for line in cont[2:]], dtype=float)
        symbols = np.array([line.strip().split()[0] for line in cont[2:]], dtype=str)

        x_force = self.x.reshape([-1,3])
        y_force = self.y.reshape([-1,3])

        dx = self.sigma * np.cos(self.theta_s + np.pi)
        dy = self.sigma * np.sin(self.theta_s + np.pi)
        
        min_tilt_force = dx * x_force + dy * y_force
        min_tilt_force /= np.linalg.norm(min_tilt_force)

        mirrored_force = -dx * x_force + dy * y_force
        mirrored_force /= np.linalg.norm(mirrored_force)

        # print(y_force)

        with open(f'{directory}/{xyz_file.replace(".xyz", "_x.xyz")}', 'w') as vecfile:

            vecfile.write(header[0])
            vecfile.write(f'displaced in x rescaled by {rescaling}\n')
            coordinates = x_force * rescaling + np.copy(ref_coordinates)
            
            for index, coordinate in enumerate(coordinates):
                symbol = symbols[index]
                coord = coordinate

                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)}\n')

        with open(f'{directory}/{xyz_file.replace(".xyz", "_minus_x.xyz")}', 'w') as vecfile:

            vecfile.write(header[0])
            vecfile.write(f'displaced in -x rescaled by {rescaling}\n')
            coordinates = -x_force * rescaling + np.copy(ref_coordinates)
            
            for index, coordinate in enumerate(coordinates):
                symbol = symbols[index]
                coord = coordinate

                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)}\n')

        with open(f'{directory}/{xyz_file.replace(".xyz", "_y.xyz")}', 'w') as vecfile:

            vecfile.write(header[0])
            vecfile.write(f'displaced in y rescaled by {rescaling}\n')
            coordinates = y_force * rescaling + np.copy(ref_coordinates)
            
            for index, coordinate in enumerate(coordinates):
                symbol = symbols[index]
                coord = coordinate

                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)}\n')

        with open(f'{directory}/{xyz_file.replace(".xyz", "_minus_y.xyz")}', 'w') as vecfile:

            vecfile.write(header[0])
            vecfile.write(f'displaced in -y rescaled by {rescaling}\n')
            coordinates = -y_force * rescaling + np.copy(ref_coordinates)
            
            for index, coordinate in enumerate(coordinates):
                symbol = symbols[index]
                coord = coordinate

                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)}\n')

        with open(f'{directory}/{xyz_file.replace(".xyz", "_s_ab.xyz")}', 'w') as vecfile:

            vecfile.write(header[0])
            vecfile.write(f'displaced in -s_ab rescaled by {rescaling}\n')
            coordinates = min_tilt_force * rescaling + np.copy(ref_coordinates)
            
            for index, coordinate in enumerate(coordinates):
                symbol = symbols[index]
                coord = coordinate

                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)}\n')

        with open(f'{directory}/{xyz_file.replace(".xyz", "_mirror_s_ab.xyz")}', 'w') as vecfile:

            vecfile.write(header[0])
            vecfile.write(f'displaced in mirrored (in x) s_ab rescaled by {rescaling}\n')
            coordinates = mirrored_force * rescaling + np.copy(ref_coordinates)
            
            for index, coordinate in enumerate(coordinates):
                symbol = symbols[index]
                coord = coordinate

                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)}\n')

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # j = Javi.from_molcas('ci-trampa-caspt2.input.out')

    args = parser.parse_args()

    # try:
    if args.molcas_file != '0.0':
        j = Javi.from_molcas(args.molcas_file)
    else:
        print('Manually selected files')
        j = Javi(
            args.g0,
            args.g1,
            args.cdv
        )
    # except:
    #     parser.print_help(sys.stderr)
    #     print('\n\n\nIt requires g0, g1 and cdv OR molcas output file')
    #     exit()

    # print(f'scalar    = {j.x@j.y:8.4f}')
    # print(f'beta      = {j.beta%np.pi:8.4f}')
    # print(f'pitch     = {j.pitch:8.4f}')
    # print(f'asymmetry = {j.asymmetry:8.4f}')
    # print(f'tilt      = {j.sigma:8.4f}')
    # print(f'theta     = {j.theta_s/2/np.pi*360:8.1f}')
    # print(f'is_rot_ne = {j.is_rotation_needed}')

    # print(np.array2string(j.x.reshape(-1, 3), precision=6, suppress_small=True))
    # print(np.array2string(j.y.reshape(-1, 3), precision=6, suppress_small=True))
    # print(j.x.reshape([-1,3]))
    # print(j.y.reshape([-1,3]))
    if args.ref_xyz != 'None':
        if args.generate_forces:
            j.generate_force_file(args.ref_xyz)
        if args.generate_displacements:
            os.makedirs('displaced_geoms', exist_ok=True)
            j.displace_geoms(
                    directory='displaced_geoms',
                    xyz_file=args.ref_xyz,
                    rescaling=args.disp_mod
                )
    
    if args.plot_lower:
        j.plot_2d(surf='a')

    if args.plot_upper:
        j.plot_2d(surf='b')

    if args.interactive:
        j.plot_CI()
    
    # TESTING
    # j = Javi(
    #     'engrad_0_gradient.dat', 
    #     'engrad_1_gradient.dat', 
    #     'y_minus_one.dat',
    # )
    # print(j.g_ab)
    # print(j.beta)
    # if args.interactive == True:
    #     j.plot_CI()

