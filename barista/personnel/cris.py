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
    description="""Takes the lower and upper state gradients in a CI with the CDV vector and characterizes the type oc CI.\
        \nResults can be plotted interactively using plotly.""",
    epilog='It requires g0, g1 and cdv OR molcas output file',
)

#############################################
#       Core data for characterization      #
#############################################

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
    action='store_true',
    required=False,
    help="Generate xyz file with  x, y, and -s_ab vectors for visualization."
)

forces_args.add_argument(
    "--generate_displacements",
    action='store_true',
    required=False,
    help="Generate displaced geometries in x, -x, y, -y -s_ab and mirrored -s_ab. A folder will be generated with these geometries.",
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
    """
    Class for analyzing conical intersections (CIs) using gradients and CDV vectors.
    Provides methods for parsing Molcas output, calculating CI properties, and visualization.
    """

    def __init__(
        self,
        ga_filename: str,
        gb_filename: str,
        hab_filename: str,
        ci_energy:float = 0.0
    ) -> None:
        """
        Initialize Javi object with filenames for gradients and CDV vector.
        Args:
            ga_filename: Path to lower state gradient file.
            gb_filename: Path to higher state gradient file.
            hab_filename: Path to CDV vector file.
            ci_energy: CI energy (default 0.0).
        """
        self._ga_filename = ga_filename
        self._gb_filename = gb_filename
        self._hab_filename = hab_filename
        self._ci_energy = ci_energy
        self.rotation_counter = 0

        self._init()

    def _init(self) -> None:
        """
        Load vectors and perform initial rotation checks.
        """
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
        """
        Parse a Molcas output file to extract gradients and CDV vector, then create a Javi instance.
        Args:
            output_filename: Path to Molcas .out file.
        Returns:
            Javi instance initialized with parsed data.
        """
        with open(f'{output_filename}', 'r') as output_file:
            output_list = output_file.readlines()
    
            output_list =[line.strip() for line in output_list]
            
            gradient_pre_indices = [i for i, line in enumerate(output_list) if 'Molecular gradients' in line]
            pt_gradient_pre_indices = [i for i, line in enumerate(output_list) if 'Numerical gradient     ' in line]
    
            if len(gradient_pre_indices) > 2:
                print('More than 2 gradients were calculated, assuming the last two are the highest theory level ones.')

            gradient_indices = []
            gradient_ends = []

            for index in gradient_pre_indices[-2:]:
                lines = []
                for linenumber, line in enumerate(output_list[index:index +20]):
                    if '-----' in line:
                        lines.append(int(linenumber + index))
                
                gradient_indices.append(lines[-2]+1)
                gradient_ends.append(lines[-1])

            cdv_pre_index = [i for i, line in enumerate(output_list) if 'CI derivative coupling' in line][-1]
            lines = []
            for linenumber, line in enumerate(output_list[cdv_pre_index:cdv_pre_index + 20]):
                if '-----' in line:
                    lines.append(int(linenumber + cdv_pre_index))
            
            cdv_start, cdv_end = lines[-2]+1 , lines[-1]

            grad_0 = np.array([
                line.strip().split()[1:] for line in output_list[gradient_indices[0]:gradient_ends[0]]
                ], dtype=float)

            grad_1 = np.array([
                line.strip().split()[1:] for line in output_list[gradient_indices[1]:gradient_ends[1]]
                ], dtype=float)

            cdv = np.array([
                line.strip().split()[1:] for line in output_list[cdv_start:cdv_end]
                ], dtype=float)

            np.savetxt('engrad_0.dat', grad_0)
            np.savetxt('engrad_1.dat', grad_1)
            np.savetxt('cdv.dat', cdv)
    
            return cls('engrad_0.dat', 'engrad_1.dat', 'cdv.dat') 
 


    @property
    def ci_energy(self) -> float:
        return float(self._ci_energy)

    def _load_vectors(self) -> None:
        """
        Load gradient and CDV vectors from files.
        """
        self._ga = np.loadtxt(self._ga_filename).reshape(-1)
        self._gb = np.loadtxt(self._gb_filename).reshape(-1)
        self._h_ab = np.loadtxt(self._hab_filename).reshape(-1)

    @cached_property
    def g_ab(self) -> np.ndarray:
        """
        Branching space vector g_ab.
        Returns:
            np.ndarray: g_ab vector.
        """
        return 0.5 * np.copy(self._gb - self._ga)

    @cached_property
    def s_ab(self) -> np.ndarray:
        """
        Average gradient vector s_ab.
        Returns:
            np.ndarray: s_ab vector.
        """
        return 0.5 * np.copy(self._gb + self._ga)
    
    @cached_property
    def h_ab(self) -> np.ndarray:
        """
        CDV vector h_ab.
        Returns:
            np.ndarray: h_ab vector.
        """
        return np.copy(self._h_ab)

    @cached_property
    def _pre_beta(self) -> np.array:
        """
        Calculate possible beta angles for rotation.
        Returns:
            list: Four possible beta values.
        """
        beta = 0.5 * np.arctan2(2 * np.dot(self.g_ab, self.h_ab), (np.dot(self.g_ab, self.g_ab) - np.dot(self.h_ab, self.h_ab)))


        beta_2 = np.arctan2(
            2 * np.dot(self.g_ab, self.h_ab),
            np.dot(self.g_ab, self.g_ab) - np.dot(self.h_ab, self.h_ab)
        )
        beta = beta_2/2

        return [beta + i*np.pi/2 for i in range(0, 4)]
        # return 0.5 * np.arctan2(2 * np.dot(self.g_ab, self.h_ab), (np.dot(self.g_ab, self.g_ab) - np.dot(self.h_ab, self.h_ab)))

    def _rotate_for_beta(self):
        """
        Increment rotation counter for beta angle.
        """
        self.rotation_counter += 1 

    @property
    def beta(self) -> float:
        """
        Get current beta angle after rotation.
        Returns:
            float: Beta angle.
        """
        return self._pre_beta[self.rotation_counter % 4] # 2 works for x and y unit vector sign in ethylene a and b  

    def is_rotation_needed(self, asymmetry:float, theta_s:float) -> bool:
        """
        Check if rotation is needed based on asymmetry and theta_s.
        Args:
            asymmetry: Asymmetry value.
            theta_s: Theta_s value.
        Returns:
            bool: True if rotation is needed, False otherwise.
        """
        return not (asymmetry > 0 and theta_s > -10**-12  and theta_s < np.pi/2)

        # return not (self.asymmetry > 0 and self.theta_s > -10**-12  and self.theta_s < np.pi/2)
        

    @property
    def _g_tilde(self) -> np.ndarray:
        """
        Rotated g_ab vector.
        Returns:
            np.ndarray: Rotated g_ab.
        """
        return self.g_ab * np.cos(self.beta) + self.h_ab * np.sin(self.beta)

    @property
    def _h_tilde(self) -> np.ndarray:
        """
        Rotated h_ab vector.
        Returns:
            np.ndarray: Rotated h_ab.
        """
        return self.h_ab * np.cos(self.beta) - self.g_ab * np.sin(self.beta) # original

    @property
    def x(self) -> np.ndarray:
        """
        Normalized x vector in branching plane.
        Returns:
            np.ndarray: x vector.
        """
        return np.copy(self._g_tilde / np.linalg.norm(self._g_tilde))

    @property
    def y(self) -> np.ndarray:
        """
        Normalized y vector in branching plane.
        Returns:
            np.ndarray: y vector.
        """
        return np.copy(self._h_tilde / np.linalg.norm(self._h_tilde))

    @property
    def pitch(self) -> float:
        """
        Pitch (delta_gh) of the CI.
        Returns:
            float: Pitch value.
        """
        return np.sqrt( 1/2 * (np.dot(self._g_tilde, self._g_tilde) + np.dot(self._h_tilde, self._h_tilde)))

    @property
    def asymmetry(self) -> float:
        """
        Asymmetry (Delta_gh) of the CI.
        Returns:
            float: Asymmetry value.
        """
        asym = (np.dot(self._g_tilde, self._g_tilde) - np.dot(self._h_tilde, self._h_tilde)) / (np.dot(self._g_tilde, self._g_tilde) + np.dot(self._h_tilde, self._h_tilde))

        # return abs(float(asym))
        return float(asym) # original

    def energy_difference(self, x:float, y:float) -> float:
        """
        Calculate energy difference between CI surfaces at (x, y).
        Args:
            x: X coordinate.
            y: Y coordinate.
        Returns:
            float: Energy difference.
        """
        return 2 * self.pitch * np.sqrt((x**2 + y**2) + self.asymmetry * (x**2 - y**2))
 
    def average_energy(self, x:float, y:float) -> float:
        """
        Calculate average energy at (x, y).
        Args:
            x: X coordinate.
            y: Y coordinate.
        Returns:
            float: Average energy.
        """
        return self.ci_energy + x * np.dot(self.s_ab, self.x) + y * np.dot(self.s_ab, self.y)

    def E_A(self, x:float, y:float) -> float:
        """
        Calculate lower state energy at (x, y).
        Args:
            x: X coordinate.
            y: Y coordinate.
        Returns:
            float: Lower state energy.
        """
        s_x = np.dot(self.s_ab, self.x) / self.pitch 
        s_y = np.dot(self.s_ab, self.y) / self.pitch
        return self.ci_energy + self.pitch * (x * s_x + y * s_y - ((x**2 + y **2) + self.asymmetry * (x**2 -y**2))**0.5)


        diff = self.energy_difference(x, y)
        aver = self.average_energy(x,y)
        return aver - diff

    def E_B(self, x:float, y:float) -> float:
        """
        Calculate upper state energy at (x, y).
        Args:
            x: X coordinate.
            y: Y coordinate.
        Returns:
            float: Upper state energy.
        """
        diff = self.energy_difference(x, y)
        aver = self.average_energy(x,y)
        return aver + diff

    @property
    def theta_s(self, n_points:int = 1800, radius:float=1):
        """
        Calculate the angle of maximum tilt in the branching plane.
        Args:
            n_points: Number of points to sample (default 1800).
            radius: Radius for sampling (default 1).
        Returns:
            float: Angle theta_s in radians.
        """
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

    @property
    def sigma(self):
        """
        Calculate the tilt (sigma) of the CI.
        Returns:
            float: Sigma value.
        """
        s_x = np.dot(self.s_ab, self.x) / self.pitch
        s_y = np.dot(self.s_ab, self.y) / self.pitch
        
        return np.sqrt(s_x**2 + s_y**2)

    @property
    def p(self) -> Tuple[float, str]:

        """
        Calculate the P parameter and its type (Peaked/Sloped).
        Returns:
            Tuple[float, str]: (P value, type string)
        """
        p = self.sigma**2 / (1 - self.asymmetry**2) * (1 - self.asymmetry * np.cos(2*self.theta_s))

        p_type = 'Peaked' if p < 1 else 'Sloped'

        return (p, p_type)

    @property
    def b(self) -> Tuple[float, str]:

        """
        Calculate the B parameter and its type (Bifurcating/Single Path).
        Returns:
            Tuple[float, str]: (B value, type string)
        """
        b = (self.sigma**2/(4*self.asymmetry**2)) **(1/3) * (((1+self.asymmetry)*np.cos(self.theta_s)**2)**(1/3) + ((1-self.asymmetry)*np.sin(self.theta_s)**2)**(1/3))

        b_type = 'Bifurcating' if b < 1 else 'Single Path'

        return (b, b_type)

    def plot_CI(self, max_grid:float=1, filename:str = 'NONE.html'):
        """
        Plot the conical intersection surfaces interactively using Plotly.
        Args:
            max_grid: Range for x and y axes.
            filename: Output HTML file for the plot.
        """
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
        """
        Generate an XYZ file with force vectors for visualization.
        Args:
            xyz_file: Path to reference XYZ file.
        """
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

                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces)}\n')

        with open('vectors.xyz', 'a') as vecfile:
            vecfile.write(header[0])
            vecfile.write('y vector \n')
            for index, coordinate in enumerate(coordinates):
                symbol = symbols[index]
                coord = coordinate
                forces = y_force[index]
 
                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces)} \n')

        with open('vectors.xyz', 'a') as vecfile:
            vecfile.write(header[0])
            vecfile.write('min tilt force vector (normalized) \n')
            for index, coordinate in enumerate(coordinates):
                symbol = symbols[index]
                coord = coordinate
                forces = min_tilt_force[index]
 
                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces)} \n')

        other_quadrant = - dx * x_force + dy * y_force
        other_quadrant /= np.linalg.norm(other_quadrant)

        with open('vectors.xyz', 'a') as vecfile:
            vecfile.write(header[0])
            vecfile.write('Other quadrant for bifurcating purposes force vector (normalized) \n')
            for index, coordinate in enumerate(coordinates):
                symbol = symbols[index]
                coord = coordinate
                forces = other_quadrant[index]
 
                vecfile.write(f'{symbol} {" ".join(f"{x:12.8f}" for x in coord)} {" ".join(f"{f:12.8f}" for f in forces)} \n')

    def _pre_plot_2d(self, max_grid:float = 1, surf:str='a'):
        """
        Prepare a 2D contour plot of the CI surface.
        Args:
            max_grid: Range for x and y axes.
            surf: Surface to plot ('a' for lower, 'b' for upper).
        """
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
        """
        Show a 2D contour plot of the CI surface.
        Args:
            surf: Surface to plot ('a' for lower, 'b' for upper).
        """
        self._pre_plot_2d(surf=surf)
        plt.show()
    
    def saveplot_2d(self, filename, dpi:int=800, surf='a'):
        """
        Save a 2D contour plot of the CI surface to a file.
        Args:
            filename: Output file name.
            dpi: Dots per inch for the saved figure.
            surf: Surface to plot ('a' for lower, 'b' for upper).
        """
        self._pre_plot_2d(surf=surf)
        plt.savefig(filename, dpi=dpi)

    def displace_geoms(self, directory:str='', xyz_file:str='', rescaling:float=0.1):
        """
        Generate displaced geometries along branching plane vectors and save as XYZ files.
        Args:
            directory: Output directory for displaced geometries.
            xyz_file: Reference XYZ file.
            rescaling: Displacement magnitude.
        """
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

    def plot_polar(self, n_points:int=360):
        
        n_points = n_points

        theta = np.linspace(0, 2*np.pi, n_points)
        x = np.array([np.cos(t) for t in theta])
        y = np.array([np.sin(t) for t in theta])

        coordinate_pairs = np.zeros([n_points,2])
        print(coordinate_pairs)
        for i in range(n_points):
            coordinate_pairs[i] = x[i],y[i] 

        print(coordinate_pairs)

        polar_energy = [self.average_energy(x,y) for x,y in coordinate_pairs]
        a_energy = [self.E_A(x,y) for x,y in coordinate_pairs]
        b_energy = [self.E_B(x,y) for x,y in coordinate_pairs]

        print(polar_energy)

        plt.plot(theta, polar_energy)
        plt.plot(theta, a_energy)
        plt.plot(theta, b_energy)
        plt.plot(theta, np.zeros_like(theta))
        xticks = [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi]

        xtick_labels = ['x', 'y', '-x', '-y', 'x']
        plt.xticks(xticks, xtick_labels)


        plt.axvline(x=self.theta_s +np.pi, color='red', linestyle='--', label=r'$\theta = \frac{\pi}{2}$')
        plt.axvline(x=self.theta_s, color='red', linestyle='--', label=r'$\theta = \frac{\pi}{2}$')

        plt.show()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # j = Javi.from_molcas('ci-trampa-caspt2.input.out')

    args = parser.parse_args()

    try:
        if args.molcas_file != '0.0':
            j = Javi.from_molcas(args.molcas_file)
        else:
            print('Manually selected files')
            j = Javi(
                args.g0,
                args.g1,
                args.cdv
            )
    except FileNotFoundError as e:
        parser.print_help(sys.stderr)
        print( '\n###################################')
        print(f'# File {e.filename} was not found #')
        print( '###################################\n')
        exit()
    
    except IndexError:
        parser.print_help(sys.stderr)
        print('\n###############################################')
        print('# Something went wrong. Check that g0, g1 and #')
        print('# h_ab are provided in the out file or in the #')
        print('# data files.                                 #')
        print('###############################################\n')
        exit()

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
   
    j.plot_polar(n_points=360)
