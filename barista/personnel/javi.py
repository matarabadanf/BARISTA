#!/usr/bin/env python3
from functools import cached_property
import argparse
import os
from typing import Tuple
import numpy as np


# Parser is used to input via terminal the required arguments for Emma
parser = argparse.ArgumentParser(
    description="""Takes the ground state geometry of a molecule and compares it to optimized TDDFT geometries run in ORCA.\
        \nIt compares the energy, RMSD and final root, plotting the results in an image. """,
    epilog="""It should run appropriately with ORCA 5.0.X""",
)

parser.add_argument(
    "-fc",
    type=str,
    required=True,
    default=42,
    help="Optimized ground state .xyz file",
)
parser.add_argument(
    "-ex",
    type=str,
    required=True,
    nargs="+",
    default=42,
    help='Excited state optimization .out "filenames" (it is important that the name ends with *_N.in.out)',
)
parser.add_argument(
    "-O",
    type=str,
    default=42,
    help='Output file extension (default=".in.out")',
)

parser.add_argument("-o", type=str, default=42, help="Output image filename")

parser.add_argument(
    "-r", type=str, default=42, help="Generate a report table in output file"
)

parser.add_argument(
    "-md",
    type=str,
    default=True,
    help="Generate report in MD format (default is True)",
)
parser.add_argument(
    "-i",
    type=bool,
    default=False,
    help="Interactive plot mode (default is False)",
)

# display help message if there is no arguments


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

        self._init()

    def _init(self) -> None:
        self._load_vectors()

    @property
    def ci_energy(self) -> float:
        return float(self._ci_energy)

    def _load_vectors(self) -> None:
        self._ga = np.loadtxt(self._ga_filename).reshape(-1)
        self._gb = np.loadtxt(self._gb_filename).reshape(-1)
        self._h_ab = np.loadtxt(self._hab_filename).reshape(-1)

    @cached_property
    def g_ab(self) -> np.ndarray:
        return 0.5 * np.copy(self._ga + self._gb)


    @cached_property
    def s_ab(self) -> np.ndarray:
        return 0.5 * np.copy(self._gb + self._ga)
    
    @cached_property
    def h_ab(self) -> np.ndarray:
        return np.copy(self._h_ab)

    @cached_property
    def beta(self) -> float:
        tan_2_beta = 2 * (np.dot(self.g_ab, self.h_ab)) / (np.linalg.norm(self.g_ab) - np.linalg.norm(self.h_ab))

        return np.atan(tan_2_beta)

    @cached_property
    def _g_tilde(self) -> np.ndarray:
        return self.g_ab * np.cos(self.beta) + self.h_ab - np.sin(self.beta)

    @cached_property
    def _h_tilde(self) -> np.ndarray:
        return self.g_ab * np.cos(self.beta) - self.h_ab - np.sin(self.beta)

    @cached_property
    def x(self) -> np.ndarray:
        return np.copy(self._g_tilde / np.linalg.norm(self._g_tilde))

    @cached_property
    def y(self) -> np.ndarray:
        return np.copy(self._h_tilde / np.linalg.norm(self._h_tilde))

    @cached_property
    def pitch(self) -> float:
        ''' Pitch \\delta_gh. '''

        return np.sqrt( 1/2 * (np.linalg.norm(self._g_tilde) + np.linalg.norm(self._h_tilde)))

    @cached_property
    def asymmetry(self) -> float:
        ''' Asymmetry \\Delta_gh. '''

        asym = (np.linalg.norm(self._g_tilde) - np.linalg.norm(self._h_tilde)) / (np.linalg.norm(self._g_tilde) + np.linalg.norm(self._h_tilde))

        return float(asym)

    def energy_difference(self, x:float, y:float) -> float:
        return 2 * self.pitch * np.sqrt((x**2 + y**2) + self.asymmetry * (x**2 - y**2))
 
    def average_energy(self, x:float, y:float) -> float:
        return self.ci_energy + x * np.dot(self.s_ab, self.x) + y * np.dot(self.s_ab, self.y)

    def E_A(self, x:float, y:float) -> float:
        diff = self.energy_difference(x, y)
        aver = self.average_energy(x,y)
        return aver - diff

    def E_B(self, x:float, y:float) -> float:
        diff = self.energy_difference(x, y)
        aver = self.average_energy(x,y)
        return aver + diff

    @cached_property
    def theta_s(self, max_grid = .1):
        x = np.linspace(-max_grid, max_grid, 101)
        y = np.linspace(-max_grid, max_grid, 101)
        
        X, Y = np.meshgrid(x, y)
        b_p = self.ci_energy + np.zeros_like(X)
        e_mean = self.average_energy(X, Y)
    
    
        max_index = np.unravel_index(np.argmax(e_mean), e_mean.shape)
        
        max_x = X[max_index]
        max_y = Y[max_index]

        vec = np.array([max_x, max_y])
        x = np.array([1,0])

        theta = np.arccos(np.dot(vec, x) / (np.linalg.norm(vec) * np.linalg.norm(x)))
        
        return theta

    @cached_property
    def sigma(self):
        s_x = np.dot(self.s_ab, self.x) / self.pitch
        s_y = np.dot(self.s_ab, self.y) / self.pitch
        
        return np.sqrt(s_x**2 + s_y**2)

    @cached_property
    def p(self) -> Tuple[float, str]:

        p = self.sigma**2 / (1 - self.asymmetry**2) * (1 - self.asymmetry * np.cos(2*self.theta_s))

        p_type = 'Peaked' if p < 1 else 'Sloped'

        return (p, p_type)

    @cached_property
    def b(self) -> Tuple[float, str]:

        b = (self.sigma**2/(4*self.asymmetry**2)) **(1/3) * (((1+self.asymmetry)*np.cos(self.theta_s)**2)**(1/3) + ((1-self.asymmetry)*np.sin(self.theta_s)**2)**(1/3))

        b_type = 'Bifurcating' if b < 1 else 'Single Path'

        return (b, b_type)

    def plot_CI(self, max_grid:float=0.1):
        try:
            import plotly.graph_objects as go
        except:
            print('Plotly is required')
    
        x = np.linspace(-max_grid, max_grid, 101)
        y = np.linspace(-max_grid, max_grid, 101)
        
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
    
        max_index = np.unravel_index(np.argmax(e_mean), e_mean.shape)
        
        # Find the x, y values corresponding to the maximum
        max_x = X[max_index]
        max_y = Y[max_index]
        
        # Find the maximum value
        max_value = e_mean[max_index]
        
        print(f"The maximum value is {max_value} at x = {max_x}, y = {max_y}")
    
        # rgba(100, 120, 140, 0.7)
        # Create a 3D surface plot
        fig = go.Figure(data=[
            go.Surface(z=e_a,    x=X, y=Y, showscale=False, name='Surface 1', opacity=0.8, colorscale=[[0, 'rgba(100, 120, 140, 0.7)'],[1, 'rgba(100, 120, 140, 0.7)']]),
            go.Surface(z=e_b,    x=X, y=Y, showscale=False, name='Surface 2', opacity=0.8, colorscale=[[0, 'rgba(173, 216, 230, 0.7)'],[1, 'rgba(173, 216, 230, 0.7)']]),
            go.Surface(z=e_mean, x=X, y=Y, showscale=False, name='Branching plane', opacity=0.5, colorscale=[[0, 'rgba(255, 0, 0, 0.7)'],[1, 'rgba(255, 0, 0, 0.7)']]),
            go.Surface(z=b_p,    x=X, y=Y, showscale=False, name='Mean energy plane', opacity=0.5, colorscale=[[0, 'rgba(0, 255, 0, 0.7)'],[1, 'rgba(0, 255, 0, 0.7)']]),
    
            # Vectors
            go.Scatter3d(x=x_intersect, y=y_intersect, z=z_intersect, mode='lines+text', line=dict(color='rgba(0, 0, 0, 0.5)', width=4), name='BP and ME intersection'),
            go.Scatter3d(x=[0, max_x], y=[0, max_y], z=[0, max_value], mode='lines+text',line=dict(color='black', width=4), name = r'Max tilt direction', text=['',  'Max tilt direction'], textfont=dict(color='black')),
            go.Scatter3d(x=[0, 1*max_grid], y=[0, 0], z=[0, 0], mode='lines+text',line=dict(color='black', width=4), name = r'$\hat{\mathbf{x}}$', text=['',  'x'], textfont=dict(color='black')),
            go.Scatter3d(x=[0, 0], y=[0, 1*max_grid], z=[0, 0], mode='lines+text',line=dict(color='black', width=4), name = r'$\hat{\mathbf{y}}$', text=['', 'y'], textfont=dict(color='black')),
    
            # Dummy scatter traces for the legend to reflect the surface color
            go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(100, 120, 140, 0.7)', size=10), name='Surface 1'),
            go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(173, 216, 230, 0.7)', size=10), name='Surface 2'),
            go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(255, 0, 0, 0.7)', size=10), name='Branching plane'),
            go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(0, 255, 0, 0.7)', size=10), name='Mean energy plane'),
    
    
        ])
    
        fig.update_layout(
            title='Conical intersection representation',
            title_font=dict(size=28),
            annotations=[
                dict(
                    x=1.05,  # X position of the text box (relative to the plot)
                    y=0.5,  # Y position of the text box (relative to the plot)
                    xref='paper',  # Use 'paper' for positioning relative to the whole figure
                    yref='paper',  # Use 'paper' for positioning relative to the whole figure
                    text=f"P = {self.p[0]:5.3f} -> {self.p[1]} conical.<br>B = {self.b[0]:5.3f} -> {self.b[1]} conical.",  # Text to display in the box
                    showarrow=False,  # Don't show an arrow
                    font=dict(size=16, color="black"),  # Font size and color
                    bordercolor="black",  # Border color of the box
                    borderwidth=2  # Border width
                )
            ],
            scene=dict(
                xaxis_title='X',
                yaxis_title='Y',
                zaxis_title='Energy / Hartree',
                aspectratio=dict(x=1, y=1, z=2)
            ),
            showlegend=True,
            legend=dict(
                font=dict(size=22),  # Item font size
                traceorder='normal',  # Order of legend items
            )
        )
    
        fig.show()

    

if __name__ == "__main__":
    # if len(sys.argv) == 1:
    #     parser.print_help(sys.stderr)
    #     sys.exit(1)

    # args = parser.parse_args()

    
    # TESTING
    j = Javi(
        'engrad_0_gradient.dat', 
        'engrad_1_gradient.dat', 
        'y_minus_one.dat',
    )

    print(j.g_ab)
    print(j.beta)

    j.plot_CI()



