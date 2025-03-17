#!/usr/bin/env python3
#from functools import property
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
)

parser.add_argument(
    "-g0",
    type=str,
    required=True,
    help="Lower state gradient file.",
)

parser.add_argument(
    "-g1",
    type=str,
    required=True,
    help="Higher state gradient file.",
)

parser.add_argument(
    "-nac",
    type=str,
    required=True,
    help="NAC vector file.",
)

parser.add_argument(
    "--interactive",
    action='store_true',
    required=False,
    help="Open interactive plot. Default is False.",
)

parser.add_argument(
    "--savefig",
    type=str,
    required=False,
    default='Default',
    help="Save CI figure with selected name.",
)

parser.add_argument(
    "--full_plot",
    action='store_true',
    required=False,
    help="Plot the full Figure, including vectors and directions. Default is False.",
)

class Javi:
    """
    A class to represent a conical intersection (CI) and perform various calculations related to it.
    Attributes
    ----------
    ga_filename : str
        Filename for the alpha component of the gradient vector.
    gb_filename : str
        Filename for the beta component of the gradient vector.
    hab_filename : str
        Filename for the Hamiltonian matrix elements.
    ci_energy : float
        The CI energy value.
    Methods
    -------
    from_xy(cls):
        Class method to create an instance from x and y coordinates.
    ci_energy() -> float:
        Returns the CI energy value.
    g_ab() -> np.ndarray:
        Returns the difference vector between gb and ga.
    s_ab() -> np.ndarray:
        Returns the sum vector of gb and ga.
    h_ab() -> np.ndarray:
        Returns the Hamiltonian matrix elements.
    beta() -> float:
        Returns the beta angle.
    beta(beta: float):
        Sets the beta angle.
    g_tilde() -> np.ndarray:
        Returns the g tilde vector.
    h_tilde() -> np.ndarray:
        Returns the h tilde vector.
    x() -> np.ndarray:
        Returns the normalized g tilde vector.
    y() -> np.ndarray:
        Returns the normalized h tilde vector.
    pitch() -> float:
        Returns the pitch (δ_gh).
    asymmetry() -> float:
        Returns the asymmetry (Δ_gh).
    energy_difference(x: float, y: float) -> float:
        Calculates the energy difference for given x and y coordinates.
    average_energy(x: float, y: float) -> float:
        Calculates the average energy for given x and y coordinates.
    E_A(x: float, y: float) -> float:
        Calculates the energy E_A for given x and y coordinates.
    E_B(x: float, y: float) -> float:
        Calculates the energy E_B for given x and y coordinates.
    theta_s(n_points: int = 360, radius: float = 1) -> float:
        Calculates the angle theta_s for the maximum tilt direction.
    sigma() -> float:
        Returns the sigma value.
    p() -> Tuple[float, str]:
        Returns the p value and its type.
    b() -> Tuple[float, str]:
        Returns the b value and its type.
    plot_CI(max_grid: float = 1):
        Plots the conical intersection representation.
    """

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
        self._calculate_beta()

   
    @classmethod
    def from_xy(cls):
        pass 

    @property
    def ci_energy(self) -> float:
        return float(self._ci_energy)

    def normalize(self, vector: np.ndarray) -> np.ndarray:
        return vector/np.linalg.norm(vector) 

    def _load_vectors(self) -> None:
        self._ga = np.loadtxt(self._ga_filename).flatten()
        self._gb = np.loadtxt(self._gb_filename).flatten()
        self._h_ab = np.loadtxt(self._hab_filename).flatten()

    @property
    def g_ab(self) -> np.ndarray:
        #return self.normalize(0.5 * np.copy(self._gb - self._ga))
        return 0.5 * np.copy(self._gb - self._ga)
    
    @property
    def s_ab(self) -> np.ndarray:
        #return self.normalize(0.5 * np.copy(self._gb + self._ga))
        return 0.5 * np.copy(self._gb + self._ga)
    
    @property
    def h_ab(self) -> np.ndarray:
        #return self.normalize(np.copy(self._h_ab))
        return np.copy(self._h_ab)

    def _calculate_beta(self) -> None:
        # tan_2_beta = 2 * (np.dot(self.g_ab, self.h_ab)) / (np.linalg.norm(self.g_ab) - np.linalg.norm(self.h_ab))
        # return np.arctan(tan_2_beta)
        tan_2_beta = 2 * (np.dot(self.g_ab, self.h_ab)) / (np.dot(self.g_ab, self.g_ab) - np.dot(self.h_ab, self.h_ab))

        beta_value =  0.5 * np.arctan2(2 * np.dot(self.g_ab, self.h_ab), (np.dot(self.g_ab, self.g_ab) - np.dot(self.h_ab, self.h_ab)))

        self._beta = beta_value 

    @property
    def beta(self):
        return self._beta

    @beta.setter
    def beta(self, beta):
        self._beta = beta


    @property
    def g_tilde(self) -> np.ndarray:
        return self.g_ab * np.cos(self.beta) + self.h_ab * np.sin(self.beta)

    @property
    def h_tilde(self) -> np.ndarray:
        return self.h_ab * np.cos(self.beta) - self.g_ab * np.sin(self.beta)

    @property
    def x(self) -> np.ndarray:
        return np.copy(self.g_tilde / np.linalg.norm(self.g_tilde))

    @property
    def y(self) -> np.ndarray:
        return np.copy(self.h_tilde / np.linalg.norm(self.h_tilde))

    @property
    def pitch(self) -> float:
        ''' Pitch \\delta_gh. '''

        return np.sqrt( 1/2 * (np.dot(self.g_tilde, self.g_tilde) + np.dot(self.h_tilde, self.h_tilde)))

    @property
    def asymmetry(self) -> float:
        ''' Asymmetry \\Delta_gh. '''

        asym = (np.dot(self.g_tilde, self.g_tilde) - np.dot(self.h_tilde, self.h_tilde)) / (np.dot(self.g_tilde, self.g_tilde) + np.dot(self.h_tilde, self.h_tilde))

        return float(asym)

    def energy_difference(self, x:float, y:float) -> float:
        return 2 * self.pitch * np.sqrt((x**2 + y**2) + self.asymmetry * (x**2 - y**2))
 
    def average_energy(self, x:float, y:float) -> float:
        # return x + y  
        return self.ci_energy + x * np.dot(self.s_ab, self.x) + y * np.dot(self.s_ab, self.y)

    def E_A(self, x:float, y:float) -> float:
        diff = self.energy_difference(x, y)
        aver = self.average_energy(x,y)
        return aver - diff

    def E_B(self, x:float, y:float) -> float:
        diff = self.energy_difference(x, y)
        aver = self.average_energy(x,y)
        return aver + diff

    @property
    def theta_s(self, n_points:int = 360, radius:float=1):

        angles = np.linspace(0, 2 * np.pi, n_points, endpoint=False)
        x_points = radius * np.cos(angles)
        y_points = radius * np.sin(angles)

        pairs = [[x_points[i], y_points[i]] for i in range(n_points)]

        samples = []

        for pair in pairs:
            x, y = pair
            sample = [x, y, self.average_energy(x, y)]
            #print(f'\n\n\n The xyz values of tilt is: {x:5.2f}, {y:5.2f}, {self.average_energy(x,y):5.2f}')
            samples.append(sample)

        x,y,z = max(samples, key=lambda item: item[2])
        #print(f'\n\n\n The xyz values of the max tilt is: {x:5.2f}, {y:5.2f}, {z:5.2f}')

        vec = np.array([x,y])

        #print(vec)
        cos_theta = np.dot(vec, np.array([1,0])) / np.linalg.norm(vec)
        cos_theta = np.clip(cos_theta, -1.0, 1.0)

        theta = np.arccos(cos_theta)
 
#        theta = -theta if vec[0] < 0 else theta

        #print(f'The angle theta is {theta} or {(theta+np.pi/2)%np.pi}')
        

        #print(f'The angle of the maximum tilt is {theta:5.3} Radians or {theta*360/np.pi:3.5} degrees')
        
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

    def _preplot_CI(self, max_grid:float=1, full:bool= False) -> go.Figure:
    
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

        # Locate the direction of the tilt 


        # max_x = max_grid * np.cos(self.theta_s)
        # max_y = max_grid * np.sin(self.theta_s)
        # max_value = self.average_energy(max_x,max_y)
    
        # rgba(100, 120, 140, 0.7)
        # Create a 3D surface plot
        fig = go.Figure(data=[
            go.Surface(z=e_a,    x=X, y=Y, showscale=False, name='Surface 1', opacity=0.8, colorscale=[[0, 'rgba(100, 120, 140, 0.7)'],[1, 'rgba(100, 120, 140, 0.7)']]),
            go.Surface(z=e_b,    x=X, y=Y, showscale=False, name='Surface 2', opacity=0.8, colorscale=[[0, 'rgba(173, 216, 230, 0.7)'],[1, 'rgba(173, 216, 230, 0.7)']]),
            go.Surface(z=e_mean, x=X, y=Y, showscale=False, name='Branching plane', opacity=0.5, colorscale=[[0, 'rgba(255, 0, 0, 0.7)'],[1, 'rgba(255, 0, 0, 0.7)']]),
            go.Surface(z=b_p,    x=X, y=Y, showscale=False, name='Mean energy plane', opacity=0.5, colorscale=[[0, 'rgba(0, 255, 0, 0.7)'],[1, 'rgba(0, 255, 0, 0.7)']]),

    
            # Dummy scatter traces for the legend to reflect the surface color
            go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(100, 120, 140, 0.7)', size=10), name='Surface 1'),
            go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(173, 216, 230, 0.7)', size=10), name='Surface 2'),
            go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(255, 0, 0, 0.7)', size=10), name='Branching plane'),
            go.Scatter3d(x=[None], y=[None], z=[None], mode='markers', marker=dict(color='rgba(0, 255, 0, 0.7)', size=10), name='Mean energy plane'),
    
    
        ])
    
        if full:
            fig.add_traces([
                # Vectors
                go.Scatter3d(x=x_intersect, y=y_intersect, z=z_intersect, mode='lines+text', line=dict(color='rgba(0, 0, 0, 0.5)', width=4), name='BP and ME intersection'),
                # go.Scatter3d(x=[0, max_x], y=[0, max_y], z=[0, max_value], mode='lines+text',line=dict(color='black', width=4), name = r'Max tilt direction', text=['',  'Max tilt direction'], textfont=dict(color='black')),
                # go.Scatter3d(x=[0, x_max], y=[0, y_max], z=[0, z_max], mode='lines+text',line=dict(color='black', width=4), name = r'Max tilt direction 2', text=['',  'Max tilt direction 2'], textfont=dict(color='black')),
                go.Scatter3d(x=[0, 1*max_grid], y=[0, 0], z=[0, 0], mode='lines+text',line=dict(color='black', width=4), name = r'$\hat{\mathbf{x}}$', text=['', 'x'], textfont=dict(color='black')),
                go.Scatter3d(x=[0, 0], y=[0, 1*max_grid], z=[0, 0], mode='lines+text',line=dict(color='black', width=4), name = r'$\hat{\mathbf{y}}$', text=['', 'y'], textfont=dict(color='black')),
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
        return fig

    def plot_CI(self, full:bool=True):
        fig = self._preplot_CI(full=full)
        fig.show()

    def save_CI_fig(self, name:str='CI_plot.jpg', full:bool=True) -> None:
        fig = self._preplot_CI(full=full)
        fig.write_image(name, width=800, height=600)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    j = Javi(
        args.g0,
        args.g1,
        args.nac
    )

    print(j.x)
    print(j.y)

    if args.interactive:
        j.plot_CI()
    
    if args.savefig != 'Default':
        j.save_CI_fig(name=args.savefig, full=args.full_plot)
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


