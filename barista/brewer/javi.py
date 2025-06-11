#!/usr/bin/env python3
from functools import cached_property
import argparse
import os
from typing import Tuple
import numpy as np
import sys

# Parser is used to input via terminal the required arguments for Javi
parser = argparse.ArgumentParser(
    description="""Takes the lower and upper state gradients in a CI with the CDV vector and characterizes the type oc CI.""",
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

