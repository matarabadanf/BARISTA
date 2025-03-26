from barista.personnel.maxime import Maxime
import os 

script_dir = os.path.dirname(os.path.abspath(__file__))

xyz_filename = os.path.join(script_dir, 'test.xyz')


m = Maxime(xyz_filename)

m.extract_radius_atoms(r=10.1, new_filename='kxcmvbnmkolcvbkncvbn.xyz')

