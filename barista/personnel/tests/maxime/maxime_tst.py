from barista.personnel.maxime import Maxime
import os 

script_dir = os.path.dirname(os.path.abspath(__file__))

xyz_filename = os.path.join(script_dir, 'test.xyz')


m = Maxime(xyz_filename)

print(m.xyz_df.iloc[[1,2]])

m.extract_radius(r=1.1, new_filename='result.xyz')

