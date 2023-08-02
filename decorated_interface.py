import numpy as np
import functions
import sys
import configparser

config = configparser.ConfigParser()
config.read(sys.argv[1])

# Extract lattice vectors
lattice_vectors =  functions.extract_lattice_vectors('decorated_interface_files/bottom_slab_with_adatom.txt')
a = lattice_vectors['a']
b = lattice_vectors['b']
c = lattice_vectors['c']

# Extract atomic coordinates
atomic_coord_bottom_slab =  functions.extract_atomic_coordinates('decorated_interface_files/bottom_slab_with_adatom.txt')
atomic_coord_upper_slab =  functions.extract_atomic_coordinates('decorated_interface_files/upper_slab.txt')

# Convert direct coordinates to Cartesian coordinates
cartesian_coord_bottom_slab =  functions.direct_to_cartesian_coord(a, b, c, atomic_coord_bottom_slab)
cartesian_coord_upper_slab =  functions.direct_to_cartesian_coord(a, b, c, atomic_coord_upper_slab)


lattice_vectors_adsorption =  functions.extract_lattice_vectors('decorated_interface_files/adsorption_on_upper_slab.txt')
a_adsorption = lattice_vectors_adsorption['a']
b_adosrption = lattice_vectors_adsorption['b']
c_adsorption = lattice_vectors_adsorption['c']

# Extract atomic coordinates
atomic_coord_adsorption =  functions.extract_atomic_coordinates('decorated_interface_files/adsorption_on_upper_slab.txt')

file_path = "decorated_interface_files/upper_slab.txt"
selected_site_Cu = config.get('settings', 'selected_site_Cu')

reference_site_Cu = functions.metal_fcc_111_high_symmetry_points("decorated_interface_files/upper_slab.txt",selected_site_Cu)
print(reference_site_Cu)
upper_slab_coords_for_adatom_adsorption = functions.shift_slab_on_xy(file_path, reference_site_Cu,cartesian_coord_bottom_slab[0])


z = upper_slab_coords_for_adatom_adsorption[-1][2] if upper_slab_coords_for_adatom_adsorption else None
plane_point = np.array([0, 0,  z])
plane_normal = np.array([0, 0, 1])
reflected_coords =  functions.reflect_coord(upper_slab_coords_for_adatom_adsorption, plane_point, plane_normal)
# Convert direct coordinates to Cartesian coordinates
cartesian_coord_adsorption =  functions.direct_to_cartesian_coord(a_adsorption, b_adosrption, c_adsorption, atomic_coord_adsorption)
shift_z = cartesian_coord_bottom_slab[0][2]+functions.distance_between_highest_z_values(cartesian_coord_adsorption)-cartesian_coord_upper_slab[-1][2]
shifted_coords =  functions.shift_slab_along_z(reflected_coords,shift_z)

x_relax = config.get('settings', 'x_relax')
y_relax = config.get('settings', 'y_relax')
z_relax = config.get('settings', 'z_relax')


functions.write_POSCAR_interface('decorated_interface_files/upper_slab.txt', 'decorated_interface_files/bottom_slab_with_adatom.txt', cartesian_coord_bottom_slab, x_relax, y_relax, z_relax, shifted_coords, a, b, c)













































