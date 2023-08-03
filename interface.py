import numpy as np
import functions
import sys
import configparser

config = configparser.ConfigParser()
config.read(sys.argv[1])

# Extract lattice vectors
lattice_vectors = functions.extract_lattice_vectors('clean_interface_files/bottom_slab.txt')
a = lattice_vectors['a']
b = lattice_vectors['b']
c = lattice_vectors['c']

# Extract atomic coordinates of bottom and upper slabs
atomic_coord_bottom_slab = functions.extract_atomic_coordinates('clean_interface_files/bottom_slab.txt')
atomic_coord_upper_slab =  functions.extract_atomic_coordinates('clean_interface_files/upper_slab.txt')

# Convert direct coordinates to Cartesian coordinates
cartesian_coord_bottom_slab =  functions.direct_to_cartesian_coord(a, b, c, atomic_coord_bottom_slab)
cartesian_coord_upper_slab =  functions.direct_to_cartesian_coord(a, b, c, atomic_coord_upper_slab)

selected_site_Cu = config.get('settings', 'selected_site_Cu')
selected_site_C = config.get('settings', 'selected_site_C')
reference_site_Cu = functions.metal_fcc_111_high_symmetry_points("clean_interface_files/upper_slab.txt",selected_site_Cu)
reference_site_C = functions.C_111_high_symmetry_points("clean_interface_files/bottom_slab.txt", selected_site_C)
# Shift upper slab on xy plane so that the selected high symmetry points of upper and bottom slabs are aligned
shifted_upper_slab_on_xy = functions.shift_slab_on_xy("clean_interface_files/upper_slab.txt", reference_site_Cu,reference_site_C)

# Reflect atomic coordinates of the upper slab with respect to its highest atom along z direction
z = max(cartesian_coord_upper_slab, key=lambda x: x[2])[2]
plane_point = np.array([0, 0,  z])
plane_normal = np.array([0, 0, 1])
reflected_coords =  functions.reflect_coord(shifted_upper_slab_on_xy, plane_point, plane_normal)

# The distance between the two slabs is given by the average value of their respective interlayer distances
interlayer_distance_bottom_slab = float(config.get('settings', 'interlayer_distance_bottom_slab'))
interlayer_distance_upper_slab = float(config.get('settings', 'interlayer_distance_upper_slab'))
slabs_distance = (interlayer_distance_bottom_slab+interlayer_distance_upper_slab)/2
shift_z = max(cartesian_coord_bottom_slab, key=lambda x: x[2])[2]+slabs_distance-max(cartesian_coord_upper_slab, key=lambda x: x[2])[2]
shifted_coords =  functions.shift_slab_along_z(reflected_coords,shift_z)

x_relax = config.get('settings', 'x_relax')
y_relax = config.get('settings', 'y_relax')
z_relax = config.get('settings', 'z_relax')

functions.write_POSCAR_interface('clean_interface_files/upper_slab.txt', 'clean_interface_files/bottom_slab.txt', cartesian_coord_bottom_slab, x_relax, y_relax, z_relax, shifted_coords, a, b, c)

