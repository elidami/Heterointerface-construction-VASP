import functions
import sys
import configparser

config = configparser.ConfigParser()
config.read(sys.argv[1])

bottom_slab= config.get('paths', 'clean_bottom_slab')
upper_slab = config.get('paths', 'clean_upper_slab')

# Extract lattice vectors
lattice_vectors = functions.extract_lattice_vectors(bottom_slab)
a = lattice_vectors['a']
b = lattice_vectors['b']
c = lattice_vectors['c']

# Extract atomic coordinates of bottom and upper slabs
atomic_coord_bottom_slab = functions.extract_atomic_coordinates(bottom_slab)
atomic_coord_upper_slab =  functions.extract_atomic_coordinates(upper_slab)

# Convert direct coordinates to Cartesian coordinates
cartesian_coord_bottom_slab =  functions.direct_to_cartesian_coord(a, b, c, atomic_coord_bottom_slab)
cartesian_coord_upper_slab =  functions.direct_to_cartesian_coord(a, b, c, atomic_coord_upper_slab)

# Shift upper slab on xy plane so that the selected high symmetry points of upper and bottom slabs are aligned
selected_site_metal = config.get('settings', 'selected_site_metal')
selected_site_C = config.get('settings', 'selected_site_C')
reference_site_metal = functions.metal_fcc_111_high_symmetry_points(upper_slab,selected_site_metal)
reference_site_C = functions.C_111_high_symmetry_points(bottom_slab, selected_site_C)
shifted_upper_slab_on_xy = functions.shift_slab_on_xy(upper_slab, reference_site_metal,reference_site_C)

# Reflect atomic coordinates of the upper slab with respect to its highest atom along z direction
plane_point, plane_normal = functions.calculate_plane_parameters(cartesian_coord_upper_slab)
reflected_coords =  functions.reflect_coord(shifted_upper_slab_on_xy, plane_point, plane_normal)

# The distance between the two slabs is given by the average value of their respective interlayer distances
interlayer_distance_bottom_slab = float(config.get('settings', 'interlayer_distance_bottom_slab'))
interlayer_distance_upper_slab = float(config.get('settings', 'interlayer_distance_upper_slab'))
shift_z = functions.calculate_shift_z_clean_case(interlayer_distance_bottom_slab, interlayer_distance_upper_slab,cartesian_coord_bottom_slab, cartesian_coord_upper_slab)
shifted_coords =  functions.shift_slab_along_z(reflected_coords,shift_z)

x_relax = config.get('settings', 'x_relax')
y_relax = config.get('settings', 'y_relax')
z_relax = config.get('settings', 'z_relax')

functions.write_POSCAR_interface(upper_slab, bottom_slab, cartesian_coord_bottom_slab, x_relax, y_relax, z_relax, shifted_coords, a, b, c)

