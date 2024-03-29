import functions
import sys
import configparser

config = configparser.ConfigParser()
config.read(sys.argv[1])

bottom_slab_with_adatom = config.get('paths', 'bottom_slab_with_adatom')
upper_slab = config.get('paths', 'decorated_upper_slab')
adsorption_on_upper_slab = config.get('paths', 'adsorption_on_upper_slab')

# Extract lattice vectors for interface supercell
lattice_vectors =  functions.extract_lattice_vectors(bottom_slab_with_adatom)
a = lattice_vectors['a']
b = lattice_vectors['b']
c = lattice_vectors['c']

# Extract lattice vectors for supercell that regards adatom adsorption on upper slab
lattice_vectors_adsorption =  functions.extract_lattice_vectors(adsorption_on_upper_slab)
a_adsorption = lattice_vectors_adsorption['a']
b_adosrption = lattice_vectors_adsorption['b']
c_adsorption = lattice_vectors_adsorption['c']

# Extract atomic coordinates of bottom and upper slabs
atomic_coord_bottom_slab =  functions.extract_atomic_coordinates(bottom_slab_with_adatom)
atomic_coord_upper_slab =  functions.extract_atomic_coordinates(upper_slab)
atomic_coord_adsorption =  functions.extract_atomic_coordinates(adsorption_on_upper_slab)

# Convert direct coordinates to Cartesian coordinates
cartesian_coord_bottom_slab =  functions.direct_to_cartesian_coord(a, b, c, atomic_coord_bottom_slab)
cartesian_coord_upper_slab =  functions.direct_to_cartesian_coord(a, b, c, atomic_coord_upper_slab)
cartesian_coord_adsorption =  functions.direct_to_cartesian_coord(a_adsorption, b_adosrption, c_adsorption, atomic_coord_adsorption)

# Shift upper slab on xy plane so that its selected high symmetry point is aligned with the adatom adsorbed on bottom slab
selected_site_metal = config.get('settings', 'selected_site_metal')
reference_site_metal = functions.metal_fcc_111_high_symmetry_points(upper_slab,selected_site_metal)
upper_slab_coords_for_adatom_adsorption = functions.shift_slab_on_xy(upper_slab, reference_site_metal,max(cartesian_coord_bottom_slab, key=lambda x: x[2]))

# Reflect atomic coordinates of the upper slab with respect to its highest atom along z direction
plane_point, plane_normal = functions.calculate_plane_parameters(cartesian_coord_upper_slab)
reflected_coords =  functions.reflect_coord(upper_slab_coords_for_adatom_adsorption, plane_point, plane_normal)

# The distance between the two slabs is determined by the distance of the adatom when it is adsorbed on the upper slab  
shift_z = functions.calculate_shift_z_decorated_case(cartesian_coord_bottom_slab, cartesian_coord_adsorption, cartesian_coord_upper_slab)
shifted_coords =  functions.shift_slab_along_z(reflected_coords,shift_z)
x_relax = config.get('settings', 'x_relax')
y_relax = config.get('settings', 'y_relax')
z_relax = config.get('settings', 'z_relax')

functions.write_POSCAR_interface(upper_slab, bottom_slab_with_adatom, cartesian_coord_bottom_slab, x_relax, y_relax, z_relax, shifted_coords, a, b, c)
