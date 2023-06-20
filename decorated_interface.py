import numpy as np
from pymatgen.core.structure import Structure
import os
import functions

def distance_between_highest_z_values(coord):
    # Ordina le coordinate in base al valore di z in modo decrescente
    sorted_coord = sorted(coord, key=lambda x: x[2], reverse=True)

    # Trova i due valori di z più alti
    highest_z_values = sorted_coord[:2]
    # Calcola la distanza tra i due valori di z più alti
    distance = highest_z_values[0][2] - highest_z_values[1][2]

    return distance



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


def shift_for_adatom_adsorption_on_upper_slab(file_path, selected_site):
    original_file = file_path
    temporary_file = "CONTCAR"

    # Renaming the original file as the temporary file
    os.rename(original_file, temporary_file)

    # Crystalline structure from POSCAR
    structure = Structure.from_file(temporary_file)

    shift_x = cartesian_coord_bottom_slab[0][0]-selected_site.coords[0]
    shift_y = cartesian_coord_bottom_slab[0][1]-selected_site.coords[1]

    # Shifting of the atomic coordinates in the POSCAR file
    for site in structure:
        site.coords[0] += shift_x
        site.coords[1] += shift_y


    # Original file name restoration
    os.rename(temporary_file, original_file)
    # Restituisci la lista di coordinate shiftate
    shifted_coords = [site.coords for site in structure]

    return shifted_coords

file_path = "decorated_interface_files/upper_slab.txt"
selected_site_Cu = "hollow_fcc"

reference_site_Cu = functions.metal_fcc_111_high_symmetry_points("decorated_interface_files/upper_slab.txt",selected_site_Cu)
upper_slab_coords_for_adatom_adsorption = shift_for_adatom_adsorption_on_upper_slab(file_path, reference_site_Cu)


z = upper_slab_coords_for_adatom_adsorption[-1][2] if upper_slab_coords_for_adatom_adsorption else None
plane_point = np.array([0, 0,  z])
plane_normal = np.array([0, 0, 1])
reflected_coords =  functions.reflect_coord(upper_slab_coords_for_adatom_adsorption, plane_point, plane_normal)
# Convert direct coordinates to Cartesian coordinates
cartesian_coord_adsorption =  functions.direct_to_cartesian_coord(a_adsorption, b_adosrption, c_adsorption, atomic_coord_adsorption)
shift_z = cartesian_coord_bottom_slab[0][2]+distance_between_highest_z_values(cartesian_coord_adsorption)-cartesian_coord_upper_slab[-1][2]
shifted_coords =  functions.shift_slab_along_z(reflected_coords,shift_z)

x_relax = True
y_relax = True
z_relax = True


functions.write_POSCAR_interface('decorated_interface_files/upper_slab.txt', 'decorated_interface_files/bottom_slab_with_adatom.txt', cartesian_coord_bottom_slab, x_relax, y_relax, z_relax, shifted_coords, a, b, c)













































