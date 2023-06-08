import pandas as pd
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os

def extract_lattice_vectors(filename):
    # Read the text file into a pandas DataFrame, skipping the header row
    df = pd.read_csv(filename, skiprows=lambda x: x not in [2, 3, 4], delim_whitespace=True, header=None, engine='python')
    
    a = np.array(df.iloc[0, 0:3]).astype(float)
    b = np.array(df.iloc[1, 0:3]).astype(float)
    c = np.array(df.iloc[2, 0:3]).astype(float)

    # Return the extracted values as a dictionary
    return {'a': a, 'b': b, 'c': c}

def extract_atomic_coordinates(filename):
    df = pd.read_csv(filename, skiprows= 9, sep='\s+', header=None, engine='python')
    
    last_row_index = int(len(df.index)/2)
    coordinates_df = df.iloc[:last_row_index, 0:3]
    coordinates = coordinates_df.values.astype(float)

    return coordinates
  
def direct_to_cartesian_coord(a,b,c, direct_coord):
    cartesian_coords = []
    for coord in direct_coord:
        cartesian_coord = coord[0] * a + coord[1] * b + coord[2] * c
        cartesian_coords.append(cartesian_coord)
    
    return cartesian_coords

def reflect_coord(coord_list, plane_point, plane_normal):
    """
    Reflects a point across a plane defined by a point on the plane and its normal vector.
    
    Args:
        point: A numpy array of shape (3,) representing a point with x, y, z coordinates.
        plane_point: A numpy array of shape (3,) representing a point on the plane with x, y, z coordinates.
        plane_normal: A numpy array of shape (3,) representing the normal vector of the plane.
        
    Returns:
        A numpy array of shape (3,) representing the reflected point.
    """
    reflected_coord = []
    for coord in coord_list:
        v = coord - plane_point
        d = np.dot(v, plane_normal) * 2
        reflected_point = coord - d * plane_normal
        reflected_coord.append(reflected_point)
    return reflected_coord

def shift_coordinates(coord_list, shift_z):
    shifted_coords = []
    for coord in coord_list:
        shifted_coord = np.array([coord[0], coord[1], coord[2] + shift_z])
        shifted_coords.append(shifted_coord)
    return shifted_coords

def highest_z_value(coord):
    coord_excluded_first = coord[1:]
    # Ordina le coordinate in base al valore di z in modo crescente
    sorted_coord = sorted(coord_excluded_first, key=lambda x: x[2])

    # Trova il valore di z più alto
    highest_z_value = sorted_coord[-1][2]

    return highest_z_value





# Extract lattice vectors
lattice_vectors = extract_lattice_vectors('bottom_slab_with_adatom.txt')
a = lattice_vectors['a']
b = lattice_vectors['b']
c = lattice_vectors['c']

# Extract atomic coordinates
atomic_coord_bottom_slab = extract_atomic_coordinates('bottom_slab_with_adatom.txt')
atomic_coord_upper_slab = extract_atomic_coordinates('upper_slab.txt')

# Convert direct coordinates to Cartesian coordinates
cartesian_coord_bottom_slab = direct_to_cartesian_coord(a, b, c, atomic_coord_bottom_slab)
cartesian_coord_upper_slab = direct_to_cartesian_coord(a, b, c, atomic_coord_upper_slab)

z = cartesian_coord_upper_slab[-1][2] if cartesian_coord_upper_slab else None
plane_point = np.array([0, 0,  z])
plane_normal = np.array([0, 0, 1])
reflected_coords = reflect_coord(cartesian_coord_upper_slab, plane_point, plane_normal)


lattice_vectors_adsorption = extract_lattice_vectors('adsorption_on_upper_slab.txt')
a_adsorption = lattice_vectors_adsorption['a']
b_adosrption = lattice_vectors_adsorption['b']
c_adsorption = lattice_vectors_adsorption['c']

# Extract atomic coordinates
atomic_coord_adsorption = extract_atomic_coordinates('adsorption_on_upper_slab.txt')

# Convert direct coordinates to Cartesian coordinates
cartesian_coord_adsorption = direct_to_cartesian_coord(a_adsorption, b_adosrption, c_adsorption, atomic_coord_adsorption)
adatom_slab_distance = cartesian_coord_adsorption[0][2]-highest_z_value(cartesian_coord_adsorption)
shift_z = cartesian_coord_bottom_slab[0][2]+adatom_slab_distance-cartesian_coord_upper_slab[-1][2]
shifted_coords = shift_coordinates(reflected_coords,shift_z)



OutputFile = open("reflected_upper_slab.txt","w")
with open('upper_slab.txt', 'r') as f:
    upper_slab_lines = f.readlines()
    HeaderUpperSlab = upper_slab_lines[0].strip()
    AtomTypeUpperSlab = upper_slab_lines[5].strip()
    AtomNumberUpperSlab = upper_slab_lines[6].strip()

AtomCoordsUpperSlab = ["{}\n".format(' '.join(["{:<20.16f}".format(coord) for coord in coordinates.astype(float)])
                        ) for coordinates in shifted_coords]
HeaderAndScalingFactor = ["{}\n".format(HeaderUpperSlab),"1.0\n"]
LatticeVectors = ["{}\n".format(' '.join(["{:<20.16f}".format(value) for value in vector])) for vector in [a, b, c]]
AtomTypes = ["{}\n".format(AtomTypeUpperSlab)]
AtomNumbers = ["{}\n".format(AtomNumberUpperSlab)]

OutputFile.writelines(HeaderAndScalingFactor)
OutputFile.writelines(LatticeVectors)
OutputFile.writelines(AtomTypes)
OutputFile.writelines(AtomNumbers)
OutputFile.writelines("Cartesian\n")
OutputFile.writelines(AtomCoordsUpperSlab)
OutputFile.close()

def shift_for_adatom_adsorption_on_upper_slab(file_path, selected_site):
    original_file = file_path
    temporary_file = "POSCAR"

    # Renaming the original file as the temporary file
    os.rename(original_file, temporary_file)

    # Crystalline structure from POSCAR
    structure = Structure.from_file(temporary_file)

    # Get information about structure symmetry
    analyzer = SpacegroupAnalyzer(structure)
    symmetrized_structure = analyzer.get_symmetrized_structure()

    # Get high symmetry points of the structure
    high_symmetry_sites = symmetrized_structure.equivalent_sites

    # Select the reference site for the coordinates shift
    selected_site = selected_site 

    if selected_site == "top":
        reference_site = high_symmetry_sites[0][1]
    elif selected_site == "hollow_hcp":
        reference_site = high_symmetry_sites[1][1]
    elif selected_site == "hollow_fcc":
        reference_site = high_symmetry_sites[2][1]
    else:
        print("Reference site not valid.")
        return

    shift_x = cartesian_coord_bottom_slab[0][0]-reference_site.coords[0]
    shift_y = cartesian_coord_bottom_slab[0][1]-reference_site.coords[1]

    # Shifting of the atomic coordinates in the POSCAR file
    for site in structure:
        site.coords[0] += shift_x
        site.coords[1] += shift_y


    # Original file name restoration
    os.rename(temporary_file, original_file)
    # Restituisci la lista di coordinate shiftate
    shifted_coords = [site.coords for site in structure]

    return shifted_coords

file_path = "reflected_upper_slab.txt"
selected_site = "hollow_fcc"

upper_slab_coords_for_adatom_adsorption = shift_for_adatom_adsorption_on_upper_slab(file_path, selected_site)
print(upper_slab_coords_for_adatom_adsorption)

OutputFile = open("POSCAR","w")
with open('upper_slab.txt', 'r') as f:
    upper_slab_lines = f.readlines()
    HeaderUpperSlab = upper_slab_lines[0].strip()
    AtomTypeUpperSlab = upper_slab_lines[5].strip()
    AtomNumberUpperSlab = upper_slab_lines[6].strip()


with open('bottom_slab_with_adatom.txt', 'r') as f:
    bottom_slab_lines = f.readlines()
    HeaderBottomSlab = bottom_slab_lines[0].strip()
    AtomTypeBottomSlab = bottom_slab_lines[5].strip()
    AtomNumberBottomSlab = bottom_slab_lines[6].strip()



def write_coords(coords, x_relax, y_relax, z_relax):
    atom_coords = []
    for coord in coords:
        coord_string = "{}  {}  {}  {}\n".format(
                    ' '.join(["{:<20.16f}".format(c) for c in coord.astype(float)]),
                    "T" if x_relax else "F",
                    "T" if y_relax else "F",
                    "T" if z_relax else "F"
        )
        atom_coords.append(coord_string)
    return atom_coords  

x_relax = True
y_relax = True
z_relax = True

AtomCoordsBottomSlab = write_coords(cartesian_coord_bottom_slab, x_relax, y_relax, z_relax)
AtomCoordsUpperSlab = write_coords(upper_slab_coords_for_adatom_adsorption, x_relax, y_relax, z_relax)       

HeaderAndScalingFactor = ["INTERFACE {}/{}\n".format(HeaderBottomSlab,HeaderUpperSlab),"1.0\n"]
LatticeVectors = ["{}\n".format(' '.join(["{:<20.16f}".format(value) for value in vector])) for vector in [a, b, c]]
AtomTypes = ["{} {}\n".format(AtomTypeBottomSlab,AtomTypeUpperSlab)]
AtomNumbers = ["{} {}\n".format(AtomNumberBottomSlab,AtomNumberUpperSlab)]

OutputFile.writelines(HeaderAndScalingFactor)
OutputFile.writelines(LatticeVectors)
OutputFile.writelines(AtomTypes)
OutputFile.writelines(AtomNumbers)
OutputFile.writelines("Selective Dynamics\nCartesian\n")
OutputFile.writelines(AtomCoordsBottomSlab)
OutputFile.writelines(AtomCoordsUpperSlab)
OutputFile.close()













































