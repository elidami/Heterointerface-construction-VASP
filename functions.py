import pandas as pd
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os

def extract_lattice_vectors(filename):
    '''This method reads the text file (VASP POSCAR or CONTCAR file)
       into a pandas DataFrame and extracts the lattice vectors.

       Args:
            filename: path to the text file 

        Returns:
            The lattice vectors a,b,c as a dictionary.'''
    df = pd.read_csv(filename, skiprows=lambda x: x not in [2, 3, 4], delim_whitespace=True, header=None, engine='python')
    
    a = np.array(df.iloc[0, 0:3]).astype(float)
    b = np.array(df.iloc[1, 0:3]).astype(float)
    c = np.array(df.iloc[2, 0:3]).astype(float)

    return {'a': a, 'b': b, 'c': c}

def extract_atomic_coordinates(filename):
    '''This method reads the text file (VASP CONTCAR file)
       into a pandas DataFrame and extracts the atomic coordinates.
       
       Args:
            filename: path to the text file 
       
       Returns:
            The atomic coordinates.'''
    df = pd.read_csv(filename, skiprows= 9, sep='\s+', header=None, engine='python')
    
    last_row_index = int(len(df.index)/2) #In order to exclude the series of numbers 0 (initial velocities of atoms), division by 2 is necessary.
    coordinates_df = df.iloc[:last_row_index, 0:3]
    coordinates = coordinates_df.values.astype(float)

    return coordinates
  
def direct_to_cartesian_coord(a,b,c, direct_coord):
    '''This method converts direct atomic coordinates to cartesian
       coordinates.
    
       Args:
            a,b,c: the three lattice vectors of the cell
            direct_coord: set of direct coordinates of the atoms
        
       Returns:
            Cartesian coordinates of the atoms as a list of numpy arrays.'''
    cartesian_coords = []
    for coord in direct_coord:
        cartesian_coord = coord[0] * a + coord[1] * b + coord[2] * c
        cartesian_coords.append(cartesian_coord)
    
    return cartesian_coords

def reflect_coord(coord_list, plane_point, plane_normal):
    '''This method reflects a list of points across a plane, which is
       defined by a point on the plane and its normal vector.
    
        Args:
            coord_list: list of numpy array of shape (3,) representing atomic coordinates.
            plane_point: a numpy array of shape (3,) representing a point on the plane with x, y, z coordinates.
            plane_normal: a numpy array of shape (3,) representing the normal vector of the plane.
        
        Returns:
            List of numpy arrays of shape (3,) representing the reflected points.'''
    reflected_coord = []
    for coord in coord_list:
        v = coord - plane_point
        d = np.dot(v, plane_normal) * 2
        reflected_point = coord - d * plane_normal
        reflected_coord.append(reflected_point)
    return reflected_coord

def shift_slab_along_z(coord_list, shift_z):
    shifted_coords = []
    for coord in coord_list:
        shifted_coord = np.array([coord[0], coord[1], coord[2] + shift_z])
        shifted_coords.append(shifted_coord)
    return shifted_coords

def C_111_high_symmetry_points(file_path, selected_site):
    original_file = file_path
    temporary_file = "CONTCAR"

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
        reference_site = high_symmetry_sites[0][1].coords
    elif selected_site == "hollow_hcp":
        reference_site = high_symmetry_sites[1][1].coords
    elif selected_site == "hollow_fcc":
        reference_site = high_symmetry_sites[3][1].coords
    else:
        print("Reference site not valid.")
        return None
    # Original file name restoration
    os.rename(temporary_file, original_file)
    return reference_site

def metal_fcc_111_high_symmetry_points(file_path, selected_site):
    original_file = file_path
    temporary_file = "CONTCAR"

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
        reference_site = high_symmetry_sites[0][1].coords
    elif selected_site == "hollow_hcp":
        reference_site = high_symmetry_sites[1][1].coords
    elif selected_site == "hollow_fcc":
        reference_site = high_symmetry_sites[2][1].coords
    else:
        print("Reference site not valid.")
        return None
    # Original file name restoration

    os.rename(temporary_file, original_file)
    return reference_site

def  shift_slab_on_xy(file_path, selected_site_Cu,selected_site_C):
    original_file = file_path
    temporary_file = "CONTCAR"

    # Renaming the original file as the temporary file
    os.rename(original_file, temporary_file)

    # Crystalline structure from POSCAR
    structure = Structure.from_file(temporary_file)
    shift_x = selected_site_C[0]-selected_site_Cu[0]
    shift_y = selected_site_C[1]-selected_site_Cu[1]
 # Shifting of the atomic coordinates in the POSCAR file
    for site in structure:
        site.coords[0] += shift_x
        site.coords[1] += shift_y
 # Original file name restoration
    os.rename(temporary_file, original_file)
 # Restituisci la lista di coordinate shiftate
    shifted_coords = [site.coords for site in structure]
    return shifted_coords

def distance_between_highest_z_values(coord):
    # Ordina le coordinate in base al valore di z in modo decrescente
    sorted_coord = sorted(coord, key=lambda x: x[2], reverse=True)

    # Trova i due valori di z più alti
    highest_z_values = sorted_coord[:2]
    # Calcola la distanza tra i due valori di z più alti
    distance = highest_z_values[0][2] - highest_z_values[1][2]

    return distance

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

def write_POSCAR_interface(input_file_upper, input_file_bottom,cartesian_coord_bottom_slab, x_relax, y_relax, z_relax, cartesian_coord_upper_slab, a,b,c):
    with open(input_file_upper, 'r') as f:
        upper_slab_lines = f.readlines()
        HeaderUpperSlab = upper_slab_lines[0].strip()
        AtomTypeUpperSlab = upper_slab_lines[5].strip()
        AtomNumberUpperSlab = upper_slab_lines[6].strip()

    with open(input_file_bottom, 'r') as f:
        bottom_slab_lines = f.readlines()
        HeaderBottomSlab = bottom_slab_lines[0].strip()
        AtomTypeBottomSlab = bottom_slab_lines[5].strip()
        AtomNumberBottomSlab = bottom_slab_lines[6].strip()

    AtomCoordsBottomSlab = write_coords(cartesian_coord_bottom_slab, x_relax, y_relax, z_relax)
    AtomCoordsUpperSlab = write_coords(cartesian_coord_upper_slab, x_relax, y_relax, z_relax)

    HeaderAndScalingFactor = ["INTERFACE {}/{}\n".format(HeaderBottomSlab, HeaderUpperSlab), "1.0\n"]
    LatticeVectors = ["{}\n".format(' '.join(["{:<20.16f}".format(value) for value in vector])) for vector in [a, b, c]]
    AtomTypes = ["{} {}\n".format(AtomTypeBottomSlab, AtomTypeUpperSlab)]
    AtomNumbers = ["{} {}\n".format(AtomNumberBottomSlab, AtomNumberUpperSlab)]

    with open('POSCAR', 'w') as f:
        f.writelines(HeaderAndScalingFactor)
        f.writelines(LatticeVectors)
        f.writelines(AtomTypes)
        f.writelines(AtomNumbers)
        f.writelines("Selective Dynamics\nCartesian\n")
        f.writelines(AtomCoordsBottomSlab)
        f.writelines(AtomCoordsUpperSlab)
