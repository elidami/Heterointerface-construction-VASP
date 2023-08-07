import pandas as pd
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
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
            a,b,c: the three lattice vectors of the cell.
            direct_coord: set of direct coordinates of the atoms.
        
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
    '''This method shifts a list of points in the z direction.
    
        Args:
            coord_list: list of numpy array of shape (3,) representing atomic coordinates.
            shift_z: the amount of shift to be applied along the z-axis.
    
        Returns:
            List of numpy arrays of shape (3,) representing the shifted points along z.'''
    shifted_coords = []
    for coord in coord_list:
        shifted_coord = np.array([coord[0], coord[1], coord[2] + shift_z])
        shifted_coords.append(shifted_coord)
    return shifted_coords

def C_111_high_symmetry_points(file_path, selected_site):
    '''This method searches the high symmetry points of (1x1)C(111) slab with single 
        dangling bond termination.
    
        Args:
            file_path: The path to the input file (in POSCAR or CONTCAR format) 
                       containing the structure information.
            selected_site: a string indicating the selected high symmetry point. 
                           It can be "top", "hollow_hcp", or "hollow_fcc".

        Returns:
            A numpy array of shape (3,) representing the coordinates of the selected 
            high symmetry point. This will play the role of reference site for the coordinate 
            shift in subsequent calculations.'''
    original_file = file_path
    temporary_file = "CONTCAR" #Pymatgen's methods expect the files to be named according to a specific convention.
    os.rename(original_file, temporary_file)

    #Get crystalline structure from the file and information about symmetries
    structure = Structure.from_file(temporary_file)
    analyzer = SpacegroupAnalyzer(structure)
    symmetrized_structure = analyzer.get_symmetrized_structure()

    #Get high symmetry points of the structure
    high_symmetry_sites = symmetrized_structure.equivalent_sites

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
   
    os.rename(temporary_file, original_file)
    return reference_site

def metal_fcc_111_high_symmetry_points(file_path, selected_site):
    '''This method searches the high symmetry points of (111) surface of fcc metals.
    
        Args:
            file_path: The path to the input file (in POSCAR or CONTCAR format) 
                   containing the structure information.
            selected_site: a string indicating the selected high symmetry point. 
                       It can be "top", "hollow_hcp", or "hollow_fcc".
        Returns:
            A numpy array of shape (3,) representing the coordinates of the selected 
            high symmetry point. This will play the role of reference site for the coordinate 
            shift in subsequent calculations.'''
    original_file = file_path
    temporary_file = "CONTCAR" #Pymatgen's methods expect the files to be named according to a specific convention.
    os.rename(original_file, temporary_file)

    #Get crystalline structure from the file and information about symmetries
    structure = Structure.from_file(temporary_file)
    analyzer = SpacegroupAnalyzer(structure)
    symmetrized_structure = analyzer.get_symmetrized_structure()  

    #Get high symmetry points of the structure
    high_symmetry_sites = symmetrized_structure.equivalent_sites

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

    os.rename(temporary_file, original_file)
    return reference_site

def  shift_slab_on_xy(file_path, selected_site_metal, selected_site_C):
    '''This method shfits the atomic coordinates in the x and y directions.
    
        Args:
            file_path: The path to the input file (in POSCAR or CONTCAR format) 
                   containing the structure information.
            selected_site_metal: A numpy array of shape (3,) representing the coordinates of a 
                              selected reference site of the metal(111) slab.
            selected_site_C: A numpy array of shape (3,) representing the coordinates of a 
                             selected reference site of (1x1)C(111) slab.
                  
        Returns:
            A list of numpy arrays with shape (3,) representing the shifted atomic coordinates 
            along the x and y directions.
            The atomic coordinates are modified based on the difference between the selected reference
            site for metal and C slabs.
            The function restores the original file name after the shift operation. '''
    original_file = file_path
    temporary_file = "CONTCAR"
    os.rename(original_file, temporary_file)

    #Crystalline structure from file
    structure = Structure.from_file(temporary_file)
    shift_x = selected_site_C[0]-selected_site_metal[0]
    shift_y = selected_site_C[1]-selected_site_metal[1]
    #Shifting of the atomic coordinates in the file
    for site in structure:
        site.coords[0] += shift_x
        site.coords[1] += shift_y

    os.rename(temporary_file, original_file)

    shifted_coords = [site.coords for site in structure]
    return shifted_coords

def distance_between_highest_z_values(coord):
    '''This method calculates the distance between the two highest z-values 
       in a list of coordinates.

       Args:
            coord: A list of numpy arrays with shape (3,) representing the atomic coordinates.

       Returns:
            The distance between the two highest z-values among the coordinates. 
            The function first sorts the coordinates based on their z-values in descending order. 
            It then selects the two coordinates with the highest z-values and calculates the 
            difference between these two values to obtain the distance.'''
    sorted_coord = sorted(coord, key=lambda x: x[2], reverse=True)
    highest_z_values = sorted_coord[:2]
    distance = highest_z_values[0][2] - highest_z_values[1][2]

    return distance

def write_coords(coords, x_relax, y_relax, z_relax):
    '''This method converts a list of atomic coordinates and relaxation options 
       into a formatted string representation.

       Args:
            coords: A list of numpy arrays with shape (3,) representing the atomic coordinates.
            x_relax: A boolean value indicating whether the relaxation is allowed along the x-direction.
            y_relax: A boolean value indicating whether the relaxation is allowed along the y-direction.
            z_relax: A boolean value indicating whether the relaxation is allowed along the z-direction.

       Returns:
            A list of strings, each representing a formatted line of atomic coordinates along with 
            the relaxation options. The function iterates over the input coordinates and converts each 
            coordinate to a string with a specific format. The relaxation options (x_relax, y_relax, z_relax) 
            are represented by "T" (True) or "F" (False) and appended to the string.'''
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

def write_POSCAR_interface(input_file_upper, input_file_bottom,cartesian_coord_bottom_slab, 
                           x_relax, y_relax, z_relax, cartesian_coord_upper_slab, a,b,c):
    '''This method creates a new POSCAR file representing the interface between two slabs.

        Args:
            input_file_upper: the path to the input file containing the structure information of the upper slab.
            input_file_bottom: the path to the input file containing the structure information of the bottom slab.
            cartesian_coord_bottom_slab: a list of numpy arrays representing the Cartesian coordinates of the bottom slab atoms.
            x_relax, y_relax, z_relax: boolean values indicating whether relaxation is allowed along x, y, and z directions, respectively.
            cartesian_coord_upper_slab: a list of numpy arrays representing the Cartesian coordinates of the upper slab atoms.
            a, b, c: three lattice vectors representing the lattice of the interface structure.

        Returns:
            A new POSCAR file that represents the interface between the two slabs. 
            It combines the header, lattice vectors, atom types, and atom numbers from the input files for the bottom and upper slabs. 
            The Cartesian coordinates of the atoms in the interface supercell are written along with their relaxation options.'''
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
