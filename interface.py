import pandas as pd
import numpy as np

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

def shift_coordinates(coord_list, shift_x, shift_y, shift_z):
    shifted_coords = []
    for coord in coord_list:
        shifted_coord = np.array([coord[0] + shift_x, coord[1] + shift_y, coord[2] + shift_z])
        shifted_coords.append(shifted_coord)
    return shifted_coords


# Extract lattice vectors
lattice_vectors = extract_lattice_vectors('bottom_slab.txt')
a = lattice_vectors['a']
b = lattice_vectors['b']
c = lattice_vectors['c']

# Extract atomic coordinates
atomic_coord_bottom_slab = extract_atomic_coordinates('bottom_slab.txt')
atomic_coord_upper_slab = extract_atomic_coordinates('upper_slab.txt')

# Convert direct coordinates to Cartesian coordinates
cartesian_coord_bottom_slab = direct_to_cartesian_coord(a, b, c, atomic_coord_bottom_slab)
cartesian_coord_upper_slab = direct_to_cartesian_coord(a, b, c, atomic_coord_upper_slab)

z = cartesian_coord_upper_slab[-1][2] if cartesian_coord_upper_slab else None
plane_point = np.array([0, 0,  z])
plane_normal = np.array([0, 0, 1])
reflected_coords = reflect_coord(cartesian_coord_upper_slab, plane_point, plane_normal)

shift_x = 0
shift_y = 0
interlayer_distance_bottom_slab =1.545855346
interlayer_distance_upper_slab =2.098396221
slabs_distance = (interlayer_distance_bottom_slab+interlayer_distance_upper_slab)/2

shift_z = cartesian_coord_bottom_slab[-1][2]+slabs_distance-cartesian_coord_upper_slab[-1][2]
shifted_coords = shift_coordinates(reflected_coords,shift_x,shift_y,shift_z)
