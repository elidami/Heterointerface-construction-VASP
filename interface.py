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
