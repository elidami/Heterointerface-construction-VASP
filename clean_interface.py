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


OutputFile = open("POSCAR","w")
with open('upper_slab.txt', 'r') as f:
    upper_slab_lines = f.readlines()
    HeaderUpperSlab = upper_slab_lines[0].strip()
    AtomTypeUpperSlab = upper_slab_lines[5].strip()
    AtomNumberUpperSlab = upper_slab_lines[6].strip()


with open('bottom_slab.txt', 'r') as f:
    bottom_slab_lines = f.readlines()
    HeaderBottomSlab = bottom_slab_lines[0].strip()
    AtomTypeBottomSlab = bottom_slab_lines[5].strip()
    AtomNumberBottomSlab = bottom_slab_lines[6].strip()

x_relax = True
y_relax = True
z_relax = True
fixed_atoms = 0

def lowest_z_rows(coord,fixed_atoms):
    sorted_coord = sorted(coord, key=lambda x: x[2])
    fixed_atoms = fixed_atoms
    # Search lowest value of z
    lowest_z_values = sorted(set([row[2] for row in sorted_coord]))[:fixed_atoms]

    # Takes all the lines that have the lowest z values
    lowest_z_rows = [row for row in sorted_coord if row[2] in lowest_z_values]

    return lowest_z_rows

lowest_z_lines = lowest_z_rows(cartesian_coord_bottom_slab,fixed_atoms)

print(lowest_z_lines)



AtomCoordsBottomSlab = []
for coordinates in cartesian_coord_bottom_slab:
    if fixed_atoms != 0:
        if any(np.array_equal(coordinates, lowest) for lowest in lowest_z_lines):
            coord_string = "{}  {}\n".format(
                '\t'.join(["{:<20.16f}".format(coord) for coord in coordinates.astype(float)]),
                "F\tF\tF"
            )
        else:
            coord_string = "{}  {}  {}  {}\n".format(
                '\t'.join(["{:<20.16f}".format(coord) for coord in coordinates.astype(float)]),
                "T" if x_relax else "F",
                "T" if y_relax else "F",
                "T" if z_relax else "F"
            )
    else:
        coord_string = "{}  {}  {}  {}\n".format(
            '\t'.join(["{:<20.16f}".format(coord) for coord in coordinates.astype(float)]),
            "T" if x_relax else "F",
            "T" if y_relax else "F",
            "T" if z_relax else "F"
)
    AtomCoordsBottomSlab.append(coord_string)

AtomCoordsUpperSlab = ["{}  {}  {}  {}\n".format('\t'.join(["{:<20.16f}".format(coord) for coord in coordinates.astype(float)]),
                        "T" if x_relax else "F",
                        "T" if y_relax else "F",
                        "T" if z_relax else "T"
                         ) for coordinates in shifted_coords]

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
