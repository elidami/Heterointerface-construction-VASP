import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os
import functions

# Extract lattice vectors
lattice_vectors = functions.extract_lattice_vectors('bottom_slab.txt')
a = lattice_vectors['a']
b = lattice_vectors['b']
c = lattice_vectors['c']

# Extract atomic coordinates
atomic_coord_bottom_slab = functions.extract_atomic_coordinates('bottom_slab.txt')
atomic_coord_upper_slab =  functions.extract_atomic_coordinates('upper_slab.txt')

# Convert direct coordinates to Cartesian coordinates
cartesian_coord_bottom_slab =  functions.direct_to_cartesian_coord(a, b, c, atomic_coord_bottom_slab)
cartesian_coord_upper_slab =  functions.direct_to_cartesian_coord(a, b, c, atomic_coord_upper_slab)

z = cartesian_coord_upper_slab[-1][2] if cartesian_coord_upper_slab else None
plane_point = np.array([0, 0,  z])
plane_normal = np.array([0, 0, 1])
reflected_coords =  functions.reflect_coord(cartesian_coord_upper_slab, plane_point, plane_normal)

interlayer_distance_bottom_slab =1.545855346
interlayer_distance_upper_slab =2.098396221
slabs_distance = (interlayer_distance_bottom_slab+interlayer_distance_upper_slab)/2

shift_z = cartesian_coord_bottom_slab[-1][2]+slabs_distance-cartesian_coord_upper_slab[-1][2]
shifted_coords =  functions.shift_coordinates(reflected_coords,shift_z)

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

def C_high_symmetry_points(file_path, selected_site):
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
        reference_site = high_symmetry_sites[0][1]
    elif selected_site == "hollow_hcp":
        reference_site = high_symmetry_sites[1][1]
    elif selected_site == "hollow_fcc":
        reference_site = high_symmetry_sites[3][1]
    else:
        print("Reference site not valid.")
        return None
    # Original file name restoration
    os.rename(temporary_file, original_file)
    return reference_site

def  shift_slab_on_xy(file_path, selected_site_Cu,selected_site_C):
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
    selected_site_Cu = selected_site_Cu 

    if selected_site_Cu == "top":
        reference_site = high_symmetry_sites[0][1]
    elif selected_site_Cu == "hollow_hcp":
        reference_site = high_symmetry_sites[1][1]
    elif selected_site_Cu == "hollow_fcc":
        reference_site = high_symmetry_sites[2][1]
    else:
        print("Reference site not valid.")
        return None

    shift_x = selected_site_C.coords[0]-reference_site.coords[0]
    shift_y = selected_site_C.coords[1]-reference_site.coords[1]


    # Shifting of the atomic coordinates in the POSCAR file
    for site in structure:
        site.coords[0] += shift_x
        site.coords[1] += shift_y


    # Original file name restoration
    os.rename(temporary_file, original_file)
    # Restituisci la lista di coordinate shiftate
    shifted_coords = [site.coords for site in structure]

    return shifted_coords


selected_site_Cu = "top"
selected_site_C = "top"

reference_site_C = C_high_symmetry_points("bottom_slab.txt", selected_site_C)
shifted_upper_slab_on_xy = shift_slab_on_xy("reflected_upper_slab.txt", selected_site_Cu,reference_site_C)


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

AtomCoordsBottomSlab = write_coords(cartesian_coord_bottom_slab, x_relax, y_relax, z_relax)
AtomCoordsUpperSlab = write_coords(shifted_upper_slab_on_xy, x_relax, y_relax, z_relax)       

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
