import functions
import numpy as np
import os
from pymatgen.core.structure import Structure, Lattice
from pymatgen.core.surface import SlabGenerator


def test_extract_lattice_vectors():
    # Creare un file di esempio con i vettori del reticolo
    lattice_file = 'lattice_vectors.csv'
    with open(lattice_file, 'w') as f:
        f.write('# Comments or header information\n')
        f.write('# More comments or header information\n')
        f.write('3.45 0.00 0.00\n')
        f.write('0.00 3.45 0.00\n')
        f.write('0.00 0.00 5.00\n')

    # Eseguire la funzione da testare
    result = functions.extract_lattice_vectors(lattice_file)

    # Verificare se i risultati sono conformi alle aspettative
    expected_result = {
        'a': np.array([3.45, 0.00, 0.00]),
        'b': np.array([0.00, 3.45, 0.00]),
        'c': np.array([0.00, 0.00, 5.00])
    }
    assert np.array_equal(result['a'], expected_result['a']), "I risultati ottenuti per 'a' non corrispondono alle aspettative"
    assert np.array_equal(result['b'], expected_result['b']), "I risultati ottenuti per 'b' non corrispondono alle aspettative"
    assert np.array_equal(result['c'], expected_result['c']), "I risultati ottenuti per 'c' non corrispondono alle aspettative"
 
    os.remove(lattice_file)

# Chiamare la funzione di test
test_extract_lattice_vectors()

def test_extract_atomic_coordinates():
    # Creare un file di esempio con le coordinate atomiche
    coordinates_file = 'atomic_coordinates.csv'
    with open(coordinates_file, 'w') as f:
        other_info = '#Other information\n'
        f.writelines([other_info] * 9)
        f.write('1.23 4.56 7.89\n')
        f.write('2.34 5.67 8.90\n')
        f.write('3.45 6.78 9.01\n')
        f.write('0 0 0\n')
        f.write('0 0 0\n')
        f.write('0 0 0\n')
    # Eseguire la funzione da testare
    result = functions.extract_atomic_coordinates(coordinates_file)

    # Verificare se i risultati sono conformi alle aspettative
    expected_result = np.array([[1.23, 4.56, 7.89],
                               [2.34, 5.67, 8.90],
                               [3.45, 6.78, 9.01]])
    assert np.array_equal(result, expected_result), "I risultati ottenuti non corrispondono alle aspettative"
    
    # Cancellare il file di test
    os.remove(coordinates_file)

# Chiamare la funzione di test
test_extract_atomic_coordinates()

def test_direct_to_cartesian_coord():
    # Definire i vettori del reticolo
    a = np.array([1.00, 0.00, 0.00])
    b = np.array([0.00, 2.00, 0.00])
    c = np.array([0.00, 0.00, 3.00])

    # Definire le coordinate dirette
    direct_coord = np.array([[0.25, 0.50, 0.75], [0.10, 0.20, 0.30]])

    # Eseguire la funzione da testare
    result = functions.direct_to_cartesian_coord(a, b, c, direct_coord)

    # Verificare se i risultati sono conformi alle aspettative
    expected_result = [
        np.array([0.25, 1, 2.25]),
        np.array([0.10, 0.40, 0.90])
    ]
    for i in range(len(result)):
        assert np.allclose(result[i], expected_result[i]), f"I risultati ottenuti per la coordinata {i+1} non corrispondono alle aspettative"


# Chiamare la funzione di test
test_direct_to_cartesian_coord()

def test_reflect_coord():
    # Definire i dati di input per il test
    coord_list = [
        np.array([1.0, 2.0, 3.0]),
        np.array([4.0, 5.0, 6.0]),
        np.array([7.0, 8.0, 9.0])
    ]
    plane_point = np.array([0.0, 0.0, 0.0])
    plane_normal = np.array([1.0, 0.0, 0.0])

    # Eseguire la funzione da testare
    result = functions.reflect_coord(coord_list, plane_point, plane_normal)

    # Definire i risultati attesi
    expected_result = [
        np.array([-1.0, 2.0, 3.0]),
        np.array([-4.0, 5.0, 6.0]),
        np.array([-7.0, 8.0, 9.0])
    ]

    # Verificare se i risultati sono conformi alle aspettative
    for i in range(len(result)):
        assert np.array_equal(result[i], expected_result[i]), f"I risultati ottenuti per il punto {i+1} non corrispondono alle aspettative"

# Chiamare la funzione di test
test_reflect_coord()

def test_shift_slab_along_z():
    # Definire i dati di input per il test
    coord_list = [
        np.array([1.0, 2.0, 3.0]),
        np.array([4.0, 5.0, 6.0]),
        np.array([7.0, 8.0, 9.0])
    ]
    shift_z = 2.0

    # Eseguire la funzione da testare
    result = functions.shift_slab_along_z(coord_list, shift_z)

    # Definire i risultati attesi
    expected_result = [
        np.array([1.0, 2.0, 5.0]),
        np.array([4.0, 5.0, 8.0]),
        np.array([7.0, 8.0, 11.0])
    ]

    for i in range(len(result)):
        assert np.array_equal(result[i], expected_result[i]), f"I risultati ottenuti per il punto {i+1} non corrispondono alle aspettative"

# Chiamare la funzione di test
test_shift_slab_along_z()

def test_C_111_high_symmetry_points():
    # Creazione della superficie 111 del diamante
    lattice_constant = 3.57  # Costante di reticolo del diamante
    lattice = Lattice.cubic(lattice_constant)
    structure = Structure.from_spacegroup(227, lattice, ["C"], [[0, 0, 0]])
    print(structure)
    #miller_index = (1, 1, 1)
    #min_slab_size = 10
    #min_vacuum_size = 15
    #slabgen = SlabGenerator(structure, miller_index, min_slab_size, min_vacuum_size)
    #slab = slabgen.get_slab()

test_C_111_high_symmetry_points()










































