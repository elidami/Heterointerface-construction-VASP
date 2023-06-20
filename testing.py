import functions
import numpy as np
import os

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

    try:
        # Eseguire la funzione da testare
        result = functions.extract_atomic_coordinates(coordinates_file)

        # Verificare se i risultati sono conformi alle aspettative
        expected_result = np.array([[1.23, 4.56, 7.89],
                                   [2.34, 5.67, 8.90],
                                   [3.45, 6.78, 9.01]])
        assert np.array_equal(result, expected_result), "I risultati ottenuti non corrispondono alle aspettative"
    finally:
        # Cancellare il file di test
        os.remove(coordinates_file)

# Chiamare la funzione di test
test_extract_atomic_coordinates()


