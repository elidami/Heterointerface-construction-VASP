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

