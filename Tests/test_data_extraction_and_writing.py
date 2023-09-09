import functions
import numpy as np
import os
import unittest

class TestLatticeVectorsExtraction(unittest.TestCase):
    def test_extract_lattice_vectors(self):
        '''This test ensures that the "extract_lattice_vectors" function correctly reads the three generic lattice 
           vectors from the input file and returns them in the expected format, represented as a dictionary 
           with keys 'a', 'b', and 'c' mapping to their respective numpy array values.

           Test Steps:
            -Create a temporary example file named "lattice_vectors.csv" with lattice vectors.
            -Execute the "extract_lattice_vectors" function, passing the file path of the created file as input.
            -Verify if the returned results match the expected lattice vectors.'''
        lattice_file = 'lattice_vectors.csv'
        with open(lattice_file, 'w') as f:
            f.write('# Comments or header information\n')
            f.write('# More comments or header information\n')
            f.write('3.45 0.00 0.00\n')
            f.write('0.00 3.45 0.00\n')
            f.write('0.00 0.00 5.00\n')

        result = functions.extract_lattice_vectors(lattice_file)

        expected_result = {
            'a': np.array([3.45, 0.00, 0.00]),
            'b': np.array([0.00, 3.45, 0.00]),
            'c': np.array([0.00, 0.00, 5.00])
        }
        assert np.array_equal(result['a'], expected_result['a']), "The obtained results for 'a' do not match the expectations."
        assert np.array_equal(result['b'], expected_result['b']), "The obtained results for 'b' do not match the expectations."
        assert np.array_equal(result['c'], expected_result['c']), "The obtained results for 'c' do not match the expectations."
    
        os.remove(lattice_file)

class TestAtomicCoordinatesExtraction(unittest.TestCase):
    def test_extract_atomic_coordinates(self):
        '''This test ensures that the "extract_atomic_coordinates" function extracts correctly 
           the three generic atomic coordinates written in the example text file.

           Test Steps:
                -Create a temporary example file named "atomic_coordinates.csv" with atomic coordinates.
                -Execute the "extract_atomic_coordinates" function, passing the file path of the created file as input.
                -Verify if the returned results match the expected atomic coordinates.'''
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
    
        result = functions.extract_atomic_coordinates(coordinates_file)

        expected_result = np.array([[1.23, 4.56, 7.89],
                                   [2.34, 5.67, 8.90],
                                   [3.45, 6.78, 9.01]])
        assert np.array_equal(result, expected_result), "The obtained results do not match the expectations."

        os.remove(coordinates_file)

class TestWriteCoords(unittest.TestCase):
    coords = [
    np.array([2.5035679999999956 ,  1.4454350000000022,   9.6615339926441077 ]),
    np.array([1.2519332746767966 ,  0.7228036064245595,   0.0583534541551990]),
    np.array([5.0071359999999903,   1.4458782853606842  , 0.4115548962189028])]

    def test_write_coords(self):
        '''This function is a unit test for the "write_coords" function.

        Test Steps:
            - Define sample input data, including a list of NumPy arrays representing atomic coordinates (coords)
              and boolean variables 'x_relax', 'y_relax', and 'z_relax'.
            - Call the "write_coords" function with the provided input data to get the result.
            - Compare the obtained result with the expected output to validate the accuracy of the function.'''
        x_relax = True
        y_relax = False
        z_relax = True

        result = functions.write_coords(self.coords, x_relax, y_relax, z_relax)

        expected_output = ['2.5035679999999956   1.4454350000000022   9.6615339926441077    T  F  T\n', 
                           '1.2519332746767966   0.7228036064245595   0.0583534541551990    T  F  T\n', 
                           '5.0071359999999903   1.4458782853606842   0.4115548962189028    T  F  T\n']

        assert result == expected_output, "Test failed: Output doesn't match the expected result."

    def test_all_false(self):
        '''This function is a unit test for the "write_coords" function.
        Test Steps:
            - Define sample input data, including a list of NumPy arrays representing atomic coordinates (co
              and boolean variables 'x_relax', 'y_relax', and 'z_relax'.
            - Call the "write_coords" function with the provided input data to get the result.
            - Compare the obtained result with the expected output to validate the accuracy of the function.'''
        x_relax = False
        y_relax = False
        z_relax = False
        result = functions.write_coords(self.coords, x_relax, y_relax, z_relax)
        expected_output = ['2.5035679999999956   1.4454350000000022   9.6615339926441077    F  F  F\n', 
                           '1.2519332746767966   0.7228036064245595   0.0583534541551990    F  F  F\n', 
                           '5.0071359999999903   1.4458782853606842   0.4115548962189028    F  F  F\n']
        assert result == expected_output, "Test failed: Output doesn't match the expected result."

    def test_all_true(self):
        '''This function is a unit test for the "write_coords" function.
        Test Steps:
            - Define sample input data, including a list of NumPy arrays representing atomic coordinates (co
              and boolean variables 'x_relax', 'y_relax', and 'z_relax'.
            - Call the "write_coords" function with the provided input data to get the result.
            - Compare the obtained result with the expected output to validate the accuracy of the function.'''
        x_relax = True
        y_relax = True
        z_relax = True
        result = functions.write_coords(self.coords, x_relax, y_relax, z_relax)
        expected_output = ['2.5035679999999956   1.4454350000000022   9.6615339926441077    T  T  T\n', 
                           '1.2519332746767966   0.7228036064245595   0.0583534541551990    T  T  T\n', 
                           '5.0071359999999903   1.4458782853606842   0.4115548962189028    T  T  T\n']
        assert result == expected_output, "Test failed: Output doesn't match the expected result."


    def test_no_coordinates(self):
        coords = []  
        x_relax = True
        y_relax = True
        z_relax = True

        result = functions.write_coords(coords, x_relax, y_relax, z_relax)

        expected_output = []  

        assert result == expected_output, "Test failed: Output doesn't match the expected result."