import functions
import numpy as np
import unittest


class TestDirectToCartesian(unittest.TestCase):
    '''This test class ensures that the "direct_to_cartesian_coord" function correctly converts
       three sets of direct coordinates into the corresponding cartesian coordinates. 

        Each test function is characterized by these steps:
            -Defined the three lattice vectors 'a', 'b', and 'c' representing the lattice of the structure
            -Define the direct coordinates 'direct_coord' as a numpy array containing two sets of direct coordinates.
            -Execute the "direct_to_cartesian_coord" function, passing the lattice vectors and direct coordinates as input.
            -Verify if the returned results match the expected Cartesian coordinates.
            
        Methods:
            -test_generic_direct_coords: test the case of positive, zero and negative coordinates.
            -test_lattice_vectors: test if changing the sign in the lattice vectors corresponds to changing signs of cartesian coordinates.'''

    def test_generic_direct_coords(self):
        a = np.array([1.00, 0.00, 0.00])
        b = np.array([0.00, 2.00, 0.00])
        c = np.array([0.00, 0.00, 3.00])
        direct_coord = np.array([[0.25, 0.50, 0.75], [0, 0, 0], [-0.10, -0.20, -0.30] ])
        result = functions.direct_to_cartesian_coord(a, b, c, direct_coord)
        expected_result = [
            np.array([0.25, 1, 2.25]),
            np.array([0, 0, 0]),
            np.array([-0.10, -0.40, -0.90])
        ]
        for i in range(len(result)):
            assert np.allclose(result[i], expected_result[i]), f"The obtained results for coordinate {i+1} do not match the expectations."

    def test_lattice_vectors(self):
        a = np.array([-1.00, 0.00, 0.00])
        b = np.array([0.00, -2.00, 0.00])
        c = np.array([0.00, 0.00, -3.00])
        direct_coord = np.array([[0.25, 0.50, 0.75], [0, 0, 0], [-0.10, -0.20, -0.30] ])
        result = functions.direct_to_cartesian_coord(a, b, c, direct_coord)
        expected_result = [
            np.array([-0.25, -1, -2.25]),
            np.array([0, 0, 0]),
            np.array([0.10, 0.40, 0.90])
        ]
        for i in range(len(result)):
            assert np.allclose(result[i], expected_result[i]), f"The obtained results for coordinate {i+1} do not match the expectations."