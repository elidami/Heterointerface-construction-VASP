import functions
import numpy as np
import unittest

class TestShiftAlongZ(unittest.TestCase):
    '''This class represents a unit test for the "shift_slab_along_z" function.

            Attributes:
                 -coord_list: list of three numpy arrays of shape (3,) representing atomic positions.

            Methods:
                 -test_generic_shift
                 -test_shift_zero
                 -test_return_to_initial_conditions'''
    coord_list = [
    np.array([1.0, 2.0, 3.0]),
    np.array([4.0, 5.0, 6.0]),
    np.array([7.0, 8.0, 9.0])]

    def test_generic_shift(self):
        '''This function test a generic positive shift along the z direction.

            Test Steps:
                -Define a value 'shift_z', representing the amount of displacement in the z-direction.
                -Execute the "shift_slab_along_z" function, passing the coordinate list and shift value as input.
                -Verify if the returned results match the expected coordinates after the z-direction shift.'''
        shift_z = 2.0
        result = functions.shift_slab_along_z(self.coord_list, shift_z)
        expected_result = [
            np.array([1.0, 2.0, 5.0]),
            np.array([4.0, 5.0, 8.0]),
            np.array([7.0, 8.0, 11.0])]

        for i in range(len(result)):
            assert np.array_equal(result[i], expected_result[i]), f"The obtained results for point {i+1} do not match the expectations."

    def test_shift_zero(self):
        '''This test verify that by applying a shift of 0 no change in the coordinates is applied.

            Test Steps:
                -Define a value for 'shift_z' of 0.
                -Execute the "shift_slab_along_z" function, passing the coordinate list and shift value as input.
                -Verify if the returned results match the initial coordinate list.'''
        shift_z = 0
        result = functions.shift_slab_along_z(self.coord_list, shift_z)
        for i in range(len(result)):
            assert np.array_equal(result[i], self.coord_list[i]), f"The obtained results for point {i+1} do not match the expectations."

    def test_return_to_initial_conditions(self):
        '''This test verify if the application of a positive shift and the subsequent negative shift of the same 
           amount lead to the inital coordinate values.

            Test Steps:
                -Define two shifts 'first_shift_z' and 'second_shift_z' that are equal except for their signs.
                -Execute the "shift_slab_along_z" function, passing the coordinate list and the first shift value as input.
                -Execute the "shift_slab_along_z" function a second time, applying the second shift value to the coordinates just obtained.
                -Verify if the returned results match the initial coordinate list.'''
        first_shift_z = 2.0
        second_shift_z = -2.0
        first_result = functions.shift_slab_along_z(self.coord_list, first_shift_z)
        final_result = functions.shift_slab_along_z(first_result, second_shift_z)
        for i in range(len(final_result)):
            assert np.array_equal(final_result[i], self.coord_list[i]), f"The obtained results for point {i+1} do not match the expectations."


class TestDistanceBetweenHighestZvalues(unittest.TestCase):
    def test_distance_between_highest_z_values(self):
        '''This function is a unit test for the "distance_between_highest_z_values" function.

            Test Steps:
                -Define a sample set of atomic coordinates (coords).
                -Define the expected value for the distance.
                -Call the distance_between_highest_z_values function with the provided coordinates to get the result.
                -Compare the obtained result with the expected distance to validate the accuracy of the function.'''
        coords = [
            [1.0, 2.0, 5.0],
            [3.0, 4.0, 8.0],
            [5.0, 6.0, 2.0],
            [7.0, 8.0, 7.0],
        ]

        expected_distance=1.0
        result = functions.distance_between_highest_z_values(coords)

        assert result == expected_distance, "Test failed: Output doesn't match the expected result."

    def test_multiple_points_highest_z(self):

        coords = [
            [1.0, 2.0, 5.0],
            [3.0, 4.0, 8.0],
            [5.0, 6.0, 8.0],
            [7.0, 8.0, 7.0],
        ]

        expected_distance=0
        result = functions.distance_between_highest_z_values(coords)

        assert result == expected_distance, "Test failed: Output doesn't match the expected result."

    def test_distance_with_negative_z_values(self):
        coords = [
            [1.0, 2.0, -5.0],
            [3.0, 4.0, -8.0],
            [5.0, 6.0, -2.0],
            [7.0, 8.0, -7.0],
        ]
        expected_distance = abs(coords[0][2] - coords[2][2])
        result = functions.distance_between_highest_z_values(coords)
        assert result == expected_distance, "Test failed: Output doesn't match the expected result."


class CalculateShiftZCleanCase(unittest.TestCase):
    coord_bottom_slab = [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0], [0.0, 0.0, 4.0]]
    coord_upper_slab = [[0.0, 0.0, 6.0], [0.0, 0.0, 8.0], [0.0, 0.0, 10.0]]

    def test_equal_interlayer_distances(self):
        interlayer_distance = 2.0
        expected_shift_z = -4

        result = functions.calculate_shift_z_clean_case(interlayer_distance, interlayer_distance, self.coord_bottom_slab, self.coord_upper_slab)

        assert result == expected_shift_z, "Test failed: Output doesn't match the expected result."

    def test_calculate_shift_z_different_distances(self):
        interlayer_distance_bottom = 2.0
        interlayer_distance_upper = 3.0
        expected_shift_z = -3.5

        result = functions.calculate_shift_z_clean_case(interlayer_distance_bottom, interlayer_distance_upper, self.coord_bottom_slab, self.coord_upper_slab)

        assert result == expected_shift_z, "Test failed: Output doesn't match the expected result."

    def test_empty_upper_slab(self):
        interlayer_distance_bottom = 2.0
        interlayer_distance_upper = 3.0
        coord_upper_slab = [] 

        with self.assertRaises(ValueError):
            functions.calculate_shift_z_clean_case(interlayer_distance_bottom, interlayer_distance_upper, self.coord_bottom_slab, coord_upper_slab)

class CalculateShiftZDecoratedCase(unittest.TestCase):
    def test_generic_case(self):
        coord_bottom_slab = [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0], [0.0, 0.0, 4.0]]
        coord_adsorption = [[0.0, 0.0, 5.0], [0.0, 0.0, 6.0]]
        coord_upper_slab = [[0.0, 0.0, 6.0], [0.0, 0.0, 8.0]]

        expected_shift_z = -3

        result = functions.calculate_shift_z_decorated_case(coord_bottom_slab, coord_adsorption, coord_upper_slab)

        assert result == expected_shift_z, "Test failed: Output doesn't match the expected result."

    def test_empty_upper_slab(self):
        coord_bottom_slab = [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0], [0.0, 0.0, 4.0]]
        coord_adsorption = [[0.0, 0.0, 5.0]]
        coord_upper_slab = [] 

        with self.assertRaises(ValueError):
            functions.calculate_shift_z_decorated_case(coord_bottom_slab, coord_adsorption, coord_upper_slab)

    def test_empty_bottom_slab(self):
        coord_bottom_slab = []  
        coord_adsorption = [[0.0, 0.0, 5.0]]
        coord_upper_slab = [[0.0, 0.0, 6.0], [0.0, 0.0, 8.0]]

        with self.assertRaises(ValueError):
            functions.calculate_shift_z_decorated_case(coord_bottom_slab, coord_adsorption, coord_upper_slab)

    def test_empty_coord_adsorption(self):
        coord_bottom_slab = [[0.0, 0.0, 0.0], [0.0, 0.0, 2.0], [0.0, 0.0, 4.0]]
        coord_adsorption = []  
        coord_upper_slab = [[0.0, 0.0, 6.0], [0.0, 0.0, 8.0]]

        with self.assertRaises(IndexError):
                functions.calculate_shift_z_decorated_case(coord_bottom_slab, coord_adsorption, coord_upper_slab)













































