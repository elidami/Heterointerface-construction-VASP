import functions
import numpy as np
import unittest

class PlaneParameters(unittest.TestCase):
    def test_z_coord_equal_zero(self):
        cartesian_coord_xy = [
            [1.0, 2.0, 0.0],
            [3.0, 4.0, 0.0],
            [5.0, 6.0, 0.0],
        ]
        plane_point, plane_normal = functions.calculate_plane_parameters(cartesian_coord_xy)
        expected_plane_point = np.array([0, 0, 0]) 
        expected_plane_normal = np.array([0, 0, 1])
        assert np.array_equal(plane_point, expected_plane_point), "Test failed: plane_point doesn't match the expected result."
        assert np.array_equal(plane_normal, expected_plane_normal), "Test failed: plane_normal doesn't match the expected result."

    def test_x_coord_equal_zero(self):
        cartesian_coord_yz = [
            [0.0, 1.0, 2.0],
            [0.0, 3.0, 4.0],
            [0.0, 5.0, 6.0],
        ]
        plane_point, plane_normal = functions.calculate_plane_parameters(cartesian_coord_yz)
        expected_plane_point = np.array([0, 0, 6.0])
        expected_plane_normal = np.array([0, 0, 1]) 
        assert np.array_equal(plane_point, expected_plane_point), "Test failed: plane_point doesn't match the expected result."
        assert np.array_equal(plane_normal, expected_plane_normal), "Test failed: plane_normal doesn't match the expected result."

    def test_general_case(self):
        cartesian_coord_variable_z = [
            [1.0, 2.0, 5.0],
            [3.0, 4.0, 8.0],
            [5.0, 6.0, 2.0],
            [7.0, 8.0, 7.0],
        ]
        plane_point, plane_normal = functions.calculate_plane_parameters(cartesian_coord_variable_z)
        expected_plane_point = np.array([0, 0, 8.0]) 
        expected_plane_normal = np.array([0, 0, 1])  
        assert np.array_equal(plane_point, expected_plane_point), "Test failed: plane_point doesn't match the expected result."
        assert np.array_equal(plane_normal, expected_plane_normal), "Test failed: plane_normal doesn't match the expected result."
   

class TestReflectCoords(unittest.TestCase):
    '''This class represents a unit test for the "reflect_coord" function.
       
       Attributes:
            -coord_list: list of three numpy arrays of shape (3,) representing atomic positions.
            
       Methods:
            -test_reflect_coord_wrt_x_axis
            -test_reflect_coord_wrt_y_axis
            -test_reflect_coord_wrt_y_axis'''
    coord_list = [
    np.array([1.0, 2.0, 3.0]),
    np.array([4.0, 5.0, 6.0]),
    np.array([7.0, 8.0, 9.0])]
    
    def test_reflect_coord_wrt_x_axis(self):
        '''This function tests the reflection of coordinates with respect to a plane perpendicular to the x direction
           which passes from the origin.

            Test Steps:
                -Define a plane point 'plane_point' and a plane normal 'plane_normal', which represent the plane used for reflection.
                -Execute the "reflect_coord" function, passing the coordinate list, plane point, and plane normal as input.
                -Verify if the returned results match the expected reflected coordinates.'''
        plane_point = np.array([0.0, 0.0, 0.0])
        plane_normal = np.array([1.0, 0.0, 0.0])

        result = functions.reflect_coord(self.coord_list, plane_point, plane_normal)

        expected_result = [
            np.array([-1.0, 2.0, 3.0]),
            np.array([-4.0, 5.0, 6.0]),
            np.array([-7.0, 8.0, 9.0])
        ]
        for i in range(len(result)):
            assert np.array_equal(result[i], expected_result[i]), f"The obtained results for point {i+1} do not match the expectations."

    def test_reflect_coord_wrt_y_axis(self):
        '''This function tests the reflection of coordinates with respect to a plane perpendicular to the y direction
            which passes from the origin. It is verified that in this case are the y atomic coordinates the ones that change 
            sign.

            Test Steps:
                -Define a plane point 'plane_point' and a plane normal 'plane_normal', which represent the plane used for reflection.
                -Execute the "reflect_coord" function, passing the coordinate list, plane point, and plane normal as input.
                -Verify if the returned results match the expected reflected coordinates.'''
        plane_point = np.array([0.0, 0.0, 0.0])
        plane_normal = np.array([0.0, 1.0, 0.0])
        result = functions.reflect_coord(self.coord_list, plane_point, plane_normal)
        expected_result = [
            np.array([1.0, -2.0, 3.0]),
            np.array([4.0, -5.0, 6.0]),
            np.array([7.0, -8.0, 9.0])
        ]

        for i in range(len(result)):
            assert np.array_equal(result[i], expected_result[i]), f"The obtained results for point {i+1} do not match the expectations."
    
    def test_invert_normals(self):
        '''Test to verify that changing sign in the components of the normal to the plane doesn't change the results.

            Test Steps:
                -Define a plane point 'plane_point' which represents the plane used for reflection.
                -Define two plane normals 'plane_normal_positive' and 'plane_normal_negative', that are different only for their signs.
                -Execute the "reflect_coord" function two times considering the two different plane normals.
                -Verify if the returned results are the same.'''
        plane_point = np.array([1.0, 0.0, 0.0])
        plane_normal_positive = np.array([1.0, 1.0, 1.0])
        plane_normal_negative = np.array([-1.0, -1.0, -1.0])
        result_positive = functions.reflect_coord(self.coord_list, plane_point, plane_normal_positive)
        result_negative = functions.reflect_coord(self.coord_list, plane_point, plane_normal_negative)

        for i in range(len(result_positive)):
            assert np.array_equal(result_positive[i], result_negative[i]), f"The obtained results for point {i+1} do not match the expectations."