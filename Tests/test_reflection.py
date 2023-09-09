import functions
import numpy as np
import unittest

'''
    In this file are defined the classes PlaneParameters and TestReflectCoords
    that contain the test methods for the functions that deal with reflection operations. '''


class PlaneParameters(unittest.TestCase):
    '''This class represents a unit test for the "calculate_plane_parameters" function. Different cases
       are investigated in each of the following methods.
   
        Methods:
             -test_z_coord_equal_zero: case with all z coordinates equal to zero.
             -test_reflect_coord_wrt_y_axis: reflection wrt plane perpendicular to x direction that passes from the origin.
             -test_invert_normals: verify that changing sign in plane normal coordinates don't change the results. '''
    def test_z_coord_equal_zero(self):
        '''Test to verify that the plane point has z coordinate equal to zero when the list of all the z atomic coordinates
           are equal to zero.
           
             Test Steps:
                 -Define a list of three numpy arrays of shape (3,) representing atomic positions with the z coordinates equal to zero.
                 -Execute the "calculate_plane_parameters" function, passing the coordinate list.
                 -Verify if the returned results match the expected plane point.'''
        cartesian_coord_xy = [
            [1.0, 2.0, 0.0],
            [3.0, 4.0, 0.0],
            [5.0, 6.0, 0.0],
        ]
        plane_point, plane_normal = functions.calculate_plane_parameters(cartesian_coord_xy)
        expected_plane_point = np.array([0, 0, 0]) 
        assert np.array_equal(plane_point, expected_plane_point), "Test failed: plane_point doesn't match the expected result."

    def test_x_coord_equal_zero(self):
        '''Test to verify that the plane_point variable is not affected by the fact that x atomic coordinates are equal to zero.
        
          Test Steps:
              -Define a list of three numpy arrays of shape (3,) representing atomic positions with the z coordinates equal to zero.
              -Execute the "calculate_plane_parameters" function, passing the coordinate list.
              -Verify if the returned results match the expected plane point.'''
        cartesian_coord_yz = [
            [0.0, 1.0, 2.0],
            [0.0, 3.0, 4.0],
            [0.0, 5.0, 6.0],
        ]
        plane_point, plane_normal = functions.calculate_plane_parameters(cartesian_coord_yz)
        expected_plane_point = np.array([0, 0, 6.0])
        assert np.array_equal(plane_point, expected_plane_point), "Test failed: plane_point doesn't match the expected result."

    def test_general_case(self):
        '''Test for a general case in which atomic coordinates are different from zero.

        Test Steps:
            -Define a list of three numpy arrays of shape (3,) representing random atomic positions.
            -Execute the "calculate_plane_parameters" function, passing the coordinate list.
            -Verify if the returned results match the expected plane point.'''
        cartesian_coord_variable_z = [
            [1.0, 2.0, 5.0],
            [3.0, 4.0, 8.0],
            [5.0, 6.0, 2.0],
            [7.0, 8.0, 7.0],
        ]
        plane_point, plane_normal = functions.calculate_plane_parameters(cartesian_coord_variable_z)
        expected_plane_point = np.array([0, 0, 8.0])  
        assert np.array_equal(plane_point, expected_plane_point), "Test failed: plane_point doesn't match the expected result."
   

class TestReflectCoords(unittest.TestCase):
    '''This class represents a unit test for the "reflect_coord" function.
       
       Attributes:
            -coord_list: list of three numpy arrays of shape (3,) representing atomic positions.
            
       Methods:
            -test_reflect_coord_wrt_x_axis: reflection wrt plane perpendicular to x direction that passes from the origin.
            -test_reflect_coord_wrt_y_axis: reflection wrt plane perpendicular to x direction that passes from the origin.
            -test_invert_normals: verify that changing sign in plane normal coordinates don't change the results. '''
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