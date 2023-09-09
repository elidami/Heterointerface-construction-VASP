import functions
import numpy as np
import os
import unittest
from pymatgen.core.structure import Structure

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

class TestDiamondHighSymmetryPoints(unittest.TestCase):
    '''This class contains unit tests for the function 'C_111_high_symmetry_points'. to verify the accuracy
        of the function's output for different high symmetry point options.

        Methods:
            -setUp: Create a temporary example file "C_slab_file" containing the structure of a (1x1)C(111) slab in POSCAR format.
            -test_hollow_fcc: Test the calculation of coordinates for the "hollow_fcc" high symmetry point.
            -test_hollow_hcp: Test the calculation of coordinates for the "hollow_hcp" high symmetry point.
            -test_top: Test the calculation of coordinates for the "top" high symmetry point.
            -test_invalid_selected_site: Test that the function raises a ValueError for an invalid selected_site.
            -tearDown: Clean up after the tests by removing the temporary file. '''
    def setUp(self):
        self.C_slab_file = 'C_slab.txt'
        with open(self.C_slab_file, 'w') as f:
            f.writelines('C\n'                                      
                         '1.00000000000000\n'     
                         ' 2.5243712088359747    0.0000000000000000    0.0000000000000000\n'
                         ' 1.2621856044179873    2.1861695954339861    0.0000000000000000\n'
                         '  0.0000000000000000    0.0000000000000000   16.6987064980000000\n'
                         'C\n'
                          ' 8\n'
                         'Selective dynamics\n'
                         'Direct\n'
                         '  0.3333333332770749  0.3333333331146875  0.0034289678025159\n'
                         '  0.6666666666566030  0.6666666666867940  0.0240215642292031\n'
                         '  0.6666666666566030  0.6666666666867940  0.1241682240429344\n'
                         ' -0.0000000000000000  0.0000000000000000  0.1531807799907116\n'
                         '  0.0000000000000000 -0.0000000000000000  0.2479704408420117\n'
                         '  0.3333333332770749  0.3333333331146875  0.2769829967897892\n'
                         '  0.3333333332770749  0.3333333331146875  0.3771296566035199\n'
                         '  0.6666666666566030  0.6666666666867940  0.3977222530302071\n')
            
    def tearDown(self):
        try:
            os.remove(self.C_slab_file)
        except FileNotFoundError:
            # If the file doesn't exist don't do anything
            pass

    def test_hollow_fcc(self):
       hollow_fcc = functions.C_111_high_symmetry_points(self.C_slab_file, "hollow_fcc")
       expected_hollow_fcc = [0.  ,       0.   ,      4.14078561] 
       hollow_fcc_x, hollow_fcc_y, hollow_fcc_z = round(hollow_fcc[0], 8), round(hollow_fcc[1], 8), round(hollow_fcc[2], 8)   
       assert (
                   hollow_fcc_x == expected_hollow_fcc[0] and
                   hollow_fcc_y == expected_hollow_fcc[1] and
                   hollow_fcc_z == expected_hollow_fcc[2]
               ), "The obtained coordinates for hollow_fcc site do not match the expectations."

    def test_hollow_hcp(self):
        hollow_hcp = functions.C_111_high_symmetry_points(self.C_slab_file, "hollow_hcp")
        expected_hollow_hcp = [1.2621856,  0.7287232,  6.29757745]
        hollow_hcp_x, hollow_hcp_y, hollow_hcp_z = round(hollow_hcp[0], 8), round(hollow_hcp[1], 8), round(hollow_hcp[2], 8) 
        assert (
                    hollow_hcp_x == expected_hollow_hcp[0] and
                    hollow_hcp_y == expected_hollow_hcp[1] and
                    hollow_hcp_z == expected_hollow_hcp[2]
                ), "The obtained coordinates for hollow_hcp site do not match the expectations."
  
    def test_top(self):
        top = functions.C_111_high_symmetry_points(self.C_slab_file, "top")
        expected_top = [2.52437121, 1.4574464,  6.64144717] 
        top_x, top_y, top_z = round(top[0], 8), round(top[1], 8), round(top[2], 8)   
        assert (
                    top_x == expected_top[0] and
                    top_y == expected_top[1] and
                    top_z == expected_top[2]
                ), "The obtained coordinates for top site do not match the expectations."
        
    def test_invalid_selected_site(self):
        with self.assertRaises(ValueError):
            functions.C_111_high_symmetry_points(self.C_slab_file, "invalid_site")

class TestMetalHighSymmetryPoints(unittest.TestCase):
    '''This class contains unit tests for the function 'metal_fcc_111_high_symmetry_points'. to verify the accuracy
    of the function's output for different high symmetry point options.
    Methods:
        -setUp: Create a temporary example file "metal_slab_file" containing the structure of a Cu(111) slab in POSCAR format.
        -test_hollow_fcc: Test the calculation of coordinates for the "hollow_fcc" high symmetry point.
        -test_hollow_hcp: Test the calculation of coordinates for the "hollow_hcp" high symmetry point.
        -test_top: Test the calculation of coordinates for the "top" high symmetry point.
        -test_invalid_selected_site: Test that the function raises a ValueError for an invalid selected_site.
        -tearDown: Clean up after the tests by removing the temporary file. '''
    def setUp(self):
        self.metal_slab_file = 'metal_fcc_slab.txt'
        with open(self.metal_slab_file, 'w') as f:
             f.writelines('Cu 111\n'                                  
                         '1.00000000000000\n'      
                           '2.5699999332000001    0.0000000000000000    0.0000000000000000\n'  
                           '1.2849999666000000    2.2256852298999998    0.0000000000000000\n'  
                           '0.0000000000000000    0.0000000000000000   24.6887716150000003\n'  
                         'Cu\n'  
                           '8\n'  
                      'Selective dynamics\n'  
                      'Direct\n'  
                      '-0.0000000000000000 -0.0000000000000000  0.0007633363954040\n'
                      ' 0.3333333432502954  0.3333333433826695  0.0846976036770734\n'
                      ' 0.6666666271921713  0.6666666863160486  0.1696670099955020\n'
                      ' 0.0000000000000000  0.0000000000000000  0.2548589898860777\n'
                      ' 0.3333333432502954  0.3333333433826695  0.3400986447852040\n'
                      ' 0.6666666271921713  0.6666666863160486  0.4252906246757862\n'
                      ' 0.0000000000000000  0.0000000000000000  0.5102600153595757\n'
                      ' 0.3333333432502954  0.3333333433826695  0.5941942355753135\n')
             
    def tearDown(self):
        try:
            os.remove(self.metal_slab_file)
        except FileNotFoundError:
            # If the file doesn't exist don't do anything
            pass

    def test_hollow_fcc(self):
       hollow_fcc = functions.metal_fcc_111_high_symmetry_points(self.metal_slab_file, "hollow_fcc")
       expected_hollow_fcc = [ 2.56999986 , 1.4837902 , 10.4999031 ]
       hollow_fcc_x, hollow_fcc_y, hollow_fcc_z = round(hollow_fcc[0], 8), round(hollow_fcc[1], 8), round(hollow_fcc[2], 8)   

       assert (
               hollow_fcc_x == expected_hollow_fcc[0] and
               hollow_fcc_y == expected_hollow_fcc[1] and
               hollow_fcc_z == expected_hollow_fcc[2]
           ), "The obtained coordinates for hollow_fcc site do not match the expectations."

    def test_hollow_hcp(self):
        hollow_hcp = functions.metal_fcc_111_high_symmetry_points(self.metal_slab_file, "hollow_hcp")
        expected_hollow_hcp = [ 0.,          0.,         12.59769298] 
        hollow_hcp_x, hollow_hcp_y, hollow_hcp_z = round(hollow_hcp[0], 8), round(hollow_hcp[1], 8), round(hollow_hcp[2], 8)

        assert (
            hollow_hcp_x == expected_hollow_hcp[0] and
            hollow_hcp_y == expected_hollow_hcp[1] and
            hollow_hcp_z == expected_hollow_hcp[2]
            ), "The obtained coordinates for hollow_hcp site do not match the expectations."

    def test_top(self):   
        top = functions.metal_fcc_111_high_symmetry_points(self.metal_slab_file, "top")
        expected_top = [ 1.285 ,      0.7418951,  14.66992578]
        top_x, top_y, top_z = round(top[0], 8), round(top[1], 8), round(top[2], 8)    

        assert (
            top_x == expected_top[0] and
            top_y == expected_top[1] and
            top_z == expected_top[2]
            ), "The obtained coordinates for top site do not match the expectations."
        
    def test_invalid_selected_site(self):
        with self.assertRaises(ValueError):
            functions.metal_fcc_111_high_symmetry_points(self.metal_slab_file, "invalid_site")
    


class TestShiftSlabOnXY(unittest.TestCase):
    def setUp(self):
        self.original_file = 'POSCAR_temp'
        with open(self.original_file, 'w') as f:
            f.writelines('Cu 111\n'
             '1.0\n'
                    '2.5699999332         0.0000000000         0.0000000000\n'
                    '1.2849999666         2.2256852299         0.0000000000\n'
                    '0.0000000000         0.0000000000        20.491980720\n'
             'Cu\n'
             '6\n'
             'Selective Dynamics\n'
             'Cartesian\n'
                '0.000000000         0.000000000         0.000000000\n' 
                '1.285000005         0.741895099         2.098396222\n' 
                '2.569999857         1.483790197         4.196792443\n' 
                '0.000000000         0.000000000         6.295188277\n' 
                '1.285000005         0.741895099         8.393584886\n' 
                '2.569999857         1.483790197        10.491980720\n' )
            
    def tearDown(self):
            os.remove(self.original_file)

    def test_shift_slab_on_xy(self):
        '''This function is a unit test for the "shift_slab_on_xy" function.

            Test Steps:
              -Create a temporary file named "POSCAR_temp" containing initial atomic coordinates.
              -Define the selected reference points for Cu and C atoms: selected_site_Cu and selected_site_C.
              -Call the shift_slab_on_xy function with the temporary file and the selected reference points to get the result.
              -Define the expected output, which represents the shifted coordinates after the x and y shifts.
              -Compare the obtained result with the expected output to validate the accuracy of the function.'''

        selected_site_Cu = [0.0, 0.0, 0.0]
        selected_site_C = [0.5, 0.5, 0.0]
        result = functions.shift_slab_on_xy(self.original_file, selected_site_Cu, selected_site_C)

        expected_output = [np.array([0.5 ,       0.5 ,  0. ]),
                           np.array([1.785 ,     1.2418951 , 2.098396222]),
                           np.array([3.06999986 ,1.9837902 ,  4.19679244]),
                           np.array([0.5     ,   0.5  ,  6.29518828]),
                           np.array([1.785     , 1.2418951 , 8.39358489]),
                           np.array([3.06999986,  1.9837902 , 10.49198072])
                           ]

        for i in range(len(result)):
            assert np.allclose(result[i], expected_output[i]), f"The obtained results for point {i+1} do not match the expectations."

    def test_shift_slab_positive_direction(self):

        selected_site_metal = np.array([1.0, 1.0, 0.0])
        selected_site_C = np.array([2.0, 2.0, 0.0])
        structure = Structure.from_file(self.original_file)
        result = functions.shift_slab_on_xy(self.original_file, selected_site_metal, selected_site_C)

        for coord, site in zip(result, structure):
                assert coord[0] > site.coords[0] and coord[1] > site.coords[1]

    def test_shift_slab_negative_direction(self):

        selected_site_metal = np.array([2.0, 2.0, 0.0])
        selected_site_C = np.array([1.0, 1.0, 0.0])
        structure = Structure.from_file(self.original_file)
        result = functions.shift_slab_on_xy(self.original_file, selected_site_metal, selected_site_C)

        for coord, site in zip(result, structure):
                assert coord[0] < site.coords[0] and coord[1] < site.coords[1]

    def test_no_shift_required(self):

        selected_site_metal = np.array([1.0, 1.0, 0.0])
        selected_site_C = np.array([1.0, 1.0, 0.0])
        structure = Structure.from_file(self.original_file)
        result = functions.shift_slab_on_xy(self.original_file, selected_site_metal, selected_site_C)

        for coord, site in zip(result, structure):
                assert coord[0] == site.coords[0] and coord[1] == site.coords[1]

    def test_shift_slab_with_nonexistent_file(self):
        # Remove input file if it exists
        if os.path.exists("nonexistent_file.txt"):
            os.remove("nonexistent_file.txt")

        # Test that the function raises a FileNotFoundError for a non-existent input file.
        selected_site_metal = np.array([1.0, 1.0, 0.0])
        selected_site_C = np.array([2.0, 2.0, 0.0])
        with self.assertRaises(FileNotFoundError):
            functions.shift_slab_on_xy("nonexistent_file.txt", selected_site_metal, selected_site_C)



def test_distance_between_highest_z_values():
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

def test_write_coords():
    '''This function is a unit test for the "write_coords" function.

    Test Steps:
        - Define sample input data, including a list of NumPy arrays representing atomic coordinates (coords)
          and boolean variables 'x_relax', 'y_relax', and 'z_relax'.
        - Call the "write_coords" function with the provided input data to get the result.
        - Compare the obtained result with the expected output to validate the accuracy of the function.'''
    coords = [
        np.array([2.5035679999999956 ,  1.4454350000000022,   9.6615339926441077 ]),
        np.array([1.2519332746767966 ,  0.7228036064245595,   0.0583534541551990]),
        np.array([5.0071359999999903,   1.4458782853606842  , 0.4115548962189028])
    ]
    x_relax = True
    y_relax = False
    z_relax = True

    result = functions.write_coords(coords, x_relax, y_relax, z_relax)

    expected_output = ['2.5035679999999956   1.4454350000000022   9.6615339926441077    T  F  T\n', 
                       '1.2519332746767966   0.7228036064245595   0.0583534541551990    T  F  T\n', 
                       '5.0071359999999903   1.4458782853606842   0.4115548962189028    T  F  T\n']

    assert result == expected_output, "Test failed: Output doesn't match the expected result."









































