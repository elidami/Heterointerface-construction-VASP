import functions
import numpy as np
import os
import unittest
from pymatgen.core.structure import Structure

'''
    In this file are defined the classes TestDiamondHighSymmetryPoints, TestMetalHighSymmetryPoints and TestShiftSlabOnXY.
    They contain the test methods for the functions that deal with the identification of the high symmetry points of the slabs
     and the subsequent shift on the xy plane based on these selected sites. '''

class TestDiamondHighSymmetryPoints(unittest.TestCase):
    '''This class contains unit tests for the function 'C_111_high_symmetry_points' to verify the accuracy
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
    '''This class contains unit tests for the function 'metal_fcc_111_high_symmetry_points' to verify the accuracy
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
    '''This class contains unit tests for the function 'shift_slab_on_xy'. Different limiting cases are investigated 
    in each of the following methods.

    Methods:
        -setUp: Create a temporary example file "original_file" containing the structure of a Cu(111) slab in POSCAR format.
        -test_shift_slab_on_xy: generic example for the shift.
        -test_shift_slab_positive_direction: case of a positive shift.
        -test_shift_slab_negative_direction: case of a negative shift.
        -test_no_shift_required: case of a shift equal to zero.
        - test_shift_slab_with_nonexistent_file: verify the raise of FileNotFound Error if the file doesn't exist.
        -tearDown: Clean up after the tests by removing the temporary file. '''
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
        '''This function verify that the application of a generic shift leads to the expected results. 

            Test Steps:
              -Define a generic value for the reference points for Cu and C atoms: selected_site_Cu and selected_site_C.
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
        '''This function verify that the application of a positive shift leads to a greater value for the coordinates. 
    
            Test Steps:
              -Define values for selected_site_Cu and selected_site_C so that a positive shift is obtained.
              -Get the atomic structure of the temporary POSCAR file with pymatgen tools.
              -Call the shift_slab_on_xy function with the temporary file and the selected reference points to get the result.
              -Check that the shifted coordinates are greater than the initial ones.'''

        selected_site_metal = np.array([1.0, 1.0, 0.0])
        selected_site_C = np.array([2.0, 2.0, 0.0])
        structure = Structure.from_file(self.original_file)
        result = functions.shift_slab_on_xy(self.original_file, selected_site_metal, selected_site_C)

        for coord, site in zip(result, structure):
                assert coord[0] > site.coords[0] and coord[1] > site.coords[1]

    def test_shift_slab_negative_direction(self):
        '''This function verify that the application of a negative shift leads to a smaller value for the coordinates. 

            Test Steps:
              -Define values for selected_site_Cu and selected_site_C so that a negative shift is obtained.
              -Get the atomic structure of the temporary POSCAR file with pymatgen tools.
              -Call the shift_slab_on_xy function with the temporary file and the selected reference points to get the result.
              -Check that the shifted coordinates are smaller than the initial ones.'''

        selected_site_metal = np.array([2.0, 2.0, 0.0])
        selected_site_C = np.array([1.0, 1.0, 0.0])
        structure = Structure.from_file(self.original_file)
        result = functions.shift_slab_on_xy(self.original_file, selected_site_metal, selected_site_C)

        for coord, site in zip(result, structure):
                assert coord[0] < site.coords[0] and coord[1] < site.coords[1]

    def test_no_shift_required(self):
        '''This function verify that the application of a shift equal to 0 doesn't change the coordinates. 

            Test Steps:
              -Define values for selected_site_Cu and selected_site_C so that the shift is 0.
              -Get the atomic structure of the temporary POSCAR file with pymatgen tools.
              -Call the shift_slab_on_xy function with the temporary file and the selected reference points to get the result.
              -Check that the shifted coordinates are equal to the initial ones.'''

        selected_site_metal = np.array([1.0, 1.0, 0.0])
        selected_site_C = np.array([1.0, 1.0, 0.0])
        structure = Structure.from_file(self.original_file)
        result = functions.shift_slab_on_xy(self.original_file, selected_site_metal, selected_site_C)

        for coord, site in zip(result, structure):
                assert coord[0] == site.coords[0] and coord[1] == site.coords[1]

    def test_shift_slab_with_nonexistent_file(self):
        '''Test that the function raises a FileNotFoundError for a non-existent input file.

            Test Steps:
                - Check if the input file "nonexistent_file.txt" exists and remove it if it does.
                - Define selected reference points for metal and C atoms.
                - Call the shift_slab_on_xy function with the non-existent input file and selected reference points.
                - Verify that the function raises a FileNotFoundError.'''
        if os.path.exists("nonexistent_file.txt"):
            os.remove("nonexistent_file.txt")

        selected_site_metal = np.array([1.0, 1.0, 0.0])
        selected_site_C = np.array([2.0, 2.0, 0.0])
        with self.assertRaises(FileNotFoundError):
            functions.shift_slab_on_xy("nonexistent_file.txt", selected_site_metal, selected_site_C)
