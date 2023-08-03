import functions
import numpy as np
import os

def test_extract_lattice_vectors():
    '''This test ensures that the "extract_lattice_vectors" function correctly reads the lattice 
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

def test_extract_atomic_coordinates():
    '''This function is a unit test for the "extract_atomic_coordinates" function.

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

def test_direct_to_cartesian_coord():
    '''This function is a unit test for the "direct_to_cartesian_coord" function.

        Test Steps:
            -Define the lattice vectors 'a', 'b', and 'c' representing the lattice of the structure.
            -Define the direct coordinates 'direct_coord' as a numpy array containing two sets of direct coordinates.
            -Execute the "direct_to_cartesian_coord" function, passing the lattice vectors and direct coordinates as input.
            -Verify if the returned results match the expected Cartesian coordinates.'''
    a = np.array([1.00, 0.00, 0.00])
    b = np.array([0.00, 2.00, 0.00])
    c = np.array([0.00, 0.00, 3.00])

    direct_coord = np.array([[0.25, 0.50, 0.75], [0.10, 0.20, 0.30]])

    result = functions.direct_to_cartesian_coord(a, b, c, direct_coord)

    expected_result = [
        np.array([0.25, 1, 2.25]),
        np.array([0.10, 0.40, 0.90])
    ]
    for i in range(len(result)):
        assert np.allclose(result[i], expected_result[i]), f"The obtained results for coordinate {i+1} do not match the expectations."

def test_reflect_coord():
    '''This function is a unit test for the "reflect_coord" function.

        Test Steps:
            -Define a list of coordinate points 'coord_list', representing atomic positions.
            -Define a plane point 'plane_point' and a plane normal 'plane_normal', which represent the plane used for reflection.
            -Execute the "reflect_coord" function, passing the coordinate list, plane point, and plane normal as input.
            -Verify if the returned results match the expected reflected coordinates.'''
    coord_list = [
        np.array([1.0, 2.0, 3.0]),
        np.array([4.0, 5.0, 6.0]),
        np.array([7.0, 8.0, 9.0])
    ]
    plane_point = np.array([0.0, 0.0, 0.0])
    plane_normal = np.array([1.0, 0.0, 0.0])

    result = functions.reflect_coord(coord_list, plane_point, plane_normal)

    expected_result = [
        np.array([-1.0, 2.0, 3.0]),
        np.array([-4.0, 5.0, 6.0]),
        np.array([-7.0, 8.0, 9.0])
    ]

    for i in range(len(result)):
        assert np.array_equal(result[i], expected_result[i]), f"The obtained results for point {i+1} do not match the expectations."

def test_shift_slab_along_z():
    '''This function is a unit test for the "shift_slab_along_z" function.

        Test Steps:
            -Define a list of coordinate points 'coord_list', representing atomic positions.
            -Define a value 'shift_z', representing the amount of displacement in the z-direction.
            -Execute the "shift_slab_along_z" function, passing the coordinate list and shift value as input.
            -Verify if the returned results match the expected coordinates after the z-direction shift.'''
    coord_list = [
        np.array([1.0, 2.0, 3.0]),
        np.array([4.0, 5.0, 6.0]),
        np.array([7.0, 8.0, 9.0])
    ]
    shift_z = 2.0

    result = functions.shift_slab_along_z(coord_list, shift_z)

    expected_result = [
        np.array([1.0, 2.0, 5.0]),
        np.array([4.0, 5.0, 8.0]),
        np.array([7.0, 8.0, 11.0])
    ]

    for i in range(len(result)):
        assert np.array_equal(result[i], expected_result[i]), f"The obtained results for point {i+1} do not match the expectations."

def test_C_111_high_symmetry_points():
   '''This function is a unit test for the "C_111_high_symmetry_points" function.

      Test Steps:
        -Create a temporary example file "C_slab_file" containing the structure of a (1x1)C(111) slab in POSCAR format.
        -Call the "C_111_high_symmetry_points" function three times, each with different high symmetry points:
          "hollow_fcc", "hollow_hcp", and "top".
        -Extract the x, y, and z coordinates from the returned results.
        -Check if the obtained coordinates match the expected values for each high symmetry point.'''
   C_slab_file = 'C_slab.txt'
   with open(C_slab_file, 'w') as f:
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
   hollow_fcc = functions.C_111_high_symmetry_points(C_slab_file, "hollow_fcc")
   hollow_hcp = functions.C_111_high_symmetry_points(C_slab_file, "hollow_hcp")
   top = functions.C_111_high_symmetry_points(C_slab_file, "top")
   
   expected_hollow_fcc = [0.  ,       0.   ,      4.14078561] 
   expected_hollow_hcp = [1.2621856,  0.7287232,  6.29757745]
   expected_top = [2.52437121, 1.4574464,  6.64144717]
 
   hollow_fcc_x, hollow_fcc_y, hollow_fcc_z = round(hollow_fcc[0], 8), round(hollow_fcc[1], 8), round(hollow_fcc[2], 8)   
   hollow_hcp_x, hollow_hcp_y, hollow_hcp_z = round(hollow_hcp[0], 8), round(hollow_hcp[1], 8), round(hollow_hcp[2], 8) 
   top_x, top_y, top_z = round(top[0], 8), round(top[1], 8), round(top[2], 8)    
   
   assert hollow_fcc_x == expected_hollow_fcc[0], "The obtained x coordinate for hollow_fcc site does not match the expectations."
   assert hollow_fcc_y == expected_hollow_fcc[1], "The obtained y coordinate for hollow_fcc site does not match the expectations."
   assert hollow_fcc_z == expected_hollow_fcc[2], "The obtained z coordinate for hollow_fcc site does not match the expectations."

   assert hollow_hcp_x == expected_hollow_hcp[0], "The obtained x coordinate for hollow_hcp site does not match the expectations."
   assert hollow_hcp_y == expected_hollow_hcp[1], "The obtained y coordinate for hollow_hcp site does not match the expectations."
   assert hollow_hcp_z == expected_hollow_hcp[2], "The obtained z coordinate for hollow_hcp site does not match the expectations."

   assert top_x == expected_top[0], "The obtained x coordinate for top site does not match the expectations."
   assert top_y == expected_top[1], "The obtained y coordinate for top site does not match the expectations."
   assert top_z == expected_top[2], "The obtained z coordinate for top site does not match the expectations."
   
   os.remove(C_slab_file)

def test_metal_fcc_111_high_symmetry_points():
   '''This function is a unit test for the "metal_fcc_111_high_symmetry_points" function.
   
      Test Steps:
        -Create a temporary example file "metal_fcc_slab_file" containing the structure of a Cu(111) slab in POSCAR format.
        -Call the "metal_fcc_111_high_symmetry_points" function three times, each with different high symmetry points:
          "hollow_fcc", "hollow_hcp", and "top".
        -Extract the x, y, and z coordinates from the returned results.
        -Check if the obtained coordinates match the expected values for each high symmetry point.'''
   metal_slab_file = 'metal_fcc_slab.txt'
   with open(metal_slab_file, 'w') as f:
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
                    
                    
   hollow_fcc = functions.metal_fcc_111_high_symmetry_points(metal_slab_file, "hollow_fcc")
   hollow_hcp = functions.metal_fcc_111_high_symmetry_points(metal_slab_file, "hollow_hcp")
   top = functions.metal_fcc_111_high_symmetry_points(metal_slab_file, "top")
   
   expected_hollow_fcc = [ 2.56999986 , 1.4837902 , 10.4999031 ]
   expected_hollow_hcp = [ 0.,          0.,         12.59769298]
   expected_top = [ 1.285 ,      0.7418951,  14.66992578]
 
   hollow_fcc_x, hollow_fcc_y, hollow_fcc_z = round(hollow_fcc[0], 8), round(hollow_fcc[1], 8), round(hollow_fcc[2], 8)   
   hollow_hcp_x, hollow_hcp_y, hollow_hcp_z = round(hollow_hcp[0], 8), round(hollow_hcp[1], 8), round(hollow_hcp[2], 8) 
   top_x, top_y, top_z = round(top[0], 8), round(top[1], 8), round(top[2], 8)    
   
   assert hollow_fcc_x == expected_hollow_fcc[0], "The obtained x coordinate for hollow_fcc site does not match the expectations."
   assert hollow_fcc_y == expected_hollow_fcc[1], "The obtained y coordinate for hollow_fcc site does not match the expectations."
   assert hollow_fcc_z == expected_hollow_fcc[2], "The obtained z coordinate for hollow_fcc site does not match the expectations."

   assert hollow_hcp_x == expected_hollow_hcp[0], "The obtained x coordinate for hollow_hcp site does not match the expectations."
   assert hollow_hcp_y == expected_hollow_hcp[1], "The obtained y coordinate for hollow_hcp site does not match the expectations."
   assert hollow_hcp_z == expected_hollow_hcp[2], "The obtained z coordinate for hollow_hcp site does not match the expectations."

   assert top_x == expected_top[0], "The obtained x coordinate for top site does not match the expectations."
   assert top_y == expected_top[1], "The obtained y coordinate for top site does not match the expectations."
   assert top_z == expected_top[2], "The obtained z coordinate for top site does not match the expectations."

   os.remove(metal_slab_file)

def test_shift_slab_on_xy():
    '''This function is a unit test for the "shift_slab_on_xy" function.

        Test Steps:
          -Create a temporary file named "POSCAR_temp" containing initial atomic coordinates.
          -Define the selected reference points for Cu and C atoms: selected_site_Cu and selected_site_C.
          -Call the shift_slab_on_xy function with the temporary file and the selected reference points to get the result.
          -Define the expected output, which represents the shifted coordinates after the x and y shifts.
          -Compare the obtained result with the expected output to validate the accuracy of the function.'''
    original_file = "POSCAR_temp"
    with open(original_file, "w") as f:
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

    selected_site_Cu = [0.0, 0.0, 0.0]
    selected_site_C = [0.5, 0.5, 0.0]
    result = functions.shift_slab_on_xy(original_file, selected_site_Cu, selected_site_C)

    expected_output = [np.array([0.5 ,       0.5 ,  0. ]),
                       np.array([1.785 ,     1.2418951 , 2.098396222]),
                       np.array([3.06999986 ,1.9837902 ,  4.19679244]),
                       np.array([0.5     ,   0.5  ,  6.29518828]),
                       np.array([1.785     , 1.2418951 , 8.39358489]),
                       np.array([3.06999986,  1.9837902 , 10.49198072])
                       ]

    for i in range(len(result)):
        assert np.allclose(result[i], expected_output[i]), f"The obtained results for point {i+1} do not match the expectations."

    os.remove(original_file)

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

test_extract_lattice_vectors()
test_extract_atomic_coordinates()
test_direct_to_cartesian_coord()
test_reflect_coord()
test_shift_slab_along_z()
test_C_111_high_symmetry_points()
test_metal_fcc_111_high_symmetry_points()
test_shift_slab_on_xy()
test_distance_between_highest_z_values()
test_write_coords()









































