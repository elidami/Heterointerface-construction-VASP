# Heterointerface construction VASP

## Table of content

- [Abstract](#abstract)
- [How to](#how-to)
- [The code](#the-code)
- [Results](#results)

## Abstract
This project aims to create the POSCAR file of clean and decorated heterogeneous interfaces from the geometries of the two relaxed slabs. The created file can then be used as input for VASP (*Vienna Ab initio Simulation Package*) relaxation calculations.  
We refer to "clean" interfaces when we create a simple interface between two slabs. While, we define "decorated" interfaces the case in which we intercalate one atom between the two considered slabs.
In particular, for the clean interface case, the code is able to handle systems built between a metal fcc(111) slab and a (1x1)C(111) slab with single dangling bond terminations. While, for what regards the case of decorated interfaces, we have considered the case in which the intercalated atom is already adsorbed on the diamond surface and a metal fcc(111) slab is approached to create the interface. 


## How to
In order to start using this code, the user has to first download the .zip folder and install the three python libraries on which the script is based: `numpy`, `pandas` and `pymatgen`. The installation of the packages can be done via pip (`pip install <package name>`) or, if a conda environment is used, using the command `conda install <package name>`.  

If the user wants to create the POSCAR file of a **CLEAN INTERFACE** he/she has to follow these steps:  
1. Copy and paste the CONTCAR files of the relaxed upper and bottom slabs in the *./clean_interface_files/upper_slab.txt* and *./clean_interface_files/bottom_slab.txt*, respectively. It is necessary that the two CONTCAR files have the same lattice vectors, which will be the ones that describe the interface supercell.
2. The user has to define the input variables written in the *configuration.txt* file:
   - `selected_site_metal`  and `selected_site_C` representing the high symmetry points of metal and diamond slabs, respectively. In this way the geometry of the interface is determined: the code will perform a shift on the xy plane so that this two sites will be aligned one on top of the other. The two variables can take three possible strings: `top`, `hollow_hcp`, `hollow_fcc`.
   - `interlayer_distance_bottom_slab` and `interlayer_distance_upper_slab`: float values representing the interlayer distance between atoms in the bottom and upper slabs, respectively.
   - `x_relax`, `y_relax` and `z_relax` can take `true` or `false` values, representing relaxation options for interface optimization calculations, performed with VASP.
3. Run the script with the command `python interface.py configuration.txt`.

If the user wants to create the POSCAR file of a **DECORATED INTERFACE** he/she has to follow these steps:  
1. Copy and paste the CONTCAR files of the relaxed upper metal slab and of the bottom diamond slab with the already adsorbed atom in the *./decorated_interface_files/upper_slab.txt* and *./decorated_interface_files/bottom_slab_with_adatom.txt*, respectively. It is necessary that the two CONTCAR files have the same lattice vectors, which will be the ones that describe the interface supercell.
Copy and paste the CONTCAR file of the considered atom adsorption on the upper slab in the *./decorated_interface_files/adsorption_on_upper_slab.txt*. In this case the lattice vectors does not need to be the same as the ones describing the interface supercell.  
2. The user has to define the input variables written in the *configuration.txt* file:
   - `selected_site_metal` representing the high symmetry points of the metal slab. The code will perform a shift on the xy plane so that the intercalated atom is positioned on the metal selected site. The variable can take three possible strings: `top`, `hollow_hcp`, `hollow_fcc`.
   - `x_relax`, `y_relax` and `z_relax` can take `true` or `false` values, representing relaxation options for interface optimization calculations, performed with VASP.
3. Run the script with the command `python decorated_interface.py configuration.txt`.

In both cases, an output file named POSCAR is created in the main folder of the code. 


## The code
The project is divided in four parts:  
1. The file *interface.py* is the script that creates the POSCAR of the clean interface and is related to the input files contained in the folder *clean_interface_files*. First, the script extracts the lattice vectors and the atomic coordinates from the two input files. Then, it takes the two selected high symmetry sites on C and metal slabs and perform a shift of the upper slab on the xy plane, so that the two reference points are aligned. The upper slab is reflected with respect to its highest atom along z direction and is *z*-shifted so that the distance between the two slabs is the average of their interlayer distances.  
Finally, the lattice vectors of the supercell and the atomic coordinates that constitue the interface system are written in a output file in POSCAR format, that can be used directly for VASP calculations. 
3. The file *decorated_interface.py* is the script that creates the POSCAR of the decorated interface and is related to the input files contained in the folder *decorated_interface_files*. It works in a similar way to the previous case, with the difference that now only the high symmetry site of the metal slab needs to be selected and the *xy* shift is performed so that the metal reference site is aligned to the atom adsorbed on the diamond surface. Moreover, the shift of the upper metal slab along the *z* direction is determined by the distance that the atom has from the surface atomic plane when it is adsorbed on the upper slab. This is why the user has to provide also the file *./decorated_interface_files/adsorption_on_upper_slab.txt*.
4. The files *configuration.txt* and *functions.py* are in common for both clean and decorated interface calculations:
   - *configuration.txt* contains the input variables that needs to be specified by the user.
   - *functions.py* contains all the methods on which the two main scripts are based.
5. The file *testing.py* contains the tests of the methods written in *functions.py*.

## Results
The output of the code is a text file called POSCAR, written in a readable format for VASP.  
Here is an example:
