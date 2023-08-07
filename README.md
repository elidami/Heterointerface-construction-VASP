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


## Results
