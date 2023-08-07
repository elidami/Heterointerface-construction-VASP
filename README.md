# Heterointerface construction VASP

# Table of content

- [Abstract](#sezione-1)
- [How to](#sezione-2)
- [The code](#sezione-3)
- [Results](#sezione-4)

## Abstract
This project aims to create the POSCAR file of clean and decorated heterogeneous interfaces, which can be used as input for VASP calculations.  
We refer to "clean" interfaces when we create a simple interface between two slabs. While, we define "decorated" interfaces the case in which we intercalate one atom between the two considered slabs.
In particular, for the clean interface case, the code is able to handle systems built between a metal fcc(111) slab and a (1x1)C(111) slab with single dangling bond terminations. While, for what regards the case of decorated interfaces, we have considered the case in which the intercalated atom is already adsorbed on the diamond surface and a metal fcc(111) slab is approached to create the interface. 


## How to
In order to start using this code, the user has to download the .zip folder and install the three python libraries on which the script is based: numpy, pandas and pymatgen. The installation of the packages can be done via pip ('pip install <package name>') or, if a conda environment is used, using the command 'conda install <package name>'.


## The code


## Results