# Moire Electronic Bands

Code to calculate the tight binding band structure of Moire systems.  
The supported Moire structures right now are in LAMMPS output format.  

## Prerequisites

The prerequisites to run this code are:   
```
HDF5 (checked for v1.12.0.0 and v1.12.1)  
Scalapack
```

## Instructions for usage

Build the code using:  
```
make -f Makefile.<dist> (supported dists = Intel MKL, AOCL)
```

While compilation the user can turn on the `-D__DEBUG` flag to get certain debugging informations.  
K-pools are not supported yet.


A sample input file would look like this:

```
# -----------------------------------------------------------------
# Input file for electronic Band Structure Calculations in Moire
# -----------------------------------------------------------------


lammps file location  : ./Examples
lammps file name      : lammps_tblg_1.12
natom                 : 10444
atom types            : 1
atom style            : molecular


k file location       : ./Examples
k file name           : k_points.dat
nkpt                  : 31


num kpools            : 3


compute eigvecs       : no
range                 : I
min eigval            : -4.25
max eigval            : 1.98
min index             : 5212
max index             : 5231

num neighbours        : 1
e_field_z             : 50         # in meV  
mb                    : 32
nb                    : 32

output file name      : bands_1.12_E10.hdf5
output_file_location  : ./Examples
```

Run the code as 
```
mpiexec -np <no. of cores> <path to executable> <path to input file> 
```
