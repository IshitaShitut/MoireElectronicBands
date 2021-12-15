# Moire Electronic Bands

Code to calculate the tight binding band structure of Moire systems.  
The supported Moire structures right now are in LAMMPS output format.  

## Prerequisites

The prerequisites to run this code are:   

>
> Parallel HDF5 (checked for v1.12.0.0 and v1.12.1)  
> Scalapack (checked with Intel MKL)
>


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


lammps file location  : 
lammps file name      :
natom                 : 
atom types            : 
atom style            :


k file location       : 
k file name           : 
nkpt                  : 


num kpools            : 


compute eigvecs       : yes/no
range                 : A/I/V    # A = compute all evals
                                 # I = compute eigenvalues from min index to max index
                                 # V = compute eigenvalues from a min-max range 
min eigval            : 
max eigval            : 
min index             : 
max index             : 

num neighbours        : 
e_field_z             :          # in meV  
onsite energy         :          # in meV
mb                    : 
nb                    : 

output file name      : 
output_file_location  : 
```

Run the code as 
```
mpiexec -np <no. of cores> <path to executable> <path to input file> 
```
