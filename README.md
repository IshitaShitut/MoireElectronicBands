# Moire Electronic Bands

Code to calculate the tight binding band structure of Moire systems.  
The Moire structures should be in LAMMPS output format.  

## Prerequisites

The prerequisites to run this code are:   

>
> Parallel HDF5 (checked for v1.12.0.0 and v1.12.1)  
> Scalapack (checked with Intel MKL, AMD Scalapack)
>


## Instructions for usage

Build the code using:  
```
make -f Makefile.<dist> (supported dists = Intel MKL, AOCL)
```

While compilation the user can turn on the `-D__DEBUG` flag to get certain debugging informations.  


A sample input file would look like this:

```
# -----------------------------------------------------------------
# Input file for electronic Band Structure Calculations in Moire
# -----------------------------------------------------------------


lammps file location  : 
lammps file name      :
natom                 : 
atom types            : 	     
atom style            :		     # LAMMPS atom style used in the data file
				     # Supported styles are: atomic and molecular


k file location       : 
k file name           : 
nkpt                  :		     # Number of k points 


num kpools            : 	     # Number of k-pools --- Not implemented in v1.0


compute eigvecs       : < yes/no >
group velocity        : < yes/no >
range                 : < A/I/V >    # A = compute all evals
                                     # I = compute eigenvalues from min index to max index
                                     # V = compute eigenvalues from a min-max range 

min eigval            :              # minimum value of the eigenvalue to be found
				     # Ignored if range = A/I 
max eigval            :		     # maximum value of the eigenvalue to be found 
				     # Ignored if range = A/I 
min index             :		     # lowest index of the eigenvalue to be found
			             # Ignored if range = A/V  
max index             : 	     # highest index of the eigenvalue to be found 
				     # Ignored if range = A/V

mb                    :              # (mb, nb) are the block sizes for block cyclic
nb                    : 	     # distribution of the matrix. 

abstol		      :		     # Tolerance for eigenvalue convergence
orfac		      :	             # Tolerance for eigenvector orthogonalization


local normals         : < yes/no >   # Compute local normals. If set to no, local normals
                                     # are set along the z axis for all atoms.
num neighbours        : 	     # Number of neigbouring cells to scan for convergence
e field z             :              # in meV  
onsite energy         :              # in meV (For each type of atom, the onsite energy
                                     #         has to be specified, seperated with , )

output file name      : 
output file location  : 
```

Run the code as 
```
mpiexec -np <no. of cores> <path to executable> <path to input file> 
```

## Known Issues

K-pools are not supported yet. Will be supported in v2.0

