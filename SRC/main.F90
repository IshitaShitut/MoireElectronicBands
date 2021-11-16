program compute_electronic_bands

   implicit none

   call initialize_mpi()
   call read_input()
   call blacs_grid_initialization() 
   call read_lammps_data()
   call read_k_file()
   call compute_local_normals()
   call diagonalize_and_write()
   call close_mpi()

end program
