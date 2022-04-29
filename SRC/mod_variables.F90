module global_variables
    
    use hdf5

    integer, parameter :: RegInt_8 = selected_int_kind(8)
    integer, parameter :: RegInt_4 = selected_int_kind(4)

    integer , parameter :: char_len=2500
    character(len=char_len) :: err_msg, debug_str

    type blacs_info
#ifdef __KPOOL
        integer , allocatable, dimension(:) :: context
#else
        integer  :: context
#endif
        integer  :: nprow, npcol, myprow, mypcol, rank, size_
    end type blacs_info

    type lammps
        character(len=char_len) :: location, name_, atom_style
        integer  :: at_types
    end type lammps

    type bz_points
        character(len=char_len) :: location, name_
        integer  :: npt, start, finish
        double precision, allocatable, dimension(:,:) :: points
    end type bz_points

    type scalapack_variables
        integer  :: mb, nb, il, iu, comp_num_eval, comp_num_evec
        double precision :: vl, vu, abstol, orfac
        character(1) :: range_, comp_evec
    end type scalapack_variables

    type system
        integer  :: natom
        double precision, allocatable, dimension(:) :: mass
        integer, allocatable, dimension(:,:) :: sup_cell_info
        logical :: compute_normal
        double precision, allocatable, dimension(:,:) :: real_pos, crys, normal
        integer , allocatable, dimension(:) :: at_types_i
        integer , allocatable, dimension(:) :: lay_types
        double precision, dimension(3,3) :: lat , rec_lat
        double precision, allocatable, dimension(:) :: onsite_en(:)
    end type

    type complex_arrays
        double complex, allocatable, dimension(:) :: mat
        integer , dimension(9) :: desca
        integer :: size_, lld, locq
    end type complex_arrays





    type(blacs_info) :: grid
    type(lammps) :: lammps_file
    type(bz_points) :: k_file
    type(system) :: moire
    type(complex_arrays) :: hamiltonian, evec, vel
    logical :: evec_comp, comp_vel
    double precision, allocatable, dimension(:) :: eval
    type(scalapack_variables) :: pzheevx_vars
    integer , parameter :: BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,     &
                          CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,    &
                          RSRC_ = 7, CSRC_ = 8, LLD_ = 9
    integer  :: no_neigh
    double precision :: E_field
    


    type mpi_group
        integer  :: comm
        integer  :: color
        integer  :: key
        integer  :: size_
        integer  :: rank
    end type mpi_group

    type(mpi_group) :: mpi_global
    integer  :: mpierr

#ifdef __KPOOL
    type(mpi_group) :: mpi_local
    integer  :: num_pools
#endif

    character(500) :: output_file_name, output_file_location

    integer , dimension(8) :: date_time
    character(len=500), parameter :: date_format = '(A,X,2(I0,A),I4,3(A,I2.2),A,SP,I0,A,SS,I0)'

end module
