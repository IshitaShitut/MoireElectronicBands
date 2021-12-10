subroutine write_output(k_indx)

    use hdf5
    use global_variables
    use mpi

    implicit none
    
    integer, intent(in) :: k_indx
    character(len=char_len) :: file_name, group_name
    logical :: file_exists
    
    ! HDF5 Variables
    ! --------------

    integer(4) :: hdf5_error
    integer(4) :: comm_
    integer(hid_t) :: plist_id, glist_id, dlist_id
    integer(hid_t) :: file_id
    integer(hid_t) :: dset_id
    integer(hid_t) :: group_id
    integer(hid_t) :: filespace
    integer(hsize_t), dimension(1) :: dim_eval, dim_k
    integer(hsize_t), dimension(2) :: dim_evec
    
    ! open hdf5 interface
    ! -------------------

    call h5open_f(hdf5_error)

    ! ------------------------------------ !
    ! check if file exists                 !
    ! if file exists open file in parallel !
    ! else create new file                 !
    ! ------------------------------------ !

    write(file_name, '(2A)') trim(adjustl(output_file_location)), &
                             trim(adjustl(output_file_name))

    inquire(file=file_name, exist=file_exists)
       
    write(group_name,'(I10)') k_indx 

    if (file_exists) then
        if (mpi_global%rank == 0) then
            call h5fopen_f(trim(adjustl(file_name)), H5F_ACC_RDWR_F, file_id, hdf5_error)
        end if
    else
        if (mpi_global%rank==0) then
            call h5fcreate_f(trim(adjustl(file_name)), H5F_ACC_TRUNC_F, file_id, hdf5_error)
        end if
    end if

    if (mpi_global%rank==0) then
        call h5gcreate_f(file_id, group_name, group_id, hdf5_error)
        dim_k(1) = 3
        call h5screate_simple_f(1, dim_k, filespace, hdf5_error)
        call h5dcreate_f(group_id, 'k_vec', H5T_IEEE_F64LE, filespace, dset_id, hdf5_error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, k_file%points(k_indx,:),dim_k,hdf5_error)
        call h5dclose_f(dset_id,hdf5_error)
        call h5sclose_f(filespace, hdf5_error)
        call h5gclose_f(group_id,hdf5_error)
        call h5fclose_f(file_id, hdf5_error)
    end if

    
    comm_ = mpi_global%comm

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf5_error)
    call h5pset_fapl_mpio_f(plist_id, comm_ , MPI_INFO_NULL, hdf5_error)
    call h5fopen_f(trim(adjustl(file_name)), H5F_ACC_RDWR_F, file_id, hdf5_error, &
                   access_prp = plist_id)    

    call h5pcreate_f(H5P_GROUP_ACCESS_F, glist_id, hdf5_error)
    call h5pset_all_coll_metadata_ops_f(glist_id, .true. , hdf5_error)
    call h5gopen_f(file_id, group_name, group_id, hdf5_error, gapl_id = glist_id)


    ! -----------------------
    ! Store eigenvalues
    ! -----------------------

    dim_eval(1) = pzheevx_vars%comp_num_eval

    call h5screate_simple_f(1, dim_eval, filespace, hdf5_error)
    call h5dcreate_f(group_id, 'eigenvalues', H5T_IEEE_F64LE, filespace, dset_id, hdf5_error)

    call h5pcreate_f(H5P_DATASET_XFER_F, dlist_id, hdf5_error)
    call h5pset_dxpl_mpio_f(dlist_id, H5FD_MPIO_COLLECTIVE_F, hdf5_error)

    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eval, dim_eval, hdf5_error, xfer_prp = dlist_id)
    
    call h5sclose_f(filespace,hdf5_error)
    call h5dclose_f(dset_id, hdf5_error)

    ! --------------------------
    ! Store Eigenvectors
    ! --------------------------

    ! Since the local arrays contain elements distributed in block-cyclic fashion
    ! we have to reverse this distribution and select the coordinates in the 
    ! global file that correspond to the elements in the  

!    if (pzheevx_vars%comp_evec=='V') then
!        dim_evec(1) = moire%natom
!        dim_evec(2) = pzheevx_vars%comp_num_evec
!        call h5screate_simple_f(2, dim_evec, filespace, hdf5_error)
!        write(*,*) "Writing eigenvectors"
!    end if

    call h5gclose_f(group_id,hdf5_error)
    call h5pclose_f(plist_id, hdf5_error)
    call h5fclose_f(file_id, hdf5_error)

    call mpi_barrier(mpi_global%comm, mpierr)
    
    return

end subroutine
