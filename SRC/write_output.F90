subroutine write_output(k_indx)

    use hdf5
    use global_variables
    use mpi

    implicit none
    
    integer, intent(in) :: k_indx
    character(10) :: output_file = 'results.h5'
    character(500) :: group_name
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
    integer(hsize_t), dimension(1) :: dim_eval

    
    ! open hdf5 interface
    ! -------------------

    call h5open_f(hdf5_error)

    comm_ = mpi_global%comm
    ! ------------------------------------ !
    ! check if file exists                 !
    ! if file exists open file in parallel !
    ! else create new file                 !
    ! ------------------------------------ !

    inquire(file=trim(adjustl(output_file)), exist=file_exists)
       
    write(group_name,'(I010)') k_indx 

    if (file_exists) then
        if (mpi_global%rank == 0) then
            call h5fopen_f(trim(adjustl(output_file)), H5F_ACC_RDWR_F, file_id, hdf5_error)
        end if
    else
        if (mpi_global%rank==0) then
            call h5fcreate_f(trim(adjustl(output_file)), H5F_ACC_TRUNC_F, file_id, hdf5_error)
        end if
    end if

    if (mpi_global%rank==0) then
        call h5gcreate_f(file_id, group_name, group_id, hdf5_error)
        call h5gclose_f(group_id,hdf5_error)
        call h5fclose_f(file_id, hdf5_error)
    end if

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf5_error)
    call h5pset_fapl_mpio_f(plist_id, comm_, MPI_INFO_NULL, hdf5_error)
    call h5fopen_f(trim(adjustl(output_file)), H5F_ACC_RDWR_F, file_id, hdf5_error, &
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
    call h5gclose_f(group_id,hdf5_error)
 

    call h5pclose_f(plist_id, hdf5_error)

    call h5fclose_f(file_id, hdf5_error)

    call mpi_barrier(mpi_global%comm, mpierr)
    
    return

end subroutine
