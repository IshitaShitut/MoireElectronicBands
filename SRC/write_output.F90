subroutine write_output(k_indx)

    use hdf5
    use global_variables
    use mpi

    implicit none
    
    integer, intent(in) :: k_indx
    character(len=char_len) :: file_name, group_name, done_line
    logical :: file_exists
    
    ! HDF5 Variables
    ! --------------

    integer(4) :: hdf5_error
    integer(hid_t) :: plist_id, glist_id, dlist_id
    integer(hid_t) :: file_id
    integer(hid_t) :: dset_id
    integer(hid_t) :: group_id
    integer(hid_t) :: filespace
    integer(hsize_t), dimension(1) :: dim_eval, dim_k
    integer(hsize_t), dimension(2) :: dim_evec

   
    double complex :: work_prnt(10000)
    integer :: i


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

    
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf5_error)
    call h5pset_fapl_mpio_f(plist_id, mpi_global%comm , MPI_INFO_NULL, hdf5_error)
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
    call h5pclose_f(dlist_id, hdf5_error)
    call h5sclose_f(filespace,hdf5_error)
    call h5dclose_f(dset_id, hdf5_error)

    write(done_line,"(I0,3A)") pzheevx_vars%comp_num_eval, " eigenvalues written to ", &
                              trim(adjustl(file_name)), " on "
    call date_time_message(trim(done_line))

    ! --------------------------
    ! Store Eigenvectors
    ! --------------------------

    ! Since the local arrays contain elements distributed in block-cyclic fashion
    ! we have to reverse this distribution and select the coordinates in the 
    ! global file that correspond to the elements in the distributed array 

    if (pzheevx_vars%comp_evec=='V') then
        call pzlaprnt(moire%natom, pzheevx_vars%comp_num_evec, evec%mat, 1, 1, &
                      evec%desca, 0, 0, 'Z', 6, work_prnt) 
        dim_evec(1) = moire%natom
        dim_evec(2) = pzheevx_vars%comp_num_evec
        call h5screate_simple_f(2, dim_evec, filespace, hdf5_error)
        call h5dcreate_f(group_id, 'eigvec.real', H5T_IEEE_F64LE, filespace, &
                         dset_id, hdf5_error)
        call h5pcreate_f(H5P_DATASET_XFER_F, dlist_id, hdf5_error)
        call h5pset_dxpl_mpio_f(dlist_id, H5D_MPIO_COLLECTIVE_F, hdf5_error)        


        allocate(evec_selection_arr(2 , evec%size_)) 

        call reverse_block_cyclic_dist()

        deallocate(evec_selection_arr)

!        write(*,*) "Writing eigenvectors"
    end if

    call h5gclose_f(group_id,hdf5_error)
    call h5pclose_f(plist_id, hdf5_error)
    call h5fclose_f(file_id, hdf5_error)

    call mpi_barrier(mpi_global%comm, mpierr)
    
    return

end subroutine





subroutine reverse_block_cyclic_dist()

    use global_variables

    implicit none

    integer :: ia_first, ja_first, iastart, jastart, iaend, jaend , rsrc, csrc
    integer :: lroffset, lcoffset, ia, ja, lrindx, lcindx, ipos

    rsrc = 0
    csrc = 0

    if (grid%myprow.ge.evec%desca(RSRC_)) then
        ia_first = (grid%myprow - evec%desca(RSRC_))*evec%desca(MB_)+1
    else
        ia_first = (grid%myprow + (grid%nprow - evec%desca(RSRC_)))* &
                    evec%desca(MB_) + 1
    endif

    if (grid%mypcol.ge.evec%desca(CSRC_)) then
        ja_first = (grid%mypcol - evec%desca(CSRC_))*evec%desca(NB_)+1
    else
        ja_first = (grid%mypcol + (grid%npcol - evec%desca(CSRC_)))* &
                    evec%desca(NB_) + 1
    endif


    do jastart = ja_first, pzheevx_vars%comp_num_evec, grid%npcol*evec%desca(NB_)
!    do jastart = ja_first, evec%desca(N_), grid%npcol*evec%desca(NB_)
        do iastart = ia_first, evec%desca(M_), grid%nprow*evec%desca(MB_)
            iaend = min(evec%desca(M_), iastart+evec%desca(MB_)-1)
            jaend = min(evec%desca(N_), jastart+evec%desca(NB_)-1)

            ia = iastart
            ja = jastart

            call infog2l (ia, ja, evec%desca, grid%nprow, grid%npcol, &
                          grid%myprow, grid%mypcol, lroffset, lcoffset, rsrc, csrc)

            do ja=jastart,jaend
                do ia=iastart,iaend
                    lrindx = lroffset + (ia-iastart)
                    lcindx = lcoffset + (ja-jastart)
                    ipos = lrindx + (lcindx-1)*evec%desca(LLD_)
                    
                    evec_selection_arr(1,ipos) = ia
                    evec_selection_arr(2,ipos) = ja
                    
                    write(*,'(6I10)') mpi_global%rank, grid%myprow, grid%mypcol, ia, ja, ipos

                end do
            end do

        end do
    end do
    
    return

end subroutine
