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
    integer(hid_t) :: dset_id, dset_id2
    integer(hid_t) :: group_id
    integer(hid_t) :: filespace, memspace, dataspace_id, dataspace_id2
    integer(hsize_t), dimension(1) :: dim_eval, dim_k, dim_mem
    integer(hsize_t), dimension(2) :: dim_evec
    integer(hsize_t), dimension(3) :: dim_vel

    integer :: ia_first, ja_first, iastart, jastart, iaend, jaend, ia, ja
    integer :: lrindx, lcindx, lroffset, lcoffset, rsrc, csrc
    integer(hsize_t), dimension(1) :: ipos
    integer(hsize_t), dimension(2,1) :: coord
    integer(hsize_t), dimension(3,1) :: coord_vel
    integer(size_t), parameter :: ONE_ = 1

    integer :: proc_should_write, i
    
    integer(4) :: comm_
   
    integer, external :: numroc

    double precision, allocatable, dimension(:) :: temp

    double complex, parameter :: alpha = cmplx(1,0)
    double complex, parameter :: beta = cmplx(0,0)

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
       
    write(group_name,'(I12.12)') k_indx 

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
    call h5pclose_f(dlist_id, hdf5_error)

    write(done_line,"(I0,3A)") pzheevx_vars%comp_num_eval, " eigenvalues written to ", &
                              trim(adjustl(file_name)), " on "
    call date_time_message(trim(done_line))

    ! --------------------------
    ! Store Eigenvectors
    ! --------------------------

    ! Since the local arrays contain elements distributed in block-cyclic fashion
    ! we have to reverse this distribution and select the coordinates in the 
    ! global file that correspond to the elements in the distributed array 
    ! 
    ! If j eigenvectors are computed, only the first j columns of the distributed
    ! matrix store the eigenvectors.
    ! The processes containing the relevant columns writes them out to the file.


    if (evec_comp) then
        
        rsrc = 0
        csrc = 0

        dim_evec(1) = moire%natom
        dim_evec(2) = pzheevx_vars%comp_num_evec

        if (grid%myprow.ge.evec%desca(RSRC_)) then
            ia_first = (grid%myprow - evec%desca(RSRC_))*evec%desca(MB_)+1
        else
            ia_first = (grid%myprow + (grid%nprow - evec%desca(RSRC_)))* &
                        evec%desca(MB_) +1
        endif
        if (grid%mypcol.ge.evec%desca(CSRC_)) then
            ja_first = (grid%mypcol - evec%desca(CSRC_))*evec%desca(NB_)+1
        else
            ja_first = (grid%mypcol + (grid%npcol - evec%desca(CSRC_)))* &
                      evec%desca(NB_) +1
        endif

        proc_should_write = 0

        if (ja_first.le.pzheevx_vars%comp_num_evec) then
            proc_should_write = 1
        end if
        
        call h5screate_simple_f(2, dim_evec, dataspace_id, hdf5_error)
        call h5dcreate_f(group_id,'evec_real',H5T_IEEE_F64LE,dataspace_id,dset_id,hdf5_error)
        call h5screate_simple_f(2, dim_evec, dataspace_id2, hdf5_error)
        call h5dcreate_f(group_id,'evec_imag',H5T_IEEE_F64LE,dataspace_id2,dset_id2,hdf5_error)

        if (proc_should_write==1) then
            call h5pcreate_f(H5P_DATASET_XFER_F, dlist_id, hdf5_error)
            call h5pset_dxpl_mpio_f(dlist_id, H5FD_MPIO_INDEPENDENT_F, hdf5_error)

            dim_mem(1) = evec%size_
            call h5screate_simple_f(1, dim_mem, memspace, hdf5_error)
            
            do jastart = ja_first, pzheevx_vars%comp_num_evec, grid%npcol*evec%desca(NB_)
                do iastart = ia_first, evec%desca(M_), grid%nprow*evec%desca(MB_)
                    iaend = min(evec%desca(M_), iastart+evec%desca(MB_)-1)
                    jaend = min(evec%desca(N_), jastart+evec%desca(NB_)-1)
                    ia = iastart
                    ja = jastart
                    call infog2l (ia, ja, evec%desca, grid%nprow, grid%npcol, &
                                  grid%myprow, grid%mypcol, lroffset, lcoffset, &
                                  rsrc, csrc)
                    do ja=jastart,jaend
                        do ia=iastart,iaend
                            lrindx = lroffset + (ia-iastart)
                            lcindx = lcoffset + (ja-jastart)
                            ipos(1) = lrindx + (lcindx-1)*evec%desca(LLD_)
                            
                            coord(1,1) = ia
                            coord(2,1) = ja
                            
                            if (ja.le.pzheevx_vars%comp_num_evec) then
                             call h5sselect_elements_f(memspace, H5S_SELECT_APPEND_F, 1, &
                                    ONE_, ipos, hdf5_error) 
                             call h5sselect_elements_f(dataspace_id,H5S_SELECT_APPEND_F,2,&
                                    ONE_, coord, hdf5_error)  
                             call h5sselect_elements_f(dataspace_id2,H5S_SELECT_APPEND_F,2, &
                                    ONE_, coord, hdf5_error)
                            end if
                        end do
                    end do
                end do
            end do
            
            allocate(temp(evec%size_))
            temp = real(evec%mat)
            call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, &
                            dim_evec, hdf5_error, mem_space_id = &
                            memspace, file_space_id = dataspace_id, xfer_prp=dlist_id) 
            temp = aimag(evec%mat)
            call h5dwrite_f(dset_id2, H5T_IEEE_F64LE, temp, &
                            dim_evec, hdf5_error, mem_space_id = &
                            memspace, file_space_id = dataspace_id2, xfer_prp=dlist_id)  
            deallocate(temp)
            call h5pclose_f(dlist_id, hdf5_error)
            call h5sclose_f(memspace,hdf5_error)
        end if
        call h5dclose_f(dset_id, hdf5_error)
        call h5sclose_f(dataspace_id, hdf5_error)
        call h5dclose_f(dset_id2, hdf5_error)
        call h5sclose_f(dataspace_id2, hdf5_error)
        write(done_line,"(I0,3A)") pzheevx_vars%comp_num_evec, " eigenvectors written to ", &
                                      trim(adjustl(file_name)), " on "
        call date_time_message(trim(done_line))
    end if
    

    ! ------------------------------------------------------------------------------
    !
    !   Store velocities
    !   ----------------
    !   
    !   To store the velocities, a dataspace of dimension (natom,natom) is created
    !   Then the derivatives of the hamiltonians are constructed and the velocity matrix 
    !   is computed as
    !
    !               v_{\alpha,k}^{m,n} = < \psi_{m,k} | dH/d\alpha | \psi_{n,k} > 
    !
    !   This generates a distributed velocity matrix of which (num_evec,num_evec)
    !   elements are populated.
    !
    ! --------------------------------------------------------------------------------
    
    
    if (comp_vel) then
        
        rsrc = 0
        csrc = 0

        dim_vel(3) = pzheevx_vars%comp_num_evec 
        dim_vel(2) = pzheevx_vars%comp_num_evec
        dim_vel(1) = 3

        if (grid%myprow.ge.vel%desca(RSRC_)) then
            ia_first = (grid%myprow - vel%desca(RSRC_))*vel%desca(MB_)+1
        else
            ia_first = (grid%myprow + (grid%nprow - evec%desca(RSRC_)))* &
                        evec%desca(MB_) +1
        endif
        if (grid%mypcol.ge.vel%desca(CSRC_)) then
            ja_first = (grid%mypcol - vel%desca(CSRC_))*vel%desca(NB_)+1
        else
            ja_first = (grid%mypcol + (grid%npcol - vel%desca(CSRC_)))* &
                      vel%desca(NB_) +1
        endif

        proc_should_write = 0

        if (ja_first.le.pzheevx_vars%comp_num_evec) then
            proc_should_write = 1
        end if
        
        call h5screate_simple_f(3, dim_vel, dataspace_id, hdf5_error)
        call h5dcreate_f(group_id,'vel_real',H5T_IEEE_F64LE,dataspace_id,dset_id,hdf5_error)
        call h5screate_simple_f(3, dim_vel, dataspace_id2, hdf5_error)
        call h5dcreate_f(group_id,'vel_imag',H5T_IEEE_F64LE,dataspace_id2,dset_id2,hdf5_error)

        do i=1,3
          call create_hamiltonian(k_indx,i)
          call pzgemm('N','N',moire%natom, moire%natom, moire%natom, alpha, &
                       hamiltonian%mat, 1, 1, hamiltonian%desca,            &
                       evec%mat,1,1,evec%desca,beta, vel%mat,1,1,vel%desca)

          call pzgemm('C','N',moire%natom, moire%natom, moire%natom, alpha,  &
                       evec%mat, 1, 1, evec%desca,vel%mat,1,1,vel%desca,beta,& 
                       vel%mat, 1, 1, vel%desca)
          if (proc_should_write==1) then
              call h5pcreate_f(H5P_DATASET_XFER_F, dlist_id, hdf5_error)
              call h5pset_dxpl_mpio_f(dlist_id, H5FD_MPIO_INDEPENDENT_F, hdf5_error)

              dim_mem(1) = vel%size_
              call h5screate_simple_f(1, dim_mem, memspace, hdf5_error)
              
              do jastart = ja_first, pzheevx_vars%comp_num_evec, grid%npcol*vel%desca(NB_)
                  do iastart = ia_first, pzheevx_vars%comp_num_evec, grid%nprow*vel%desca(MB_)
                      iaend = min(vel%desca(M_), iastart+vel%desca(MB_)-1)
                      jaend = min(vel%desca(N_), jastart+vel%desca(NB_)-1)
                      ia = iastart
                      ja = jastart
                      call infog2l (ia, ja, vel%desca, grid%nprow, grid%npcol, &
                                    grid%myprow, grid%mypcol, lroffset, lcoffset, &
                                    rsrc, csrc)
                      do ja=jastart,jaend
                          do ia=iastart,iaend
                              lrindx = lroffset + (ia-iastart)
                              lcindx = lcoffset + (ja-jastart)
                              ipos(1) = lrindx + (lcindx-1)*vel%desca(LLD_)
                              
                              coord_vel(3,1) = ja
                              coord_vel(2,1) = ia
                              coord_vel(1,1) = i

                              if (ja.le.pzheevx_vars%comp_num_evec) then
                               call h5sselect_elements_f(memspace, H5S_SELECT_APPEND_F, 1, &
                                      ONE_, ipos, hdf5_error) 
                               call h5sselect_elements_f(dataspace_id,H5S_SELECT_APPEND_F,3,&
                                      ONE_, coord_vel, hdf5_error)  
                               call h5sselect_elements_f(dataspace_id2,H5S_SELECT_APPEND_F,3,&
                                      ONE_, coord_vel, hdf5_error)
                              end if
                          end do
                      end do
                  end do
              end do
              
              allocate(temp(vel%size_))
              temp = real(vel%mat)
              call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, &
                              dim_vel, hdf5_error, mem_space_id = &
                              memspace, file_space_id = dataspace_id, xfer_prp=dlist_id) 
              temp = aimag(vel%mat)
              call h5dwrite_f(dset_id2, H5T_IEEE_F64LE, temp, &
                              dim_vel, hdf5_error, mem_space_id = &
                              memspace, file_space_id = dataspace_id2, xfer_prp=dlist_id)  
              deallocate(temp)
              call h5sselect_none_f(memspace, hdf5_error)
              call h5sselect_none_f(dataspace_id,hdf5_error)
              call h5sselect_none_f(dataspace_id2,hdf5_error)
              call h5pclose_f(dlist_id, hdf5_error)
              call h5sclose_f(memspace,hdf5_error)
          end if
        end do
        call h5dclose_f(dset_id, hdf5_error)
        call h5sclose_f(dataspace_id, hdf5_error)
        call h5dclose_f(dset_id2, hdf5_error)
        call h5sclose_f(dataspace_id2, hdf5_error)
        write(done_line,"(3A)") "Velocity matrix elements written to ", trim(adjustl(file_name)), " on "
        call date_time_message(trim(done_line))
    end if
    
    
    
    
    call h5gclose_f(group_id,hdf5_error)
    call h5pclose_f(glist_id, hdf5_error)
    call h5pclose_f(plist_id, hdf5_error)
    call h5fclose_f(file_id, hdf5_error)

    call mpi_barrier(mpi_global%comm, mpierr)

    
    return

end subroutine
