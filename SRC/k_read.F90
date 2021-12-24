subroutine read_k_file()
    use global_variables
    implicit none
    integer :: i, error
    character(len=char_len) :: file_name_ , temp
#ifdef __DEBUG
    integer :: j
#endif
#ifdef __KPOOL

    call distribution_length(k_file%npt, mpi_local%color-1, num_pools, &
                             k_file%start, k_file%finish)
#else
    call distribution_length(k_file%npt, 0, 1, k_file%start, k_file%finish)
#endif


#ifdef __DEBUG
    write(debug_str, '(A)') '\r\nRank      Start K read        End K read'
    call debug_output(0)
    do i=0,mpi_global%size_-1
        if (mpi_global%rank==i) then
            write(*,'(3I8)') i, k_file%start, k_file%finish
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
#endif


    write(file_name_,'(2A)') trim(adjustl(k_file%location)),trim(adjustl(k_file%name_))

    open(unit=25, file=trim(adjustl(file_name_)),action='read',iostat=error)
    
    if (error.ne.0) then
        write(err_msg, '(2A)') 'Error reading ', trim(adjustl(file_name_))
        call error_message()
        call exit
    end if

    allocate(k_file%points(k_file%finish-k_file%start+1,3))

    if (k_file%start == 1) then
        do i=k_file%start,k_file%finish
            read(25,*) k_file%points(i,1), k_file%points(i,2), k_file%points(i,3)
        end do
    else
        do i=1,k_file%start-1
            read(25,*) temp
        end do
        do i=1,k_file%finish-k_file%start+1
            read(25,*) k_file%points(i,1), k_file%points(i,2), k_file%points(i,3)
        end do
    end if
    close(unit=25)

#ifdef __DEBUG
#ifdef __KPOOL
    debug_str = '\r\n\r\nLocal K points read in each group'
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do i=1,num_pools
        if (mpi_local%color == i) then
            if (mpi_local%rank == 0) then
                do j=1,k_file%finish-k_file%start+1
                    write(*,'(I0,3F16.6)') i, k_file%points(j,1), &
                                              k_file%points(j,2), &
                                              k_file%points(j,3)
                end do
            end if
        end if
    end do
#else
    debug_str = '\r\n\r\nK points read from file'
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    if (mpi_global%rank==0) then
        do j = 1, k_file%finish-k_file%start+1
            write(*,'(I0,3F16.6)')  j, k_file%points(j,1), &
                                       k_file%points(j,2), &
                                       k_file%points(j,3)
        end do
    end if
    call mpi_barrier(mpi_global%comm, mpierr)
    
#endif
#endif


    return

end subroutine
