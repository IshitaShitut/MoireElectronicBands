subroutine initialize_mpi()

    use mpi
    use global_variables

    implicit none
    integer :: k

    call mpi_init(mpierr)
    mpi_global%comm = MPI_COMM_WORLD
    call mpi_comm_size(mpi_global%comm, mpi_global%size_, mpierr)
    call mpi_comm_rank(mpi_global%comm, mpi_global%rank, mpierr)

    call date_time_message("Program started on")

    if (mpierr .ne. 0) then
        debug_str = 'MPI failed to initialize'    
        call debug_output(mpierr)
        call exit
#ifdef __DEBUG
    else
        debug_str = 'MPI initialized successfully'    
        call debug_output(mpierr)
        debug_str = 'MPI Rank          MPI Size'
        call debug_output(mpierr)
        call mpi_barrier(mpi_global%comm, mpierr)   
        do k = 0,mpi_global%size_-1
            if (mpi_global%rank==k) then
                write(6,'( I6,12X,I6 )') mpi_global%rank, mpi_global%size_
            end if
            call mpi_barrier(mpi_global%comm, mpierr)
        end do
        call mpi_barrier(mpi_global%comm, mpierr)
#endif
    end if

    return

end subroutine
