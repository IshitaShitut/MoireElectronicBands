subroutine initialize_mpi()

    use mpi
    use global_variables

    implicit none

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
        write(debug_str,'(A,I0,A)') 'MPI initialized successfully on ', mpi_global%size_, &
                                   ' processes.'    
        call debug_output(mpierr)
#endif
    end if

    return

end subroutine
