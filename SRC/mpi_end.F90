subroutine close_mpi()

    use mpi
    use global_variables

    implicit none

    call mpi_barrier(mpi_global%comm, mpierr)
!    call date_and_time(VALUES=date_time)
!
!    write(debug_str,trim(date_format)) "\r\nCalculations completed on ", &
!                  date_time(3), "/", date_time(2), "/", date_time(1), &
!                  " at ", date_time(5), ':',date_time(6), ':',date_time(7), &
!                  " UTC: ", date_time(4)/60, ':', mod(date_time(4),60)
!    call debug_output(0)

    call date_time_message('Calculations completed on ')     

    call mpi_finalize(mpierr)
    
    if (mpierr .ne. 0) then
        debug_str = 'MPI failed to close'    
        call debug_output(mpierr)
        call exit
    end if

    return

end subroutine
