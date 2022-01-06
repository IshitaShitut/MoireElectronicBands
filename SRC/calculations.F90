subroutine diagonalize_and_write()
    use global_variables
    use mpi
    implicit none  
    
    integer :: k_loc, allstat    

    call setup_arrays()

    allocate(hamiltonian%mat(hamiltonian%size_),stat=allstat)
    if (pzheevx_vars%comp_evec=='V') then
        allocate(evec%mat(evec%size_))
    else
        allocate(evec%mat(1))
    end if
    allocate(eval(moire%natom)) 

    do k_loc = k_file%start, k_file%finish
        call create_hamiltonian(k_loc)
        call diagonalize_hamiltonian()
        call write_output(k_loc)    
    end do  


    deallocate(hamiltonian%mat)
    deallocate(eval)
    if (allocated(evec%mat)) then
        deallocate(evec%mat)
    end if

    call mpi_barrier(mpi_global%comm, mpierr)
    return

end subroutine

subroutine setup_arrays()
    use global_variables
    implicit none
  
    integer :: rsrc, csrc, info
    integer, external :: numroc
#ifdef __DEBUG
    integer :: i
#endif
    rsrc = 0
    csrc = 0

    hamiltonian%locq = numroc(moire%natom, pzheevx_vars%mb, grid%mypcol, csrc, grid%npcol)
    hamiltonian%locq = max(hamiltonian%locq,1)
    hamiltonian%lld = numroc(moire%natom, pzheevx_vars%nb, grid%myprow, rsrc, grid%nprow)
    hamiltonian%lld = max(hamiltonian%lld,1)

    hamiltonian%size_ = hamiltonian%locq*hamiltonian%lld
    call descinit(hamiltonian%desca, moire%natom, moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, hamiltonian%lld, info)

#ifdef __DEBUG
    write(debug_str,'(A)') "\r\nHamiltonian Parameters: "
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A,I8,3(A,I0))') "Rank: ", mpi_global%rank, &
                                   "    |  Size: ", hamiltonian%size_, &
                                   " Local Rows: ", hamiltonian%lld, &
                                   " Local Cols: ", hamiltonian%locq
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)

    write(debug_str, '(A)') "\r\nHamiltoninan descriptor:"
    call debug_output(0)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A)') " ---------"
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  DTYPE_: ", hamiltonian%desca(DTYPE_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CTXT_: ", hamiltonian%desca(CTXT_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  M_: ", hamiltonian%desca(M_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  N_: ", hamiltonian%desca(N_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  MB_: ", hamiltonian%desca(MB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  NB_: ", hamiltonian%desca(NB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  RSRC_: ", hamiltonian%desca(RSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CSRC_: ", hamiltonian%desca(CSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  LLD_: ", hamiltonian%desca(LLD_)
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do


#endif

    evec%locq = numroc(moire%natom, pzheevx_vars%mb, grid%mypcol, csrc, grid%npcol)
    evec%locq = max(evec%locq,1)
    evec%lld = numroc(moire%natom, pzheevx_vars%nb, grid%myprow, rsrc, grid%nprow)
    evec%lld = max(evec%lld,1)

    evec%size_ = evec%locq*evec%lld
    call descinit(evec%desca, moire%natom, moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, evec%lld, info)    

#ifdef __DEBUG
    write(debug_str,'(A)') "\r\nEigenvector Parameters: "
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A,I8,3(A,I0))') "Rank: ", mpi_global%rank, &
                                   "    |  Size: ", evec%size_, &
                                   " Local Rows: ", evec%lld, &
                                   " Local Cols: ", evec%locq
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)

    write(debug_str, '(A)') "\r\nEigenvector descriptor:"
    call debug_output(0)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A)') " ---------"
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  DTYPE_: ", evec%desca(DTYPE_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CTXT_: ", evec%desca(CTXT_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  M_: ", evec%desca(M_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  N_: ", evec%desca(N_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  MB_: ", evec%desca(MB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  NB_: ", evec%desca(NB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  RSRC_: ", evec%desca(RSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CSRC_: ", evec%desca(CSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  LLD_: ", evec%desca(LLD_)
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do

#endif


    return

end subroutine
