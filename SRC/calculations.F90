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

    rsrc = 0
    csrc = 0

    hamiltonian%locq = numroc(moire%natom, pzheevx_vars%mb, grid%mypcol, csrc, grid%npcol)
    hamiltonian%locq = max(hamiltonian%locq,1)
    hamiltonian%lld = numroc(moire%natom, pzheevx_vars%nb, grid%myprow, rsrc, grid%nprow)
    hamiltonian%lld = max(hamiltonian%lld,1)

    hamiltonian%size_ = 1+hamiltonian%locq*hamiltonian%lld

    call descinit(hamiltonian%desca, moire%natom, moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, hamiltonian%lld, info)
        
    evec%locq = numroc(moire%natom, pzheevx_vars%mb, grid%mypcol, csrc, grid%npcol)
    evec%locq = max(evec%locq,1)
    evec%lld = numroc(moire%natom, pzheevx_vars%nb, grid%myprow, rsrc, grid%nprow)
    evec%lld = max(evec%lld,1)

    evec%size_ = 1+evec%locq*evec%lld

    call descinit(evec%desca, moire%natom, moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, evec%lld, info)    

    return

end subroutine
