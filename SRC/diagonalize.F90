subroutine diagonalize_hamiltonian()

    use global_variables
    use mpi

    implicit none

    integer :: ia, ja
    double precision, allocatable, dimension(:) :: gap, rwork
    double complex, allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork, ifail, iclustr
    integer :: lwork, liwork, lrwork, lgap, lifail, liclustr, info
    character(len=char_len) :: done_line

    ia = 1
    ja = 1
    lgap = grid%nprow*grid%npcol
    lifail = moire%natom
    liclustr = 2*lgap

    allocate(gap(lgap))
    allocate(ifail(lifail))
    allocate(iclustr(liclustr))
    lwork = -1
    lrwork = -1
    liwork = -1
    allocate(work(1))
    allocate(rwork(1))
    allocate(iwork(1))


    call pzheevx(pzheevx_vars%comp_evec, pzheevx_vars%range_, 'U', moire%natom,  &
    hamiltonian%mat, ia, ja, hamiltonian%desca, pzheevx_vars%vl,pzheevx_vars%vu, &
    pzheevx_vars%il, pzheevx_vars%iu, pzheevx_vars%abstol, pzheevx_vars%comp_num_eval, &
    pzheevx_vars%comp_num_evec, eval, pzheevx_vars%orfac, evec%mat, 1, 1, evec%desca, work, &
    lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)

    lwork  = int(abs(work(1)))+1
    lrwork = int(abs(rwork(1)))+1
    liwork = int(abs(iwork(1)))+1
    
    deallocate(work)
    deallocate(rwork)
    deallocate(iwork)

    allocate(work(lwork))
    allocate(rwork(lrwork))
    allocate(iwork(liwork))
  
    call pzheevx(pzheevx_vars%comp_evec, pzheevx_vars%range_, 'U', moire%natom,   &
    hamiltonian%mat, ia, ja, hamiltonian%desca, pzheevx_vars%vl, pzheevx_vars%vu, &
    pzheevx_vars%il, pzheevx_vars%iu, pzheevx_vars%abstol, pzheevx_vars%comp_num_eval, &
    pzheevx_vars%comp_num_evec,eval, pzheevx_vars%orfac, evec%mat, 1, 1, evec%desca, work, &
    lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)
    

    if (info.gt.0) then
        if (mod(info,2).ne.0) then
            write(debug_str,'(A)') "One or more eigenvalues failed to converge"
            call debug_output(info)
        end if
        if (mod(info/2, 2).ne.0) then
            write(debug_str,'(A)') "Eigenvectors corresponding to the following indices could not be orthogonalized"
            call debug_output(info)
#ifdef __KPOOL
            if (mpi_local%rank==0) then
#else
            if (mpi_global%rank==0) then
#endif
                write(*,*) iclustr
            end if
        end if
        if (mod(info/4, 2).ne.0) then
            write(debug_str,'(A)') "All eigenvectors in the specified interval could not be computed due to insufficient memory"
            call debug_output(info)
        end if
        if (mod(info/8,2).ne.0) then
            write(debug_str,'(A)') "One or more eigenvalues were not computed"
            call debug_output(info)
        end if
    end if

#ifdef __DEBUG    
#ifdef __KPOOL
    if (mpi_global%local==0) then
#else
    if (mpi_global%rank==0) then
#endif
        write(*,*) eval(1:pzheevx_vars%comp_num_eval)
    end if
#endif

    deallocate(work)
    deallocate(iwork)
    deallocate(rwork)
    deallocate(gap)
    deallocate(ifail)
    deallocate(iclustr)

#ifdef __KPOOL
    write(done_line,"(A,I0,A)") "Diagonalization completed in group : ",mpi_local%color,&
                                   " on "
    call date_time_message_local(trim(done_line))
#else
    write(done_line, "(A,I0,A)") "Diagonalization completed on "
    call date_time_message(trim(done_line))  
#endif

end subroutine
