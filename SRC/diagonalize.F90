subroutine diagonalize_hamiltonian()

    use global_variables
    use mpi

    implicit none

    integer :: ia, ja
    double precision, allocatable, dimension(:) :: gap, rwork
    double complex, allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork, ifail, iclustr
    integer :: lwork, liwork, lrwork, lgap, lifail, liclustr, info
    integer, parameter :: abstol = -1
    integer, parameter :: orfac = -1
    integer :: anb, sqnpc, nps, nhetrd_lwork
    integer, external :: numroc, iceil, pjlaenv
    integer :: nn, neig, np0, mq0
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
    pzheevx_vars%il, pzheevx_vars%iu, abstol, pzheevx_vars%comp_num_eval, &
    pzheevx_vars%comp_num_evec, eval, orfac, evec%mat, 1, 1, evec%desca, work, &
    lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)

!    call pzheevr(pzheevx_vars%comp_evec, pzheevx_vars%range_, 'U', moire%natom, &
!    hamiltonian%mat, ia, ja, hamiltonian%desca, pzheevx_vars%vl, pzheevx_vars%vu, &
!    pzheevx_vars%il, pzheevx_vars%iu, pzheevx_vars%comp_num_eval,  &
!    pzheevx_vars%comp_num_evec, eval, evec%mat, 1, 1, evec%desca, work, lwork, &
!    rwork, lrwork, iwork, liwork, info)

    ! Find lwork 
    ! http://www.netlib.org/scalapack/explore-html/d8/d3b/pzheevx_8f_source.html
    ! ----------
    lwork  = int(abs(work(1)))+1
    anb = pjlaenv(hamiltonian%desca(CTXT_), 3, 'PZHETTRD', 'L', 0, 0, 0, 0)
    sqnpc = sqrt(dble(grid%nprow*grid%npcol))
    nps = max(numroc(moire%natom, 1, 0, 0, sqnpc),2)
    nhetrd_lwork = moire%natom + 2*(anb+1)*(4*nps+2)+(nps+1)*nps
    lwork = max(lwork,nhetrd_lwork)


    ! Find lrwork
    ! http://www.netlib.org/scalapack/explore-html/d8/d3b/pzheevx_8f_source.html
    ! ----------

    lrwork = int(abs(rwork(1)))+1
    nn = max(moire%natom, hamiltonian%desca(NB_), 2)
    neig = moire%natom
    np0 = numroc(nn, hamiltonian%desca(NB_), 0,0, grid%nprow)
    mq0 = numroc(max(neig, hamiltonian%desca(NB_),2), hamiltonian%desca(NB_), 0,0,grid%npcol)
    lrwork = max(4*moire%natom + max(5*nn, np0*mq0+2*pzheevx_vars%mb*pzheevx_vars%nb) + &
                 iceil(neig, grid%nprow*grid%npcol)*nn, lrwork)

    ! Find liwork
    ! http://www.netlib.org/scalapack/explore-html/d8/d3b/pzheevx_8f_source.html
    ! ------------
    liwork = int(abs(iwork(1)))+1
    liwork = max(liwork, 6*max(moire%natom, grid%npcol*grid%nprow+1,4))
    
    deallocate(work)
    deallocate(rwork)
    deallocate(iwork)

    allocate(work(lwork))
    allocate(rwork(lrwork))
    allocate(iwork(liwork))
  
    call pzheevx(pzheevx_vars%comp_evec, pzheevx_vars%range_, 'U', moire%natom,   &
    hamiltonian%mat, ia, ja, hamiltonian%desca, pzheevx_vars%vl, pzheevx_vars%vu, &
    pzheevx_vars%il, pzheevx_vars%iu, abstol, pzheevx_vars%comp_num_eval, &
    pzheevx_vars%comp_num_evec, eval ,orfac, evec%mat, 1, 1, evec%desca, work,    &
    lwork, rwork, lrwork, iwork, liwork, ifail, iclustr, gap, info)
    
!    call pzheevr(pzheevx_vars%comp_evec, pzheevx_vars%range_, 'U', moire%natom, &
!    hamiltonian%mat, ia, ja, hamiltonian%desca, pzheevx_vars%vl, pzheevx_vars%vu, &
!    pzheevx_vars%il, pzheevx_vars%iu, pzheevx_vars%comp_num_eval,  &
!    pzheevx_vars%comp_num_evec, eval, evec%mat, 1, 1, evec%desca, work, lwork, &
!    rwork, lrwork, iwork, liwork, info)


    if (info.gt.0) then
        if (mod(info,2).ne.0) then
            write(debug_str,'(A)') "One or more eigenvalues failed to converge"
            call debug_output(0)
        end if
        if (mod(info/2, 2).ne.0) then
            write(*,*) "Eigenvectors corresponding to the following indices could not be orthogonalized"
            write(*,*) iclustr
        end if
        if (mod(info/4, 2).ne.0) then
            write(*,*) "All eigenvectors in the specified interval could not be computed due to insufficient memory"
        end if
        if (mod(info/8,2).ne.0) then
            write(*,*) "One or more eigenvalues were not computed"
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
