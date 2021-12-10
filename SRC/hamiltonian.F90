subroutine create_hamiltonian(k_indx)

    use mpi
    use global_variables

    implicit none
    integer, intent(in) :: k_indx
    integer :: ia_first, ja_first, iastart, jastart, iaend, jaend , rsrc, csrc
    integer :: lroffset, lcoffset, ia, ja, lrindx, lcindx, ipos
    double complex :: hij

    rsrc = 0
    csrc = 0
    
    if (grid%myprow.ge.hamiltonian%desca(RSRC_)) then
        ia_first = (grid%myprow - hamiltonian%desca(RSRC_))*hamiltonian%desca(MB_)+1
    else
        ia_first = (grid%myprow + (grid%nprow - hamiltonian%desca(RSRC_)))* &
                    hamiltonian%desca(MB_) + 1
    endif

    if (grid%mypcol.ge.hamiltonian%desca(CSRC_)) then
        ja_first = (grid%mypcol - hamiltonian%desca(CSRC_))*hamiltonian%desca(NB_)+1
    else
        ja_first = (grid%mypcol + (grid%npcol - hamiltonian%desca(CSRC_)))* &
                    hamiltonian%desca(NB_) + 1
    endif

    do jastart = ja_first, hamiltonian%desca(N_), grid%npcol*hamiltonian%desca(NB_)
        do iastart = ia_first, hamiltonian%desca(M_), grid%nprow*hamiltonian%desca(MB_)
            iaend = min(hamiltonian%desca(M_), iastart+hamiltonian%desca(MB_)-1)
            jaend = min(hamiltonian%desca(N_), jastart+hamiltonian%desca(NB_)-1)

            ia = iastart
            ja = jastart

            call infog2l (ia, ja, hamiltonian%desca, grid%nprow, grid%npcol, &
                          grid%myprow, grid%mypcol, lroffset, lcoffset, rsrc, csrc)

            do ja=jastart,jaend
                do ia=iastart,iaend
                    lrindx = lroffset + (ia-iastart)
                    lcindx = lcoffset + (ja-jastart)
                    ipos = lrindx + (lcindx-1)*hamiltonian%desca(LLD_)

                    call compute_hij(ia,ja,k_indx,hij)

                    hamiltonian%mat(ipos) = hij

                end do
            end do

        end do
    end do

#ifdef __KPOOL    
    call mpi_barrier(mpi_local%comm, mpierr)

    write(debug_str,'(A,I0,A)') "Hamiltonian created for k-pt no.", k_indx ," on "

    call date_time_message_local(trim(debug_str))
#else
    call mpi_barrier(mpi_global%comm, mpierr)

    write(debug_str,'(A,I0,A)') "Hamiltonian created for k-pt no. ", k_indx, " on "
    call date_time_message(trim(debug_str))
#endif

end subroutine


subroutine compute_hij(i,j,k_indx,hij)

    use global_variables

    implicit none

    integer, intent(in) :: i,j, k_indx
    double complex, intent(out) :: hij
    double precision, dimension(3) :: ri, rj, ri_crys, rj_crys, rij, rij_crys
    integer :: l, m
    double precision :: kr, t
    double complex :: phase
    double complex, parameter :: IM = (0,1.0)
    double complex, parameter :: PI = 3.141592653589793

    ri = moire%real_pos(i,:)
    ri_crys = moire%crys(i,:)

    hij = cmplx(0,0)

    if (i.ne.j) then
        do l = -no_neigh, no_neigh
           do m = -no_neigh, no_neigh
              rj_crys(1) = moire%crys(j,1) + l
              rj_crys(2) = moire%crys(j,2) + m
              rj_crys(3) = moire%crys(j,3)
              rj = matmul(transpose(moire%lat),rj_crys)
              rij = ri - rj
              rij_crys = ri_crys - rj_crys
              kr = dot_product(k_file%points(k_indx,:),rij_crys)
              phase = exp(2*IM*PI*kr)
              call T_ij(i,j,rij,t)
              hij = hij + t*phase
           end do
        end do
    else
        if (i.le.int(moire%natom/2)) then
            hij = cmplx( E_field/2000, 0)
        else 
            hij = cmplx(-E_field/2000, 0)
        end if
    end if

    return

end subroutine


subroutine T_ij(i,j,rij,t)

    use global_variables
    implicit none
    integer, intent(in) :: i, j
    double precision, dimension(3), intent(in) :: rij
    double precision, intent(out) :: t
    double precision :: t_pi, t_sig

    call V_pi(i,j,rij,t_pi)
    call V_sig(i,j,rij,t_sig)

    t = t_pi + t_sig

    return

end subroutine 


subroutine V_pi(i,j,rij,t)
    use global_variables
    implicit none

    integer, intent(in) :: i, j
    double precision, dimension(3), intent(in) :: rij
    double precision, intent(out) :: t

    double precision, dimension(3) :: normalized_rij
    double precision, parameter :: Vpi0 = -2.7 ! eV
    double precision, parameter :: qpi = 3.1377732021175317 ! Angstrom
    double precision, parameter :: api = 1.42  ! Angstrom
    double precision :: norm, nri, nrj
    integer :: m
    double precision, dimension(3) :: t1, t2

    norm = norm2(rij)

    t = Vpi0 * exp(qpi*(1-(norm/api)))

    normalized_rij = rij/norm
    nri = dot_product(moire%normal(i,:),normalized_rij)
    nrj = dot_product(moire%normal(j,:),normalized_rij)

    do m = 1,3 
        t1(m) = moire%normal(i,m) - normalized_rij(m)*nri
        t2(m) = moire%normal(j,m) - normalized_rij(m)*nrj
    end do

    t = t*dot_product(t1,t2)

    return

end subroutine


subroutine V_sig(i,j,rij,t)

    use global_variables
    implicit none

    integer, intent(in) :: i,j 
    double precision, dimension(3), intent(in) :: rij
    double precision, intent(out) :: t

    double precision, dimension(3) :: normalized_rij
    double precision, parameter :: Vsig0 = 0.48 ! eV
    double precision, parameter :: qsig = 7.402493117671642 ! Angstrom
    double precision, parameter :: asig = 3.35  ! Angstrom

    double precision :: norm

    norm = norm2(rij)

    t = Vsig0*exp(qsig*(1-(norm/asig)))
    
    normalized_rij = rij/norm
    
    t = t*dot_product(moire%normal(i,:),normalized_rij) * &
          dot_product(moire%normal(j,:),normalized_rij)

    return

end subroutine
