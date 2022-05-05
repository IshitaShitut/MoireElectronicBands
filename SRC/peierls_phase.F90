subroutine create_peierls_phase(i,j,p_phase)

    use mpi
    use global_variables

    ! ------------------------------------------------------------------------
    ! Subroutine generates Peierls phase for the case of constant
    ! magnetic field applied to a 2D material in tight binding approximation
    !
    ! The subroutine makes use of 3 gauges for vector potential:
    !      Gauge 1 :
    !      Gauge 2 :
    !      Gauge 3 :
    !
    ! All of the quantities are converted in Hartree atomic units where
    ! e = hbar = 1
    ! Peirels phase is a dimensionless quantity
    ! ------------------------------------------------------------------------

    implicit none
    integer, intent(in) :: i, j
    double complex, intent(out) :: p_phase
    double precision :: theta, autotesla, bohrtoang

    autotesla = 23505.18d0  !Need to change
    bohrtoang = 0.529177d0

    select case (B_gauge)

        case(1)
            theta = (moire%real_pos(j,1)+moire%real_pos(i,1))*(moire%real_pos(j,2)-moire%real_pos(j,1))
            theta = (B_field/(2.d0*autotesla))*(theta/bohrtoang**2)

        case(2)
            theta = (moire%real_pos(j,1)-moire%real_pos(i,1))*(moire%real_pos(j,2)+moire%real_pos(j,1))
            theta = (B_field/(2.d0*autotesla))*(theta/bohrtoang**2)

        case(3)
            theta = moire%real_pos(i,1)*moire%real_pos(j,2)-moire%real_pos(j,1)*moire%real_pos(i,2)
            theta = (B_field/(2.d0*autotesla))*(theta/bohrtoang**2)

    end select

    p_phase = exp(-cmplx(0.0d0,theta))

end subroutine

