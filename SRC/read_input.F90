subroutine read_input()
    
    ! Subroutine to read the input file
    ! ----------------------------------

    use global_variables
    implicit none
    
    character(len=char_len) :: input_file
    integer :: error
    character(len=char_len) :: data_
    character(len=char_len) :: format_ 
    integer :: pos, i
    character(char_len), dimension(2) :: args1_
    character(char_len), allocatable, dimension(:) :: args2_
    integer, parameter :: inp_unit = 20
    character(1), external :: capital

    call default_variables()

    call get_command_argument(1,input_file)
    if (len_trim(input_file)==0) then
        err_msg = "\r\n#!#!#!#!#!#!#!#!\r\nError! No input file"
        call error_message()
        call close_mpi()
        call exit
    end if

    call mpi_barrier(mpi_global%comm, mpierr)

    open(unit=inp_unit, file=trim(adjustl(input_file)), iostat=error)

    if (error.ne.0) then
        write(err_msg,'(2A)') 'Error reading input file ', trim(adjustl(input_file))
        call error_message()
        call exit
    end if

    debug_str = ' '
    call debug_output(0)
    write(debug_str,'(2A)') "Reading input file : ", trim(adjustl(input_file))
    call debug_output(0)


    ! -------------------------------------------------------------------------------- !
    !                                                                                  !
    ! Read all lines in the input file.                                                !
    ! If a "#" or "!" is encountered, the line is treated as a comment                 ! 
    ! Ignores empty lines                                                              !
    !                                                                                  !
    ! Splits the lines taking the delimiter as ':'                                     !
    ! The left hand side of the delimiter is used to identify the variable being read  !
    ! The right hand side is stored directly to the necessary variables                !
    !                                                                                  !
    ! -------------------------------------------------------------------------------- !


    do 
100     read(inp_unit,'(A)',iostat = error) data_
        data_ = trim(adjustl(data_))
        if (error.eq.0) then
            if (data_(1:1) == '#'.or.data_(1:1) == '!'.or.data_(1:1) == ' ' ) goto 100 
            
            pos = index(data_, '!')
            if (pos .ne. 0) then
                data_ = data_(1:pos-1)
            end if
            
            pos = index(data_, '#')
            if (pos .ne. 0) then
                data_ = data_(1:pos-1)
            end if

            pos = index(data_,':')
            args1_(1) = data_(1:pos-1)
            args1_(2) = data_(pos+1:)

            do i=1,len_trim(args1_(1))
                args1_(1)(i:i) = capital(args1_(1)(i:i))
            end do            
            
            select case(trim(args1_(1)))
                case ("LAMMPS FILE LOCATION", "LAMMPS_FILE_LOCATION")
                    lammps_file%location = trim(adjustl(args1_(2)))
                case ("LAMMPS_FILE_NAME", "LAMMPS FILE NAME")
                    lammps_file%name_ = trim(adjustl(args1_(2)))
                case ("NATOM")
                    read(args1_(2), *) moire%natom
                case ("ATOM_STYLE", "ATOM STYLE")
                    lammps_file%atom_style = trim(adjustl(args1_(2)))
                case ("ATOM_TYPES", "ATOM TYPES")
                    read(args1_(2),*) lammps_file%at_types
                case ("K_FILE_LOCATION", "K FILE LOCATION")
                    k_file%location = trim(adjustl(args1_(2)))
                case ("K_FILE_NAME","K FILE NAME")
                    k_file%name_ = trim(adjustl(args1_(2)))
                case ("NKPT")
                    read(args1_(2),*) k_file%npt
#ifdef __KPOOL
                case ("NUM_KPOOLS","NUM KPOOLS")
                    read(args1_(2), *) num_pools
#endif  
                case ("COMPUTE_EIGVECS","COMPUTE EIGVECS")
                    args1_(2) = trim(adjustl(args1_(2)))
                    if (args1_(2) == 'yes' .or. args1_(2) == 'true' .or. &
                        args1_(2) == 'Yes' .or. args1_(2) == 'True' .or. &
                        args1_(2) == 'YES' .or. args1_(2) == 'TRUE') then
                        pzheevx_vars%comp_evec = 'V'
                        evec_comp = .true.
                    elseif (args1_(2) == 'no' .or. args1_(2) == 'false' .or. &
                            args1_(2) == 'No' .or. args1_(2) == 'False' .or. &
                            args1_(2) == 'NO' .or. args1_(2) == 'FALSE') then
                        pzheevx_vars%comp_evec = 'N'
                        evec_comp = .false.
                    else
                        write(err_msg, '(3A)') "Could not interpret command ", &
                                               args1_(2), &
                                               "\r\nEigenvectors will not be computed."
                        call error_message()
                        pzheevx_vars%comp_evec = 'N'
                        evec_comp = .false.
                    end if
                case ("GROUP VELOCITY", "GROUP_VELOCITY")
                    args1_(2) = trim(adjustl(args1_(2)))
                    if (args1_(2) == 'yes' .or. args1_(2) == 'true' .or. &
                        args1_(2) == 'Yes' .or. args1_(2) == 'True' .or. &
                        args1_(2) == 'YES' .or. args1_(2) == 'TRUE') then
                        comp_vel = .true.
                    elseif (args1_(2) == 'no' .or. args1_(2) == 'false' .or. &
                            args1_(2) == 'No' .or. args1_(2) == 'False' .or. &
                            args1_(2) == 'NO' .or. args1_(2) == 'FALSE') then
                        comp_vel = .false.
                    else
                        write(err_msg, '(3A)') "Could not interpret command ", &
                                               args1_(2), &
                                               "\r\nVelocity will not be computed."
                        call error_message()
                        comp_vel = .false.
                    end if
                case ("RANGE")
                    if (trim(adjustl(args1_(2))) == 'A') then
                        pzheevx_vars%range_ = 'A'
                    elseif (trim(adjustl(args1_(2))) == 'V') then 
                        pzheevx_vars%range_ = 'V'
                    elseif (trim(adjustl(args1_(2))) == 'I') then   
                        pzheevx_vars%range_ = 'I'
                    else 
                        write(err_msg, '(3A)') 'Could not interpret command', &
                                                trim(adjustl(args1_(2))), &
                                                '\r\n Setting range to A'
                    end if
                case ("MAX_EIGVAL", "MAX EIGVAL")
                    read(args1_(2), *) pzheevx_vars%vu
                case ("MIN_EIGVAL", "MIN EIGVAL")
                    read(args1_(2), *) pzheevx_vars%vl
                case ("MAX_INDEX", "MAX INDEX")
                    read(args1_(2), *) pzheevx_vars%iu
                case ("MIN_INDEX", "MIN INDEX")
                    read(args1_(2), *) pzheevx_vars%il
                case ("ABSTOL")
                    read(args1_(2), *) pzheevx_vars%abstol
                case ("ORFAC")
                    read(args1_(2), *) pzheevx_vars%orfac 
                case ("COMPUTE LOCAL NORMALS", "LOCAL_NORMALS", "LOCAL NORMALS", "COMPUTE_LOCAL_NORMALS")
                    args1_(2) = trim(adjustl(args1_(2)))
                    if (args1_(2) == 'yes' .or. args1_(2) == 'true' .or. &
                        args1_(2) == 'Yes' .or. args1_(2) == 'True' .or. &
                        args1_(2) == 'YES' .or. args1_(2) == 'TRUE') then
                        moire%compute_normal = .true.
                    elseif (args1_(2) == 'no' .or. args1_(2) == 'false' .or. &
                            args1_(2) == 'No' .or. args1_(2) == 'False' .or. &
                            args1_(2) == 'NO' .or. args1_(2) == 'FALSE') then
                        moire%compute_normal = .false.
                    else
                        write(err_msg, '(3A)') "Could not interpret command ", &
                                               args1_(2), &
                                               "\r\nVelocity will not be computed."
                        call error_message()
                        moire%compute_normal = .false.
                    end if
                case ("NUM NEIGHBOURS", "NUM_NEIGHBOURS")
                    read(args1_(2), *) no_neigh
                case ("E FIELD Z", "E_FIELD_Z")
                    read(args1_(2), *) E_field
                case ("B FIELD Z", "B_FIELD_Z")
                    read(args1_(2), *) B_field
                case ("B GAUGE", "B_GAUGE")
                    read(args1_(2), *) B_gauge
                case ("ONSITE_ENERGY", "ONSITE ENERGY")
                    if (lammps_file%at_types == 0) then
                        write(err_msg, '(2A)') "Atom type cannot be 0 ", &
                            "Atom type should be specified before Onsite Energy"
                        call error_message()
                        call exit
                    endif
                    allocate(moire%onsite_en(lammps_file%at_types))
                    if (lammps_file%at_types > 1) then
                        allocate(args2_(lammps_file%at_types))
                        args2_(1) = args1_(2)
                        do i = 1, (lammps_file%at_types-1)
                            pos = index(trim(args2_(i)),',')
                            data_ = args2_(i)
                            args2_(i) = data_(1:pos-1)
                            args2_(i+1) = data_(pos+1:)
                            read(args2_(i), *) moire%onsite_en(i)
                            read(args2_(i+1),*) moire%onsite_en(i+1)
                        end do
                    else
                        read(args1_(2), *) moire%onsite_en(:)
                    endif
                case ("MB")
                    read(args1_(2),*) pzheevx_vars%mb
                case ("NB")
                    read(args1_(2),*) pzheevx_vars%nb
                case ("OUTPUT FILE NAME", "OUTPUT_FILE_NAME")
                    write(output_file_name,'(A)') trim(adjustl(args1_(2)))
                case ("OUTPUT FILE LOCATION", "OUTPUT_FILE_LOCATION")
                    write(output_file_location,'(A)') trim(adjustl(args1_(2)))
                case default
                    write(err_msg, '(3A)') "\r\n%%%%%\r\nCommand ", trim(args1_(1)), " not recognised"
                    call error_message()
                    err_msg = "Ignoring this command\r\n"
                    call error_message()  
            end select

        elseif (error.gt.0) then
            err_msg = 'Error encountered while reading file'
            call error_message()
            call exit
        else
#ifdef __DEBUG
            debug_str = "\r\nFile read successfully"
            call debug_output(0)
#endif
            exit
        end if
    end do

    close(unit=inp_unit)

    call sanitize_input()

    debug_str = '\r\n========================================================='
    call debug_output(0)
    debug_str = '\r\nCalculations started with the following parameters: '
    call debug_output(0)
    write(debug_str, '( 2A )') "\r\nLAMMPS file location:  ", trim(lammps_file%location)
    call debug_output(0)
    write(debug_str, '(2A)') "LAMMPS file name:  ", trim(lammps_file%name_)
    call debug_output(0)
    write(debug_str, '(A, I0)') "Number of atoms:  ", moire%natom
    call debug_output(0)
    write(debug_str, '(A,I0)') "Atom types:  ", lammps_file%at_types
    call debug_output(0)
    write(debug_str, '(2A)') "LAMMPS atom style:  ", trim(lammps_file%atom_style)
    call debug_output(0)
    write(debug_str, '(2A)') "K file location:  ", trim(k_file%location)
    call debug_output(0)
    write(debug_str, '(2A)') "K file name:  ", trim(k_file%name_)
    call debug_output(0)
    write(debug_str, '(A,I0)') "Number of K-points:  ", k_file%npt
    call debug_output(0)
#ifdef __KPOOL
    write(debug_str, '(A,I0)') "Number of simultaneous K-point diagonalizations:  ", num_pools
    call debug_output(0)
#endif
    write(debug_str, '(A,L)') "Compute Eigenvectors: ", evec_comp
    call debug_output(0)   
    write(debug_str, '(A,L)') "Compute Velocities: ", comp_vel
    call debug_output(0) 
    write(debug_str, '(2A)') "Range of the calculation: ", trim(pzheevx_vars%range_)
    call debug_output(0)
    write(debug_str, '(A,F0.6)') "Minimum eigenvalue : ", pzheevx_vars%vl
    call debug_output(0)
    write(debug_str, '(A,F0.6)') "Maximum eigenvalue : ", pzheevx_vars%vu
    call debug_output(0)
    write(debug_str, '(A,I0)') "Minimum eigenvalue index : ", pzheevx_vars%il
    call debug_output(0)
    write(debug_str, '(A,I0)') "Maximum eigenvalue index : ", pzheevx_vars%iu
    call debug_output(0)
    write(debug_str, '(2(A,I0),A)') "Block size for scalapack diagonalization : (", &
                                     pzheevx_vars%mb,',',pzheevx_vars%nb,')'
    call debug_output(0)
    write(debug_str, '(A,E12.4)') "Abstol : ", pzheevx_vars%abstol
    call debug_output(0)
    write(debug_str, '(A,E12.4)') "Orfac : ", pzheevx_vars%orfac
    call debug_output(0)
    write(debug_str,'(A,L)') "Compute Local Normals : ", moire%compute_normal
    call debug_output(0)
    write(debug_str, '(A,I0)') "Number of neighbour cells in each direction to scan : ", &
                                 no_neigh
    call debug_output(0)
    write(debug_str, '(A,F0.6,A)') "Electric Field applied in the Z direction : ", &
                                    E_field, " meV"
    call debug_output(0)
    write(debug_str, '(A,F0.6,A)') "Magnetic Field applied in the Z direction : ", &
                                    B_field, " Tesla"
    call debug_output(0)
    if (B_gauge .eq. 1) then
        write(debug_str, '(A)') "Magnetic Vector Potential gauge : Landau Gauge"
    else if (B_gauge .eq. 2) then
        write(debug_str, '(A)') "Magnetic Vector Potential gauge : Landau Gauge"
    else if (B_gauge .eq. 3) then
        write(debug_str, '(A)') "Magnetic Vector Potential gauge : Symmetric Gauge"
    endif
    call debug_output(0)
    write(format_,'(A,I0,A)') '(A,',lammps_file%at_types,'(F0.6,2X),A)'
    write(debug_str, trim(adjustl(format_))) "Onsite energy: ", moire%onsite_en, " meV"
    call debug_output(0)
    write(debug_str, '(2A)') "Output file name : ", trim(output_file_name)
    call debug_output(0)
    write(debug_str, '(2A)') "Output file location : ", trim(output_file_location)
    call debug_output(0)
    debug_str = '\r\n========================================================='
    call debug_output(0)

    return

end subroutine






function capital (in_char)
    implicit none
    character(len=1), intent(in) :: in_char
    character(len=1) :: capital
    character(len=26), parameter :: lower='abcdefghijklmnopqrstuvwxyz', &
                                    upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer :: i

    do i=1,26
        if (in_char == lower(i:i)) then
            capital = upper(i:i)
            return
        end if
    end do
    capital = in_char
    return
end function capital



subroutine default_variables()
    use global_variables
    implicit none

    lammps_file%location = './'
    lammps_file%name_ = 'lammps.dat'
    lammps_file%atom_style = 'A'
    lammps_file%at_types = 0
    moire%natom = 1
    k_file%location = './'
    k_file%name_ = 'k_file.dat'
    k_file%npt = 1
    k_file%start = 1
    k_file%finish = 1
    evec_comp = .false.
    pzheevx_vars%comp_evec = 'N'
    comp_vel = .false.
    pzheevx_vars%range_ = 'A'
    pzheevx_vars%il = 1
    pzheevx_vars%iu = 1
    pzheevx_vars%vl = 0.0
    pzheevx_vars%vu = 0.0
    pzheevx_vars%mb = -1
    pzheevx_vars%nb = -1
    pzheevx_vars%abstol= 1E-15
    pzheevx_vars%orfac = 1E-15
    moire%compute_normal = .false.
    no_neigh = 1
    E_field = 0.0
    B_field = 0.0
    B_gauge = 3
    moire%onsite_en = 0.0
    output_file_name = 'results'
    output_file_location = './'

    return

end subroutine


subroutine sanitize_input()

    use global_variables
    implicit none
    integer :: len_
    character(char_len) :: file_name, full_file_details
    logical :: file_exists
    integer :: counter

    len_ = len_trim(lammps_file%location)

    if (lammps_file%location(len_:len_).ne.'/') then
        write(lammps_file%location,'(2A)') trim(adjustl(lammps_file%location)),'/'
    end if

    len_ = len_trim(k_file%location)

    if (k_file%location(len_:len_).ne.'/') then
        write(k_file%location,'(2A)') trim(adjustl(k_file%location)),'/'
    end if
   
    len_ = len_trim(output_file_location)

    if (output_file_location(len_:len_).ne.'/') then
        write(output_file_location,'(2A)') trim(adjustl(output_file_location)),'/'
    end if


    if (pzheevx_vars%mb == -1) then
        pzheevx_vars%mb = min(32,moire%natom)
    end if

    if (pzheevx_vars%nb == -1) then
        pzheevx_vars%nb = min(32,moire%natom)
    end if

    select case (pzheevx_vars%range_)
        case("V")
            pzheevx_vars%il = 1
            pzheevx_vars%iu = moire%natom
        case("I")
            pzheevx_vars%vl = 0.0
            pzheevx_vars%vu = 0.0
        case("A")
            pzheevx_vars%il = 1
            pzheevx_vars%iu = moire%natom
            pzheevx_vars%vl = 0.0
            pzheevx_vars%vu = 0.0
    end select

    if ((pzheevx_vars%range_ == 'I').and.(pzheevx_vars%il.gt.pzheevx_vars%iu)) then
        err_msg = "\r\n Invalid input for eigenvalue indices to be computed"
        call error_message()
        call exit
    end if

    if ((pzheevx_vars%range_ == 'V').and.(pzheevx_vars%vl .gt. pzheevx_vars%vu)) then
        err_msg = "\r\n Invalid input for the interval in which eigenvalue is to be computed"
        call error_message()
        call exit
    end if
  
    if (comp_vel) then
        pzheevx_vars%comp_evec='V'
    end if


    write(full_file_details, '(3A)') trim(adjustl(output_file_location)), &
                                     trim(adjustl(output_file_name)),'.hdf5'

    inquire(file=full_file_details, exist=file_exists)

    if (file_exists) then
        write(err_msg,'(5A)') "\r\n File ", trim(adjustl(output_file_name)),'.hdf5', &
                              " already exists at ", trim(adjustl(output_file_location))
        call error_message()

        counter = 0

        do while (file_exists) 
            counter = counter + 1
            write(file_name,'(2A,I0)') trim(adjustl(output_file_name)),'_v',counter
            write(full_file_details, '(3A)') trim(adjustl(output_file_location)), &
                                             trim(adjustl(file_name)),'.hdf5'
            inquire(file=full_file_details, exist=file_exists)
        end do

        write(output_file_name,'(2A)') trim(adjustl(file_name)),'.hdf5'
        write(err_msg,'(2A)') "\r\n Output file being renamed as ", &
                              trim(adjustl(output_file_name))
        call error_message()
    else
        write(output_file_name,'(2A)') trim(adjustl(output_file_name)),'.hdf5'
    end if

    return

end subroutine
