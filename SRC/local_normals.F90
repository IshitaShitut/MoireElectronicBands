subroutine compute_local_normals()

    use global_variables
    use mpi

    implicit none
    integer :: i_start, i_end, counter, i, j, p, q, k, l
    double precision, dimension(3) :: ri, rj, n, v1, v2, nkl
    double precision :: dist
    double precision, dimension(3,4) :: min_list

    allocate(moire%normal(moire%natom,3))
    moire%normal = 0.0

    call distribution_length(moire%natom, mpi_global%rank, mpi_global%size_, &
                             i_start, i_end)

    do i = i_start, i_end
        
        ri = moire%real_pos(i,:)
        counter = 0
        do j = 1, moire%natom
            if (j.ne.i) then
                do p = -1,1
                    do q = -1,1
                        rj(1) = moire%crys(j,1)+p
                        rj(2) = moire%crys(j,2)+q
                        rj(3) = moire%crys(j,3)
                        rj = matmul(transpose(moire%lat),rj)
                        dist = sqrt((moire%real_pos(i,1)-rj(1))**2 + &
                                    (moire%real_pos(i,2)-rj(2))**2 + &
                                    (moire%real_pos(i,3)-rj(3))**2)
                        if (counter<3) then 
                            counter = counter+1
                            min_list(counter,1) = rj(1)
                            min_list(counter,2) = rj(2)
                            min_list(counter,3) = rj(3)
                            min_list(counter,4) = dist
                        else
                            call insert_in_list(min_list,rj,dist)
                        end if
                    end do
                end do 
            end if
        end do

        n = 0.0

        do k=1,3
            v1 = min_list(k,1:3) - ri
            v1 = v1/norm2(v1)
            do l = k+1,3
                v2 = min_list(l,1:3) - ri
                v2 = v2/norm2(v2)
                call cross_product(v1,v2,nkl)
                n = n+nkl
            end do
        end do

        n = n/norm2(n)
        
        moire%normal(i,:) = n
    end do

    call mpi_allreduce(MPI_IN_PLACE, moire%normal, moire%natom*3, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, mpi_global%comm, mpierr)


    call mpi_barrier(mpi_global%comm, mpierr)           
#ifdef __DEBUG
    write(debug_str,'(A)') '\r\n The local normals computed are:'
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)           
    do i=1,moire%natom
        write(debug_str,'(I8,3F16.6)') i, moire%normal(i,1), moire%normal(i,2), &
                                          moire%normal(i,3)
        call debug_output(0)
        call mpi_barrier(mpi_global%comm, mpierr)           
    end do
#endif
    
    call date_time_message("Local normals computed for the system computed on")
               
    return

end subroutine



subroutine insert_in_list(list, r, dist)
    implicit none
    real*8, dimension(3,4), intent(inout) :: list
    real*8, dimension(3), intent(in) :: r
    real*8, intent(in) :: dist
    integer :: i, j, k
    logical :: belong

    call sortcol(list)
    belong = .false.
    do i = 1,3
        if (dist < list(i,4)) then
            belong = .true.
            exit
        end if
    end do
    if (belong) then
        do j = 2, i, -1
            do k = 1,4
                list(j+1,k) = list(j,k)
            end do
        end do
        do j = 1,3
            list(i,j) = r(j)
        end do
        list(i,4) = dist
    end if
    return
end subroutine insert_in_list



subroutine sortcol(a)
    implicit none
    real*8, intent(inout) :: a(3, 4)
    integer :: i, j
    real*8 :: tmp(4)
    do i = 1, 3
       do j=i+1,3
           if (a(i,4).gt.a(j,4)) then
               tmp = a(i,:)
               a(i,:) = a(j,:)
               a(j,:) = tmp
           end if
       end do
    end do
end subroutine sortcol
