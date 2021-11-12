module fvdom
  use types
  implicit none

contains

  pure function absorb(temp)
    use types
    implicit none
    !$acc routine seq
    real(rk) :: absorb
    real(rks), intent(IN) :: temp
    !     H2O, CO2, CO according to:
    !
    !     Chmielewski, Maciej, and Marian Gieras. "Planck Mean Absorption
    !     Coefficients of H2O, CO2, CO and NO for radiation numerical
    !     modeling in combusting flows."  Journal of Power Technologies
    !     95  .2 (2015): 97-104
    !
    !     CH4 according to:
    !          
    !    https://www.sandia.gov/TNF/radiation.html
    !
    real(rk) :: tt
    real(rk), dimension(4) :: A_H2O, B_H2O, C_H2O
    integer(ik) :: i
    A_H2O=(/ 63.87D0, 2.629D0,  222.9D0,  0.9910D0 /)
    B_H2O=(/ 149.5D0, 493.6D0, -2110.0D0, 1700.0D0 /)
    C_H2O=(/ 174.3D0, 111.1D0,  1592.0D0, 619.20D0 /)

    tt = max(300.0_rk, min(2500.0_rk, real(temp, rk)))

    absorb = 0.0D0
    do i=1,4
       absorb=absorb + A_H2O(i) * &
            &exp(-((tt - B_H2O(i)) / C_H2O(i)) ** 2)
    end do

    return

  end function absorb

  subroutine calcia
    use parameters, only: n1, nsects, niterdo, fvtol
    use data, only: istart, iend, istep, jstart, jend, jstep, &
         &kstart, kend, kstep, nn, ghostleft, ghostright,&
         &temp, ia, iba, ntemp, u, sgn, s, fu, copy, left, right,&
         &is_wq
    use mpi
    use mpivars, only: MPIRK, MPI2RK, sbuf, rbuf, mpierr, mpirank
    use numerics, only: fourier
    implicit none
    integer(ik) :: ns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Explicit MPI method for the calculation of the radiative intensity !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk) :: err, tmperr, maxerr, gmaxerr ! maximum difference of radiative intensity between iterations
    real(rk) :: diff, tmp, val, maxga
    integer(ik) :: i, j, k, m, sd, mp ! loop counters
    real(rk), dimension(3) :: shat ! unit vector correponding to s (vector of direction cosines) 
    integer :: ierr, jerr, kerr, nserr
    integer :: nit, maxiter, maxloc, maxrank
    real(rk), dimension(1, 2) :: sbuf2, rbuf2   
    integer(ik) :: maxerr_rank


    ! copy temperature
    call copy(nn, u, fu)
    call fourier(nn, -1_ik, u, nfs=ntemp, nfe=ntemp)
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       temp(i, j, k) = u(i, j, k, ntemp)
    end do; end do ; end do
    !$acc update device(temp(1:nn(1), 1:nn(2), 1:nn(3)),  ia(1:nn(1), 1:nn(2), 0:nn(3)+1, 1:nsects))

    ! iteration loop      
    iterloop:do nit=1,niterdo

       !$acc update self(left(1:nn(1), 1:nn(2), 1:nsects), &
       !$acc& right(1:nn(1), 1:nn(2), 1:nsects))
       call ghost_nodes(left, right, ghostleft, ghostright)
       !$acc update device(ghostleft(1:nn(1), 1:nn(2), 1:nsects),&
       !$acc& ghostright(1:nn(1), 1:nn(2), 1:nsects))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TODO: define ghost cells !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       maxerr = 0.0_rk
       err = 0.0_rk
       !$acc parallel num_gangs(nsects) &
       !$acc& vector_length(n1) reduction(max:maxerr)
       do ns=1,nsects ; do j=1,nn(2) ; do i=1,nn(1)
          ia(i, j, 0_ik, ns) = ghostleft(i, j, ns)
          ia(i, j, nn(3) + 1, ns) = ghostright(i, j, ns)
       end do; end do ; end do 
       !acc loop independent collapse(4)
       do ns=1,nsects ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          iba(i, j, k, ns) = ia(i, j, k, ns)
       end do; end do ; end do ; end do
       !$acc end loop


       !$acc loop seq
       do sd=1,8 ! sweep the domain in 8 directions, one from each corner
          !$acc loop independent gang 
          do ns=1,nsects             
             ! sweep...
             if(is_wq(sd, ns)) then
                !$acc loop seq
                do k=kstart(sd),kend(sd),kstep(sd)
                   !$acc loop seq
                   do j=jstart(sd),jend(sd),jstep(sd)
                      !$acc loop independent vector
                      do i=istart(sd),iend(sd),istep(sd)
                         ! update the cell's radiative intensity
                         ! according to the step scheme
                         call cell_step_scheme(ns, i, j, k)
                      end do
                   end do
                end do
             end if
          end do
       end do
       !$acc end loop
       ! calculate maximum error and its position for each MPI process
       maxerr = 0.0
       err = 0.0
       ierr = 0
       jerr = 0
       kerr = 0
       !$acc loop independent reduction(max: maxerr)
       do ns=1,nsects
          !$acc loop independent 
          do k=1,nn(3)
             !$acc loop independent 
             do j=1,nn(2)
                !$acc loop independent
                do i=1,nn(1)
                   err = abs(ia(i, j, k, ns) - iba(i, j, k, ns))
                   if(err > maxerr) then
                      maxerr = err
                      !ierr = i
                      !jerr = j
                      !kerr = k
                      !nserr = ns
                   end if
                end do
             end do
          end do
       end do
       !$acc end loop

       !$acc loop seq collapse(3)
       do ns=1,nsects ; do j=1,nn(2) ; do i=1,nn(1)
          left(i, j, ns) = ia(i, j, 1_ik, ns)
          right(i, j, ns) = ia(i, j, nn(3), ns)
       end do; end do ; end do

       !$acc end parallel

!!$          !reduce maximum error across processes
!!$          sbuf(1) = maxerr
!!$          call mpi_allreduce(sbuf, rbuf, 1_i4b, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, mpierr)
!!$          maxerr = rbuf(1)
       ! global convergence condition


       sbuf(1) = maxerr
       call mpi_allreduce(sbuf, rbuf, 1_i4b, MPIRK, &
            &MPI_MAX, MPI_COMM_WORLD, mpierr);
       maxerr = rbuf(1)
       if(mpirank == 0) print '(i5,e15.5)', nit, maxerr
       if(maxerr < fvtol) then
          exit iterloop  
       end if

    end do iterloop

!!$
!!$    ! find the process in which the maximum error occurs
!!$    sbuf2(1,1) = maxerr
!!$    sbuf2(1,2) = mpirank
!!$    call mpi_allreduce(sbuf2, rbuf2, 1_i4b, MPI2RK, &
!!$         &MPI_MAXLOC, MPI_COMM_WORLD, mpierr); 
!!$    maxerr = rbuf2(1,1)
!!$    maxerr_rank = int(rbuf2(1,2), ik)
!!$
!!$    ! broadcast
!!$    sbuf(1) = ierr
!!$    sbuf(2) = jerr
!!$    sbuf(3) = kerr
!!$    sbuf(4) = nserr
!!$    call mpi_bcast(sbuf, 4_i4b, MPIRK, int(maxerr_rank, i4b), MPI_COMM_WORLD, mpierr)
!!$    ierr = sbuf(1)
!!$    jerr = sbuf(2)
!!$    kerr = sbuf(3)
!!$    nserr = sbuf(4)



    return

  end subroutine calcia

  subroutine ghost_nodes(left, right, ghostleft, ghostright)
    use mpi
    use mpivars, only: mpirank, MPIRKS, MPIROOT, &
         &mpisize, mpierr, sbuf, rbuf, lkstart, lksize
    use data, only: ia, sendbuf, recvbuf, nsects, nn
    implicit none
    real(rks), dimension(1:nn(1), 1:nn(2), 1:nsects), intent(IN) ::&
         & left, right
    real(rks), dimension(1:nn(1), 1:nn(2), 1:nsects), intent(OUT) ::&
         & ghostleft, ghostright
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(i4b) :: leftrank, rightrank, status(MPI_STATUS_SIZE), count, &
         &dir, rank
    integer(i4b), parameter :: L2R = 0, R2L = 1
    integer(ik) :: i, j, ns


    rank = mpirank

    
    if(rank == 0) then
       leftrank = mpisize - 1
    else
       leftrank = rank - 1
    end if
    
    if(rank == mpisize - 1) then
       rightrank = 0
    else
       rightrank = rank + 1
    end if
    
    count = int(nn(1) * nn(2) * nsects, i4b)
    

    !$omp parallel do
    do ns=1,nsects ; do j=1,nn(2) ; do i=1,nn(1)
       ghostleft(i, j, ns) = 0.0_rk
       ghostright(i, j, ns) = 0.0_rk
    end do; end do; end do
    !$omp end parallel do
    
    !$omp parallel do
    do ns=1,nsects ; do j=1,nn(2) ; do i=1,nn(1)
       sendbuf(i, j, ns) = right(i, j, ns)
    end do; end do; end do
    !$omp end parallel do
    call mpi_sendrecv(sendbuf, count, MPIRKS, rightrank, R2L,&
         & recvbuf, count, MPIRKS, leftrank, R2L, MPI_COMM_WORLD, &
         &status, mpierr)
    !$omp parallel do
    do ns=1,nsects ; do j=1,nn(2) ; do i=1,nn(1)
       ghostleft(i, j, ns) = recvbuf(i, j, ns)
    end do; end do; end do 

    !$omp end parallel do
    !$omp parallel do
    do ns=1,nsects ; do j=1,nn(2) ; do i=1,nn(1)
       sendbuf(i, j, ns) = left(i, j, ns)
    end do; end do; end do
    !$omp end parallel do
    call mpi_sendrecv(sendbuf, count, MPIRKS, leftrank, L2R,&
         & recvbuf, count, MPIRKS, rightrank, L2R, &
         &MPI_COMM_WORLD, status, mpierr)
    !$omp parallel do
    do ns=1,nsects ; do j=1,nn(2) ; do i=1,nn(1)
       ghostright(i, j, ns) = recvbuf(i, j, ns)
    end do; end do; end do
    !$omp end parallel do

    
    return

  contains

    subroutine copy_sendbuf(dir)
      use types
      implicit none
      integer(i4b) :: dir
      

      if(dir == R2L) then
         !$omp parallel do
         do ns=1,nsects ; do j=1,nn(2) ; do i=1,nn(1)
            sendbuf(i, j, ns) = ia(i, j, 1, ns)
         end do; end do; end do
         !$omp end parallel do
      else
         !$omp parallel do
         do ns=1,nsects ; do j=1,nn(2) ; do i=1,nn(1)
            sendbuf(i, j, ns) = ia(i, j, nn(3), ns)
         end do; end do; end do
         !$omp end parallel do
      end if

      return
      
    end subroutine copy_sendbuf

    subroutine copy_recvbuf(dir)
      use types
      implicit none
      integer(i4b) :: dir

      if(dir == R2L) then
         !$omp parallel do
         do ns=1,nsects ; do j=1,nn(2) ; do i=1,nn(1)
            ghostright(i, j, ns) = recvbuf(i, j, ns)
         end do; end do; end do
         !$omp end parallel do
      else
         !$omp parallel do
         do ns=1,nsects ; do j=1,nn(2) ; do i=1,nn(1)
            ghostleft(i, j, ns) = recvbuf(i, j, ns)
         end do; end do; end do 
         !$omp end parallel do
      end if

      return
      
    end subroutine copy_recvbuf

  end subroutine ghost_nodes
  
  
  subroutine calcqr
    use types
    use parameters, only: nsects
    use data, only: nn, qr, ga, ia, s, omeg
    use mpi
    use mpivars
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculates radiative heat flux and incindent radiation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: i, j, k, ns
    real(rk) :: maxga
    
    !$acc update self(ia(1:nn(1), 1:nn(2), 0:nn(3) + 1, 1:nsects))
    !radiative heat flux and incindent radiation at the interior of the domain
!    !$omp parallel do 
    do k=1,nn(3); do j=1,nn(2); do i=1,nn(1)
       ga(i, j, k) = 0.0
       qr(i, j, k, 1:3) = 0.0 
       do ns=1,nsects
          !qr = sum ia * s
          qr(i, j, k, 1:3) = qr(i, j, k, 1:3) + (ia(i, j, k, ns) * s(ns, 1:3)) ! add contribution of direction s
          !G = sum ia * omega
          ga(i, j, k) = ga(i, j, k) + (ia(i, j, k, ns) * omeg(ns)) ! same for the incindent radiation
       end do
    end do; end do; end do
 !   !$omp end parallel do

    sbuf(1) = maxval(ga)
    call mpi_allreduce(sbuf, rbuf, 1_i4b, MPIRK, &
         &MPI_MAX, MPI_COMM_WORLD, mpierr);
    maxga = rbuf(1)
    if(mpirank == 0) print '(a,e15.5)', 'max(G) = ', maxga
    return

  end subroutine calcqr
  
  subroutine cell_step_scheme(ns, i, j, k)
    use data, only: nn, ia, s, omeg, temp
    use parameters, only: PI, STEFB, LBOX
    implicit none
    !$acc routine seq
    integer(ik), intent(IN) :: ns ! number of direction
    integer(ik), intent(IN) :: i, j, k  ! number of cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Update the radiative intensity at the center of the cell according to the step scheme      !
    ! see: Modest, Michael F. Radiative heat transfer. Academic press, 2013., eq. (17.62) p. 568 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: dotprd, sumin, sumout 
    real(8) :: sp, nom, denom, fac
    real(8), dimension(6, 3) :: norm ! norm(nface, 1:3) unit vector normal to finite volume face
    real(8), dimension(6) :: surf ! surf(nface) finite volume face surface area
    real(8), dimension(6) :: faces_step ! radiative intensity on the finite volume face
    integer(8) :: nface
    real(8) :: vol

    ! set normal unit vectors for each face of the spatial control volume
    norm(1, 1:3) = (/  1.0,  0.0,  0.0 /)
    norm(2, 1:3) = (/ -1.0,  0.0,  0.0 /)
    norm(3, 1:3) = (/  0.0,  1.0,  0.0 /)
    norm(4, 1:3) = (/  0.0, -1.0,  0.0 /)
    norm(5, 1:3) = (/  0.0,  0.0,  1.0 /)
    norm(6, 1:3) = (/  0.0,  0.0, -1.0 /)

    ! initialize fluxes to zero
    faces_step(:) = 0.0

    ! initialize face areas to zero
    surf(:) = 0.0

    ! calculate intensity on faces and face surface area
    call faces_step_scheme(ns, i, j, k, faces_step, surf) 

    sumin = 0.0  ! sum of incoming intensities
    sumout = 0.0 ! sum of outgoing intensities
    !$acc loop seq
    do nface=1,6 ! loop over finite volume faces
       dotprd = dot_product(s(ns, :), norm(nface, :)) 
       if(dotprd < 0.0) then ! s is incoming
          sumin = sumin + faces_step(nface) * (-dotprd) * surf(nface)
       else ! s is outgoing
          sumout = sumout + dotprd * surf(nface)
       end if
    end do
    !$acc end loop

    vol = (LBOX / real(nn(1), rk)) ** 3
    
    sp = (STEFB / PI) * temp(i, j, k) ** 4

    fac = absorb(temp(i, j, k)) * vol * omeg(ns) ! auxiliary factor

    nom = fac * sp + sumin ! numerator of eq. (17.62) in (Modest, 2013)
    
    denom = fac + sumout ! denominator of eq. (17.62) 
    ia(i, j, k, ns) =  nom / denom ! update radiative intensity

    return

  end subroutine cell_step_scheme

  pure subroutine faces_step_scheme(ns, i, j, k, faces_step, surf)
    use data, only: nn, ia
    use parameters, only: PI, LBOX
    implicit none
    !$acc routine seq
    integer(ik), intent(IN) :: ns ! number of direction
    integer(ik), intent(IN) :: i, j, k  ! number of cell
    real(ik), dimension(6), intent(OUT) :: faces_step ! intensities at the faces of the cell
    real(ik), dimension(6), intent(OUT) :: surf ! surface area of the faces of the cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculates intensity at the faces of the cell according for the step scheme: !
    ! the intensity at the face is equal to the neighbouring intensity             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk) :: tmp
    integer(ik) :: ii, jj, kk
    faces_step(:) = 0.0
    surf(:) = 0.0

    tmp = (LBOX / real(nn(1), rk)) ** 2
    ! set surface area of cell faces
    surf(1:6) = tmp

    faces_step(1) = ia(per(i + 1), j, k, ns)
    faces_step(2) = ia(per(i - 1), j, k, ns)
    faces_step(3) = ia(i, per(j + 1), k, ns)
    faces_step(4) = ia(i, per(j - 1), k, ns)
    faces_step(5) = ia(i, j, k + 1, ns)
    faces_step(6) = ia(i, j, k - 1, ns)

    return

  end subroutine faces_step_scheme

  pure function per(i)
    use data, only: nn
    use types
    use parameters, only: PI
    implicit none
    !$acc routine seq
    integer(ik) :: per
    integer(ik), intent(in) :: i

    
    if(i == 0_ik) then
       per = nn(1)
    else if(i == nn(1) + 1) then
       per = 1
    else
       per = i
    end if
        
    return

  end function per

end module fvdom


!!$     if(mod(rank, 2) == 0) then
!!$       dir = R2L
!!$       call copy_sendbuf(dir)
!!$       call mpi_send(sendbuf, count, MPIRKS, right, dir, MPI_COMM_WORLD, mpierr)
!!$       call mpi_recv(recvbuf, count, MPIRKS, left, dir,  MPI_COMM_WORLD,&
!!$            & status, mpierr)
!!$       call copy_recvbuf(dir)
!!$       dir = L2R
!!$       call copy_sendbuf(dir)
!!$       call mpi_send(sendbuf, count, MPIRKS, left, dir, MPI_COMM_WORLD, mpierr)
!!$       call mpi_recv(recvbuf, count, MPIRKS, right, dir,  MPI_COMM_WORLD,&
!!$            & status, mpierr)
!!$       call copy_recvbuf(dir)
!!$    else
!!$       dir = R2L
!!$       call mpi_recv(recvbuf, count, MPIRKS, left, dir,  MPI_COMM_WORLD,&
!!$            & status, mpierr)
!!$       call copy_recvbuf(dir)
!!$       call copy_sendbuf(dir)
!!$       call mpi_send(sendbuf, count, MPIRKS, right, dir, MPI_COMM_WORLD, mpierr)
!!$       dir = L2R
!!$       call mpi_recv(recvbuf, count, MPIRKS, right, dir,  MPI_COMM_WORLD, &
!!$            &status, mpierr)
!!$       call copy_recvbuf(dir)
!!$       call copy_sendbuf(dir)
!!$       call mpi_send(sendbuf, count, MPIRKS, left, dir, MPI_COMM_WORLD, mpierr)
!!$    end if
