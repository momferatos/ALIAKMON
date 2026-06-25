module fvdom
  use types
  use parameters
  implicit none

contains

  pure function absorb(temp, y, p, mode)
    use types
    implicit none
    !$acc routine seq
    real(rk) :: absorb
    real(rk), intent(IN) :: temp, y, p
    integer(ik), intent(in), optional :: mode
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
    real(4) :: tt, yy
    real(4), dimension(4) :: A_H2O, B_H2O, C_H2O
    integer(ik) :: i, imode

    imode=0
    if(present(mode)) imode=mode

    A_H2O=(/ 63.87D0, 2.629D0,  222.9D0,  0.9910D0 /)
    B_H2O=(/ 149.5D0, 493.6D0, -2110.0D0, 1700.0D0 /)
    C_H2O=(/ 174.3D0, 111.1D0,  1592.0D0, 619.20D0 /)

    tt = max(TEMPMIN, min(TEMPMAX, real(temp, rk)))

    if(imode==0) then
       tt=1.0D3/tt
       absorb=-0.23093D0+tt*(-1.12390D0+tt*(9.41530D0+tt*&
            &(-2.99880D0+tt*(0.51382D0+tt*(-1.86840D-5)))))
    else
       absorb = 0.0e0
       do i=1,4
          absorb=absorb + A_H2O(i) * &
               &exp(-((tt - B_H2O(i)) / C_H2O(i)) ** 2)
       end do
    end if

    absorb = (PATM + p) * y * absorb

    return

  end function absorb

  subroutine calcia(nn, temp)
    use parameters, only: n1, nsects, niterdo, fvtol, LBOX, STEFB, PI
    use data, only: istart, iend, istep, jstart, jend, jstep, &
         &kstart, kend, kstep, ghostleft, ghostright,&
         &ia, iba, ntemp, sgn, s, copy, left, right,&
         &is_wq, omeg,copy,zero,ga,qr,press,dotprds
    use mpi
    use mpivars, only: MPIRK, MPI2RK, mpierr, mpirank
    implicit none
    integer(ik), dimension(1:4) :: nn
    real(rks), dimension(1:nn(1),1:nn(2),1:nn(3)), intent(IN) :: temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    integer(ik) :: ngangs, vlength
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: dotprd, sumin, sumout 
    real(8) :: sp, nom, denom, fac
    real(8) :: surf ! surf(nface) finite volume face surface area
    real(8), dimension(6) :: faces_step ! radiative intensity on the finite volume face
    integer(8) :: nface
    real(8) :: vol, y, T, p

    !$acc update device(temp(1:nn(1), 1:nn(2), 1:nn(3)))
    !$acc update device(press(1:nn(1), 1:nn(2), 1:nn(3)))


    ! iteration loop      
    iterloop:do nit=1,niterdo

       !$acc update self(left(1:nn(1), 1:nn(2), 1:nsects), &
       !$acc& right(1:nn(1), 1:nn(2), 1:nsects))
       call ghost_nodes(left, right, ghostleft, ghostright)
       !$acc update device(ghostleft(1:nn(1), 1:nn(2), 1:nsects),&
       !$acc& ghostright(1:nn(1), 1:nn(2), 1:nsects))

       maxerr = 0.0
       vlength=nsects
       ngangs=nsects

       !$acc data present(nn, ia, iba, is_wq, dotprds, temp, press)

       !$acc parallel copy(maxerr) &
       !$acc& private(dotprd, sumin, sumout , &
       !$acc& sp, nom, denom, fac, surf, faces_step, nface, &
       !$acc& vol, y, T, p, tmp)

#ifndef _OPENACC
       !$omp parallel do 
#endif
       !$acc loop independent collapse(3)
       do j=1,nn(2) ; do i=1,nn(1) ; do ns=1,nsects
          ia(ns, i, j, 0_ik) = ghostleft(i, j, ns)
          ia(ns, i, j, nn(3) + 1) = ghostright(i, j, ns)
       end do; end do ; end do
       !$acc end loop

#ifndef _OPENACC
       !$omp parallel do 
#endif
       !$acc loop independent collapse(4)
       do k=0,nn(3)+1 ; do j=0,nn(2)+1 ; do i=0,nn(1)+1 ; do ns=1,nsects
          iba(ns, i, j, k) = ia(ns, i, j, k)
       end do; end do ; end do ; end do
       !$acc end loop


       ! sweep the domain in 8 directions, one from each corneρ
       !$acc loop seq
       sweeploop: do sd=1,8

#ifndef _OPENACC
          !$omp parallel do 
#endif
          !$acc loop independent collapse(3)
          do k=1,nn(3) ; do j=1,nn(2) ; do ns=1,nsects
             ia(ns, 0, j, k) =  ia(ns, nn(1), j, k)  
             ia(ns, nn(1)+1, j, k) = ia(ns, 1, j, k)
          end do; end do ; end do
          !$acc end loop

#ifndef _OPENACC
          !$omp parallel do 
#endif
          !$acc loop independent collapse(3)
          do k=1,nn(3) ; do i=1,nn(1) ; do ns=1,nsects
             ia(ns, i, 0, k) = ia(ns, i, nn(2), k) 
             ia(ns, i, nn(2)+1, k) = ia(ns, i, 1, k)
          end do; end do ; end do
          !$acc end loop

#ifndef _OPENACC
          ! sweep...
          !$omp parallel do &
          !$omp& private(dotprd, sumin, sumout , &
          !$omp& sp, nom, denom, fac, surf, faces_step, nface, &
          !$omp& vol, y, T, p, tmp)
#endif
          ! sweep...
          !$acc loop independent collapse(4) &
          !$acc& private(dotprd, sumin, sumout , &
          !$acc& sp, nom, denom, fac, surf, faces_step, nface, &
          !$acc& vol, y, T, p, tmp)
          do k=kstart(sd),kend(sd),kstep(sd)
             do j=jstart(sd),jend(sd),jstep(sd)
                do i=istart(sd),iend(sd),istep(sd)
                   do ns=1,nsects
                      if(is_wq(ns,sd))  then
                         ! update the cell's radiative intensity
                         ! according to the step scheme

                         ! set surface area of cell faces
                         surf = (LBOX / real(nn(1), rk)) ** 2

                         faces_step(1) = ia(ns,i + 1, j, k)
                         faces_step(2) = ia(ns,i - 1, j, k)
                         faces_step(3) = ia(ns,i, j + 1, k)
                         faces_step(4) = ia(ns,i, j - 1, k)
                         faces_step(5) = ia(ns,i, j, k + 1)
                         faces_step(6) = ia(ns,i, j, k - 1)

                         sumin = 0.0  ! sum of incoming intensities
                         sumout = 0.0 ! sum of outgoing intensities
                         !$acc loop seq reduction(+:sumin,sumout)
                         do nface=1,6 ! loop over finite volume faces
                            dotprd = dotprds(nface,ns)
                            if(dotprd < 0.0) then ! s is incoming
                               sumin = sumin + faces_step(nface) * &
                                    &(-dotprd) * surf
                            else ! s is outgoing
                               sumout = sumout + dotprd * surf
                            end if
                         end do
                         !$acc end loop

                         vol = (LBOX / real(nn(1), rk)) ** 3

                         sp = (STEFB / PI) * temp(i, j, k) ** 4

                         T=temp(i,j,k)
                         y = (T - TEMPMIN) / (TEMPMAX - TEMPMIN)
                         y = max(0.0_rk, min(1.0_rk, real(y, rk)))
                         p = press(i,j,k)
                         fac = absorb(T, y, p, 0_ik) * &
                              &vol * omeg(ns) ! auxiliary factor

                         ! numerator of eq. (17.62) in (Modest, 2013)
                         nom = fac * sp + sumin 

                         denom = fac + sumout ! denominator of eq. (17.62)

                         ! update radiative intensity
                         ia(ns, i, j, k) =  nom / denom 
                         !call cell_step_scheme(ns, i, j, k)
                      end if
                   end do
                end do
             end do
          end do
          !$acc end loop
       end do sweeploop


       ! calculate maximum error
#ifndef _OPENACC
       !$omp parallel do reduction(max:maxerr)
#endif
       !$acc loop independent reduction(max:maxerr) collapse(4)
       do k=1,nn(3)
          do j=1,nn(2)
             do i=1,nn(1)
                do ns=1,nsects
                   ! relative change; skip cells with zero intensity to
                   ! avoid a divide-by-zero poisoning the max reduction
                   if(ia(ns, i, j, k) /= 0.0_rk) then
                      maxerr = max(maxerr, real(abs((ia(ns, i, j, k) &
                           &- iba(ns, i, j, k)) / ia(ns, i, j, k)), rk))
                   end if
                end do
             end do
          end do
       end do
       !$acc end loop

       !$acc loop independent collapse(3)
       do j=1,nn(2) ; do i=1,nn(1) ; do ns=1,nsects
          left(i, j, ns) = ia(ns, i, j, 1_ik)
          right(i, j, ns) = ia(ns, i, j, nn(3))
       end do; end do ; end do
       !$acc end loop

       !$acc end parallel

       !$acc end data

       ! reduce maximum error across processes
       call mpi_allreduce(MPI_IN_PLACE,maxerr,1_i4b,MPIRK,&
            &MPI_MAX,MPI_COMM_WORLD,mpierr)

       if(mpirank == 0) print '(i5,e15.5)', nit, maxerr

       ! comnvergence condition
       if(maxerr < FVTOL) then
          exit iterloop  
       end if

    end do iterloop

    call calcqr


    return

  end subroutine calcia

  subroutine ghost_nodes(left, right, ghostleft, ghostright)
    use mpi
    use mpivars, only: mpirank, MPIRKS, MPIROOT, &
         &mpisize, mpierr, lkstart, lksize
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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
         do j=1,nn(2) ; do i=1,nn(1) ; do ns=1,nsects
            sendbuf(i, j, ns) = ia(ns, i, j, 1)
         end do; end do; end do
         !$omp end parallel do
      else
         !$omp parallel do
         do j=1,nn(2) ; do i=1,nn(1) ; do ns=1,nsects
            sendbuf(i, j, ns) = ia(ns, i, j, nn(3))
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
    real(rk) :: maxga, tqr1, tqr2, tqr3, tga

    !$acc data present(ia, ga, qr)

    !$acc parallel
#ifndef _OPENACC
    !$omp parallel do collapse(3)
#endif
    !$acc loop independent collapse(3)
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       ga(i, j, k) = 0.0
       qr(i, j, k, 1:3) = 0.0
    end do;  end do ;  end do
    !$acc end loop

    !$acc loop independent private(tqr1,tqr2,tqr3,tga)
    do k=1,nn(3); do j=1,nn(2); do i=1,nn(1)
       tqr1=0.0
       tqr2=0.0
       tqr3=0.0
       tga=0.0
#ifndef _OPENACC
       !$omp parallel do reduction(+:tqr1,tqr2,tqr3,tga)
#endif
       !$acc loop independent reduction(+:tqr1,tqr2,tqr3,tga)
       do ns=1,nsects
          !qr = sum ia * s
          ! add contribution of direction s
          tqr1 = tqr1 + (ia(ns, i, j, k) * s(ns,1))
          tqr2 = tqr2 + (ia(ns, i, j, k) * s(ns,2))
          tqr3 = tqr3 + (ia(ns, i, j, k) * s(ns,3))
          !G = sum ia * omega
          ! same for the incindent radiation
          tga = tga + (ia(ns, i, j, k) * omeg(ns))
       end do
       !$acc end loop
       qr(i, j, k, 1) = tqr1
       qr(i, j, k, 2) = tqr2
       qr(i, j, k, 3) = tqr3
       ga(i, j, k) = tga
    end do; end do; end do
    !$acc end loop

    !$acc end parallel

    !$acc end data

    return

  end subroutine calcqr

end module fvdom
