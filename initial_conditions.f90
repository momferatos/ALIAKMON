!!$
!!$ (  / _ \\ \        /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!$Copyright (c) 2009-2020 Georgios Momferatos
module initial_conditions
  use types
  use parameters
  use data, only: zero, copy
  use numerics, only: fourier, truncate, rescale, project
  implicit none
  !
  ! subroutine for the initial conditions
  !
contains

  subroutine erandom_field(nn,fu)
    use data, only: wv
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT) :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !random field with particular energy spectrum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                             :: i,j,k,l
    real(rk)                                                :: kk,r1,r2
    complex(ck)                                             :: tmp


    fu = 0.0_rk
    do l=1,nn(4)
       do i=1,dim1(nn(1))-1,2
          do j=1,nn(2)
             do k=1,nn(3)
                kk=sqrt(real(wv(1_ik,i,j,k)**2+wv(2_ik,i,j,k)**2+&
                     &wv(3_ik,i,j,k)**2))
                call random_number(r1)
                call random_number(r2)
                !Set-up phase
                tmp=abs(r2)*exp(ii*2.0_rk*PI*r1)
                !Set up radius
                if(kk==0.0_rk) then
                   tmp=cmplx(0.0_rk,0.0_rk,ck)
                else
                   tmp=tmp*(kk**(-2.))*(kk**(-5./3.))
                end if

                tmp=sqrt(tmp)
                if (wv(1_ik,i,j,k).ne.0) tmp=tmp/2._rk

                fu(i,j,k,l)=real(tmp,rk)
                fu(i+1,j,k,l)=aimag(tmp)


             end do
          end do
       end do
    end do

    call truncate(nn,fu)

    return

  end subroutine erandom_field

  subroutine random_field(nn,fu,kmax)
    use data, only: trk1,trk2,trk3
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT)              :: fu
    real(rk), intent(IN)                                    :: kmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !random field with flat spectrum up to kmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                             :: i,j,k,l
    integer(ik)                                             :: ni,nj,nk,iii
    real(rk)                                                :: r1,r2
    complex(ck)                                             :: tmp

    ni=2.*nmodes
    nj=nmodes
    nk=nmodes


    do l=1,nn(4)
       do k=1,nn(3)
          do j=1,nn(2)
             do i=1,dim1(nn(1))-1,2
                iii=(i-1)/2+1
                if((trk1(iii)/=0.and.trk2(j)/=0.and.trk3(k)/=0).and.&
                     &(sqrt(real(trk1(iii)**2+trk2(j)**2+&
                     &trk3(k)**2,rk))<kmax)) then
                   !Set radius
                   call random_number(r1)
                   !Set phase
                   call random_number(r2)
                   tmp=abs(r1)*exp(2.0_rk*PI*r2)
                   fu(i,j,k,l)=real(tmp,rk)
                   fu(i+1,j,k,l)=aimag(tmp)
                else
                   fu(i,j,k,l)=0.0_rk
                   fu(i+1,j,k,l)=0.0_rk
                end if
             end do
          end do
       end do
    end do

    return

  end subroutine random_field


  subroutine set_initial_conditions(nn,u,fu)
    !use lagrangian, only: lset_initial_conditions
    use mpivars
    use fvdom, only: calcia
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT) :: u
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT) :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Initial conditions for the velocity, magnetic and scalar fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                             :: i, j, k, l
    
    !Set initial conditions for particles
!!$    if(PARTICLES) then
!!$       call lset_initial_conditions(np,x,vp)
!!$    end if

    !Set zero initial conditions, use only with forcing
    select case(initcond)
    case(zero_initcond)
       !zero initial conditions
       call zero(nn,fu)
    case(stochastic_initcond_flat)
       !Stochastic initial conditions with flat energy spectrum
       call random_field(nn,fu,2.0_rk)
    case(stochastic_initcond_with_spectrum)
       !Stochastic initial conditions with given energy spectrum (???)
       call erandom_field(nn,fu)
    case(orszag_tang_vortex)
       !Stochastic initial conditions based on the Orszag-Tang vortex
       call orszag_tang(nn,fu)
    case(abc)
       !Stochastic initial conditions based on the Arnold-Beltrami-Childress flow
       !Set-up large-scale component based on the ABC flow
       call abc_flow(nn,1_ik,3_ik,fu)
    case(taylor_green_vortex)
       call taylor_green(nn,fu)
    end select
        
    !Perform truncation
    call truncate(nn,fu)

    !Enforce incompressibility / zero magnetic field divergence
    call project(nn,fu)

    !Rescale to unit rms
    call rescale(nn,fu)


    !set magnetic field relative strength
    if(MHD) then
       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
          fu(i,j,k,nb1:nb3)=sqrt(BETA)*fu(i,j,k,nb1:nb3)
       end do ; end do ; end do
       !$omp end parallel do
    end if

    !If the passive scalar field is heated, set it to zero
    if(PASSIVE_SCALAR.and.HEATING) then
       !$omp parallel do
       do l=nsclf,nscll ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
          fu(i,j,k,l)=0.0_rk
       end do; end do ; end do ; end do
       !omp end parallel do
    end if

    !Copy Fourier arrays to physical-space arrays
    call copy(nn,u,fu)
    
    !Inverse Fourier transforms
    call fourier(nn,-1_ik,u)

    if(RADIATION) call calcia
    
    return

  end subroutine set_initial_conditions


  subroutine abc_flow(nn,k1,k2,u)
    use types
    use data, only: nu1,nu2,nu3,nb1,nb2,nb3,nsclf,nscll,rks1
    use mpivars
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    integer(ik), intent(IN)                                :: k1,k2
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)),intent(OUT)              :: u
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Arnol'd-Beltrami-Childress initial condition for wave-vectors with
    !radius k=k1-k2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer(ik)                                            :: i,j,k,kk
    real(rk)                                               :: a,b,c
    real(rk)                                               :: dx,dy,dz
    real(rk)                                               :: x,y,z,dzmpi
    real(rk), dimension(1:nn(4))                                 :: xx
    integer(ik), dimension(k1:k2, 1:nn(4))                   :: ii
    real(rk), dimension(1:3)                                 :: aa,r
    integer(ik)                                            :: nu,l


    !Arrays of indices used to mix the xyz components
    !for the velocity field
    ii(:,nu1)=(/1,2,3/)
    ii(:,nu2)=(/2,3,1/)
    ii(:,nu3)=(/3,1,2/)
    if(PASSIVE_SCALAR) then
       do l=nsclf,nscll
          call random_number(r)
          ii(:,l)=floor(4.0*r+1)
       end do
    end if
    !for the magnetic field
    if(MHD) then
       ii(:,nb1)=(/1,2,3/)
       ii(:,nb2)=(/1,2,3/)
       ii(:,nb3)=(/1,2,3/)
    end if
    
    !Set fields to zero
    call zero(nn,u)

    !Calculate steps
    dx=LBOX/real(n1,rk)
    dy=LBOX/real(gn2,rk)
    dz=LBOX/real(gn3,rk)
    !Calculate stochastic coefficients
    do kk=k1,k2
!!$       if(mpirank==MPIROOT) then
!!$          aa(1)=2.0*(rand()-0.5)
!!$          aa(2)=2.0*(rand()-0.5)
!!$          aa(3)=2.0*(rand()-0.5)
!!$       end if
!!$       !Broadcast   
!!$#ifdef _MPI_
!!$       call mpi_bcast(aa,3,MPIRK,MPIROOT,MPI_COMM_WORLD,mpierr)
!!$#endif
       aa(1)=0.1
       aa(2)=0.3
       aa(3)=0.4
       !Set-up stochastic coefficients
       a=aa(1)
       b=aa(2)
       c=aa(3)
       !Calculate large-scale ABC flow
       do nu=1,nn(4)
          do k=1,nn(3)
             z=(lkstart+k-1)*dz
             do j=1,nn(2)
                y=(ljstart+j-1)*dy
                do i=1,nn(1)
                   !Calculate coordinates
                   x=(i-1)*dx
                   xx(1)=a*sin(kk*z)+c*cos(kk*y)
                   xx(2)=b*sin(kk*x)+a*cos(kk*z)
                   xx(3)=c*sin(kk*y)+b*cos(kk*x)
                   !Calculate field
                   u(i,j,k,nu)=u(i,j,k,nu)+xx(ii(kk,nu))
                end do
             end do
          end do
       end do
    end do

    call fourier(nn,1_ik,u)

    !Set-up stochastic small-scale component
    call random_field(nn,rks1,2.0_rk)
    !Enforce incompressibility / zero magnetic field divergence
    call project(nn,rks1)
    !rescale to unit rms
    call rescale(nn,rks1)
    !$omp parallel do
    do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       !Add components
       u(i,j,k,l)=u(i,j,k,l)+rks1(i,j,k,l)
    end do; end do ; end do ; end do
    !$omp end parallel do

    call zero(nn,rks1)
    
    return

  end subroutine abc_flow

  subroutine taylor_green(nn,u)
    use hdf5_aliakmon, only: write_hdf5_file
    use types
    use data, only: nu1,nu2,nu3
    use mpivars
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)),intent(OUT)              :: u
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Arnol'd-Beltrami-Childress initial condition for wave-vectors with
    !radius k=k1-k2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer(ik)                                            :: i,j,k,l
    real(rk)                                               :: a,b,c
    real(rk)                                               :: dx,dy,dz
    real(rk)                                               :: x,y,z


    !Set fields to zero
    call zero(nn,u)

    !Calculate steps
    dx=LBOX/real(n1-1,rk)
    dy=LBOX/real(gn2-1,rk)
    dz=LBOX/real(gn3-1,rk)

    !Set-up stochastic coefficients
    !call random_number(a)
    !call random_number(b)
    a=0.2_rk
    b=0.3_rk
    c=1.0_rk-a-b
    !Calculate large-scale ABC flow
    do k=1,nn(3)
       z=(lkstart+k-1)*dz
       do j=1,nn(2)
          y=(ljstart+j-1)*dy
          do i=1,nn(1)
             x=(i-1)*dx
             !Calculate field
             u(i,j,k,nu1)=cos(a*x)*sin(b*y)*sin(c*z)
             u(i,j,k,nu2)=sin(a*x)*cos(b*y)*sin(c*z)
             u(i,j,k,nu3)=sin(a*x)*sin(b*y)*cos(c*z)
!!$             if(PASSIVE_SCALAR) then
!!$                do l=nsclf,nscll
!!$                   u(i,j,k,l)=cos(a*x)*sin(b*y)*sin(c*z)
!!$                end do
!!$             end if
          end do
       end do
    end do
    
    call fourier(nn,1_ik,u)

    return

  end subroutine taylor_green

  subroutine orszag_tang(nn,u)
    use types
    use data, only: nu1,nu2,nu3,rks1
    use mpivars
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)),intent(OUT)              :: u
    !
    ! Orszag-Tang vortex initial condition
    !
    integer(ik) :: i,j,k,l
    real(rk) :: dx,dy,dz,xx,yy,zz
    real(rk)                                                :: OT_RAND

    OT_RAND=1.0e-1_rk
    
    !Stochastic component
    call random_field(nn,rks1,2.0_rk)
    !rescale to unit rms
    call rescale(nn,rks1)
    call fourier(nn,-1_ik,rks1)

    !Orszag-Tang vortex large-scale component
    dx=LBOX/real(nn(1)-1,rk)
    dy=LBOX/real(gn2-1,rk)
    dz=LBOX/real(gn3-1,rk)

    do i=1,nn(1)
       do j=1,nn(2)
          do k=1,nn(3)
             xx=(i-1)*dx
             yy=(ljstart+j-1)*dy
             zz=(lkstart+k-1)*dz
             u(i,j,k,nu1)=-2.0_rk*sin(yy)
             u(i,j,k,nu2)=2.0_rk*sin(xx)
             u(i,j,k,nu3)=sin(xx)+sin(yy)
             if(MHD) then
                u(i,j,k,nb1)=-2.0_rk*sin(2.0_rk*yy)+sin(zz)
                u(i,j,k,nb2)=2.0_rk*sin(xx)+sin(zz)
                u(i,j,k,nb3)=sin(xx)+sin(yy) 
             end if
          end do
       end do
    end do

    !$omp parallel do
    do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       !Add components
       u(i,j,k,l)=u(i,j,k,l)+OT_RAND*rks1(i,j,k,l)
    end do; end do ; end do ; end do
    !$omp end parallel do

    call zero(nn,rks1)
    call fourier(nn,1_ik,u)

    return

  end subroutine orszag_tang
     
end module initial_conditions

