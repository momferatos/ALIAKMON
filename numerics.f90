!!$     ___ __                                       
!!$ (  / _ \\ \        /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!!$    Copyright (c) 2009-2020 Georgios Momferatos
!!$
module numerics
  use types
  use parameters
  use data, only: zero, copy
  implicit none
  !
  ! physics & numerics module
  !
contains

  subroutine msvalue(nn,fu,msv)
    use mpivars
    use data, only: rmsarr,nu1,nu3,nb1,nb3
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rk), dimension(1:nn(4)), intent(out)                 :: msv
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)                   :: fu
    !
    ! calculate mean-square value for vector array
    !
    integer(ik)                                                 :: i,j,k,l
    real(rk)                                                    :: fac,tmp

    
    !if(.not.allocated(rms)) allocate(rms(1:nn(4)))

    call copy(nn,rmsarr,fu)

    call fourier(nn,-1_ik,rmsarr)
    
    fac= 1.0_rk / real(nn(1)*gn2*gn3, rk)
    
    do l=1,nn(4)
       tmp=0.0_rk
       !$omp parallel do reduction(+:tmp)
       do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          tmp=tmp+fac*rmsarr(i,j,k,l)**2
       end do; end do ; end do
       !$omp end parallel do
       msv(l)=tmp
    end do
    
#ifdef _MPI_
    sbuf(1:nn(4))=msv(1:nn(4))
    call mpi_allreduce(sbuf,rbuf,int(nn(4)),MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
    msv(1:nn(4))=rbuf(1:nn(4))
#endif

    tmp=msv(nu1)+msv(nu2)+msv(nu3)
    msv(nu1:nu3)=tmp
    
!!$    if(PASSIVE_SCALAR) then
!!$       msv(nsclf:nscll)=rms(nsclf:nscll)
!!$    end if

    if(MHD) then
       tmp=msv(nb1)+msv(nb2)+msv(nb3)
       msv(nb1:nu3)=tmp
    end if

    call zero(nn,rmsarr)

    return

  end subroutine msvalue

  function mean_value(nn,fu) result(meanval)
    use mpivars
    use data, only: rmsarr,nu1,nu3,nsclf,nscll,nb1,nb3
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rk), dimension(1:nn(4))                            :: meanval
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)                   :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !calculate mean-square value for vector array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                                 :: i,j,k,l
    real(rk), dimension(:), allocatable                         :: fmean
    real(rk)                                                    :: fac,tmp

    if(.not.allocated(fmean)) allocate(fmean(1:nn(4)))
    !Calculate scale factor
    fac=1.0_rk / real(nn(1)*gn2*gn3,rk)

    !Copy Fourier-space arrays to physical-space arrays
    call copy(nn,rmsarr,fu)
    
    !Inverse Fourier transform
    call fourier(nn,-1_ik,rmsarr)

    !Calculate root-mean-square
    fmean(:)=0.0_rk    
    do l=1,nn(4)
       tmp=0.0_rk
       !$omp parallel do reduction(+:tmp)
       do k=1,nn(3)
          do j=1,nn(2)
             do i=1,nn(1)
                tmp=tmp+fac*rmsarr(i,j,k,l)
             end do
          end do
       end do
       !$omp end parallel do
       fmean(l)=tmp
    end do


    !Set temporary arrays back to zero
    call zero(nn,rmsarr)
    
    meanval(nu1:nu3)= fmean(nu1:nu3)
    if(PASSIVE_SCALAR) then
       meanval(nsclf:nscll)=fmean(nsclf:nscll)
    end if
    if(MHD) then
       meanval(nb1:nb3)=fmean(nb1:nb3)
    end if

    !Reduce root-mean-square
#ifdef _MPI_
    sbuf(1:nn(4))=meanval(1:nn(4))
    call mpi_allreduce(sbuf,rbuf,int(nn(4)),MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
    meanval(1:nn(4))=rbuf(1:nn(4))
#endif

    return

  end function mean_value


  subroutine project(nn,fu)
    use data, only:wv,nu1,nb1,isactive
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN OUT)               :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !project velocity/magnetic field to solenoidal subspace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                                 :: i,j,k,l,iii
    real(rk)                                                    :: ksq
    complex(ck)                                                 :: tmp1,p
    complex(ck), dimension(1:3)                                   :: fpgrad
    integer(ik)                                                 :: nu
    integer(ik), dimension(1:2)                                   :: ifield
    integer(ik)                                                 :: nfield


    !Do nothing if we're solving the Burgers equation
    if(BURGERS) return

    nfield=1
    ifield(1)=nu1
    if(MHD) then
       nfield=2
       ifield(2)=nb1
    end if

    !Main projection loop
    !$omp parallel do  private(nu,iii,ksq,tmp1,p,fpgrad)
    do l=1,nfield
       nu=ifield(l)
       do k=1,nn(3)
          do j=1,nn(2)
             do i=1,dim1(nn(1))-1,2
                if(isactive(i,j,k)) then
                   !Wave-vector modulus squared
                   ksq=wv(1_ik,i,j,k)**2+wv(2_ik,i,j,k)**2+wv(3_ik,i,j,k)**2
                   !Pressure term
                   if(ksq==0.0_rk) then
                      p=cmplx(0.0_rk,0.0_rk,ck)
                   else
                      tmp1= ii*wv(1_ik,i,j,k)*cmplx(fu(i,j,k,nu  ),fu(i+1,j,k,nu  ),ck)+&
                           &ii*wv(2_ik,i,j,k)*cmplx(fu(i,j,k,nu+1),fu(i+1,j,k,nu+1),ck)+&
                           &ii*wv(3_ik,i,j,k)*cmplx(fu(i,j,k,nu+2),fu(i+1,j,k,nu+2),ck)
                      p=tmp1/ksq
                   end if
                   !Pressure gradient
                   fpgrad(1)=ii*wv(1_ik,i,j,k)*p
                   fpgrad(2)=ii*wv(2_ik,i,j,k)*p
                   fpgrad(3)=ii*wv(3_ik,i,j,k)*p
                   !Projection
                   fu(i,j,k,nu:nu+2)  =fu(i,  j,k,nu:nu+2)+real( fpgrad(:),rk)
                   fu(i+1,j,k,nu:nu+2)=fu(i+1,j,k,nu:nu+2)+aimag(fpgrad(:))
                else
                   fu(i,j,k,:)  =0.0_rk
                   fu(i+1,j,k,:)=0.0_rk
                end if

             end do
          end do
       end do
    end do
    !$omp end parallel do

    return

  end subroutine project



  subroutine timestep(nn,fu,dt)
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN OUT)  :: fu
    real(rk), intent(IN)                                        :: dt
!!!!!!!!!!!!!!!!!!!!!
    !timestep advancement
!!!!!!!!!!!!!!!!!!!!!

    if(INTEGRATION_METHOD==MRUNGE_KUTTA2) then
       call runge_kutta2(nn,fu,dt)
    else if(INTEGRATION_METHOD==MRUNGE_KUTTA4) then
       call runge_kutta4(nn,fu,dt)
    else 
       stop 'timestep:: Invalid integration method.'
    end if

    return

  end subroutine timestep



  subroutine runge_kutta2(nn,fu,dt)
    use data, only: wv,rhs,rks1,isactive
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN OUT)               :: fu
    real(rk), intent(IN)                                        :: dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Second-order Runge-Kutta time integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                                 :: i,j,k,l

    !Get right-hand-side terms
    call right_hand_side(nn,rhs,fu)

    !$omp parallel do
    !Integration loop for the first step
    do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       if(isactive(i,j,k)) then
          rks1(i,j,k,l)=fu(i,j,k,l)+0.5_rk*dt*rhs(i,j,k,l)
       else
          rks1(i,j,k,l)=0.0_rk
       end if
    end do; end do ; end do ; end do
    !$omp end parallel do

    !Get right-hand-side terms
    call right_hand_side(nn,rhs,rks1)

    !$omp parallel do
    !Integration loop for the second step
    do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       if(isactive(i,j,k)) then
          fu(i,j,k,l)=fu(i,j,k,l)+dt*rhs(i,j,k,l)
       else
          fu(i,j,k,l)=0.0_rk
       end if
    end do; end do ; end do ; end do
    !$omp end parallel do

    return

  end subroutine runge_kutta2



  subroutine runge_kutta4(nn,fu,dt)
    use data, only: wv,rhs,rks1,rks2,isactive
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN OUT)               :: fu
    real(rk), intent(IN)                                        :: dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Fourth-order Runge-Kutta time integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                                 :: i,j,k,l
    real(rk), dimension(1:4), parameter                           :: rkc=&
         &(/1._rk/6._rk,1._rk/3._rk,1._rk/3._rk,1._rk/6._rk/)

    !Get right-hand-side terms
    call right_hand_side(nn,rhs,fu)

    !$omp parallel do
    !Integration loop for the first step
    do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       if(isactive(i,j,k)) then
          rks2(i,j,k,l)=rkc(1)*dt*rhs(i,j,k,l)
          rks1(i,j,k,l)=fu(i,j,k,l)+0.5_rk*dt*rhs(i,j,k,l)
       else
          rks1(i,j,k,l)=0.0_rk
          rks2(i,j,k,l)=0.0_rk
       end if
    end do; end do ; end do ; end do
    !$omp end parallel do

    !Get right-hand-side terms
    call right_hand_side(nn,rhs,rks1)

    !$omp parallel do
    !integration loop for the second step
    do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       if(isactive(i,j,k)) then
          rks2(i,j,k,l)=rks2(i,j,k,l)+rkc(2)*dt*rhs(i,j,k,l)
          rks1(i,j,k,l)=fu(i,j,k,l)+0.5_rk*dt*rhs(i,j,k,l)
       else
          rks1(i,j,k,l)=0.0_rk
          rks2(i,j,k,l)=0.0_rk
       end if
    end do; end do ; end do ; end do
    !$omp end parallel do

    !Get right-hand-side terms
    call right_hand_side(nn,rhs,rks1)

    !$omp parallel do
    !integration loop for the third step
    do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       if(isactive(i,j,k)) then
          rks2(i,j,k,l)=rks2(i,j,k,l)+rkc(3)*dt*rhs(i,j,k,l)
          rks1(i,j,k,l)=fu(i,j,k,l)+dt*rhs(i,j,k,l)
       else
          rks1(i,j,k,l)=0.0_rk
          rks2(i,j,k,l)=0.0_rk
       end if
    end do; end do ; end do ; end do
    !$omp end parallel do


    !Get right-hand-side terms
    call right_hand_side(nn,rhs,rks1)

    !$omp parallel do
    !integration loop for the fourth step
    do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       if(isactive(i,j,k)) then
          fu(i,j,k,l)=fu(i,j,k,l)+rks2(i,j,k,l)+rkc(4)*dt*rhs(i,j,k,l)
       else
          fu(i,j,k,l)=0.0_rk
       end if
    end do; end do ; end do ; end do
    !$omp end parallel do

    return

  end subroutine runge_kutta4


  subroutine right_hand_side(nn,fnl,fu)
    use data, only: isactive,u,du,psu,scratch,fnls,nu1,nu2,nu3,nb1,nb2,nb3,nsclf,nscll,ad,&
         &fsclgrads,fsclgrad,wv,arr_en_1
    use mpivars
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT) :: fnl
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN) :: fu 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Non-linear terms in rotational form
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                             :: i, j, k,l,m
    complex(ck), dimension(:),allocatable :: ff, fvisc
    real(rk) :: ksq

    if(.not.allocated(fvisc)) allocate(fvisc(1:nn(4)))
    if(.not.allocated(ff)) allocate(ff(1:nn(4)))

    call compute_hydro

    if(PASSIVE_SCALAR) call compute_passive_scalar
    
    if(MHD) call compute_mhd

    ! handle Patterson-Orszag deliasing
    if(DEALIASING==PATTERSON_ORSZAG) then
       call shift(nn,-1_ik,fnls)
       !$omp parallel do
       do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
          if(isactive(i,j,k)) then
             fnl(i,j,k,l)=0.5_rk*(fnl(i,j,k,l)+fnls(i,j,k,l))
          end if
       end do; end do ; end do ; end do
       !$omp end parallel do
    end if

    ! truncate non-linear terms to remove aliasing errors
    call truncate(nn,fnl)

    ! project to zero divergence
    if(.not.BURGERS) call project(nn,fnl)

    ! compute diffusive terms
    call compute_diffusive_terms


    ! set auxiliary arrays back to zero
    call zero(nn,u)
    call zero(nn,scratch)

    if(DEALIASING==PATTERSON_ORSZAG) then
       call zero(nn,du)
       call zero(nn,psu)
       call zero(nn,fnls)
    end if
    
    if(PASSIVE_SCALAR) call zero(nn,fsclgrads)

    return
    
  contains
    
    subroutine compute_hydro
      use types
      implicit none
      
      !initialize arrays to zero
      call zero(nn,fnl)
      call zero(nn,scratch)
      !Get velocity field in Fourier space
      call copy(nn,u,fu)

      !Calculate the vorticity field in Fourier space 
      call curl(nn,scratch,fu,nu1)

      !Patterson-Orszag dealiasing
      if(DEALIASING==PATTERSON_ORSZAG) then

         !set phase-shifted non-linear terms to zero
         call zero(nn,fnls)

         !Copy to auxiliary arrays
         call copy(nn,psu,fu)
         call copy(nn,du,scratch)

         !Perform phase shifts
         call shift(nn,1_ik,psu)
         call shift(nn,1_ik,du)

         !Inverse Fourier transforms
         call fourier(nn,-1_ik,psu,nfs=nu1,nfe=nu3)
         call fourier(nn,-1_ik,du,nfs=nu1,nfe=nu3)

      end if

      !Inverse Fourier transforms for the non-phase shifted fields
      call fourier(nn,-1_ik,u,nfs=nu1,nfe=nu3)
      call fourier(nn,-1_ik,scratch,nfs=nu1,nfe=nu3)

      maxvort=0.0_rk
      maxvel=0.0_rk
      !$omp parallel do reduction(max:maxvort,maxvel)
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         maxvel=max(maxvel,real(sqrt(u(i,j,k,nu1)**2+u(i,j,k,nu2)**2+u(i,j,k,nu3)**2),rk))
         maxvort=max(maxvort,real(sqrt(scratch(i,j,k,nu1)**2+scratch(i,j,k,nu2)**2+scratch(i,j,k,nu3)**2),rk))
         !rhs= u x omega
         fnl(i,j,k,nu1)= u(i,j,k,nu2)*scratch(i,j,k,nu3) - u(i,j,k,nu3)*scratch(i,j,k,nu2)
         fnl(i,j,k,nu2)= u(i,j,k,nu3)*scratch(i,j,k,nu1) - u(i,j,k,nu1)*scratch(i,j,k,nu3)
         fnl(i,j,k,nu3)= u(i,j,k,nu1)*scratch(i,j,k,nu2) - u(i,j,k,nu2)*scratch(i,j,k,nu1)
      end do; end do ; end do   
      !$omp end parallel do

#ifdef _MPI_
      sbuf(1)=MAXVEL
      sbuf(2)=MAXVORT
      call mpi_allreduce(sbuf,rbuf,2,MPIRK,MPI_MAX,MPI_COMM_WORLD,mpierr)
      MAXVEL=rbuf(1)
      MAXVORT=rbuf(2)
#endif


      if(DEALIASING==PATTERSON_ORSZAG) then
         !$omp parallel do
         !Non-linear term
         do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
            !rhs= u x omega
            fnls(i,j,k,nu1)=psu(i,j,k,nu2)*du(i,j,k,nu3) - psu(i,j,k,nu3)*du(i,j,k,nu2)
            fnls(i,j,k,nu2)=psu(i,j,k,nu3)*du(i,j,k,nu1) - psu(i,j,k,nu1)*du(i,j,k,nu3)
            fnls(i,j,k,nu3)=psu(i,j,k,nu1)*du(i,j,k,nu2) - psu(i,j,k,nu2)*du(i,j,k,nu1)
         end do; end do ; end do   
         !$omp end parallel do
      end if


      call fourier(nn,1_ik,fnl,nu1)

      if(DEALIASING==PATTERSON_ORSZAG) call fourier(nn,1_ik,fnls,nu1)
      
      return

    end subroutine compute_hydro

    subroutine compute_passive_scalar
      use types
      implicit none

      !Set passive scalar components to zero
      !$omp parallel
      !$omp do
      do l=nsclf,nscll ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
         fnl(i,j,k,l)=0.0_rk
      end do; end do ; end do ; end do
      !$omp end do

      !$omp do
      do m=nu1,nu3 ; do l=nsclf,nscll ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
         fsclgrads(i,j,k,l,m)=0.0_rk
      end do; end do ; end do ; end do ; end do
      !$omp end do
      !$omp end parallel

      !compute passive scalar gradients
      call gradient(nn,fsclgrads,fu,nsclf,nscll)
      do l=nu1,nu3
         fsclgrad=>fsclgrads(:,:,:,:,l)
         call fourier(nn,-1_ik,fsclgrad,nsclf,nscll)
      end do
      !$omp parallel do
      do l=nsclf,nscll ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         fnl(i,j,k,l)=u(i,j,k,nu1)*fsclgrads(i,j,k,l,nu1)+u(i,j,k,nu2)*fsclgrads(i,j,k,l,nu2)+&
              &u(i,j,k,nu3)*fsclgrads(i,j,k,l,nu3)
      end do; end do ; end do ; end do
      !$omp end parallel do
      call fourier(nn,1_ik,fnl,nsclf,nscll)

      if(DEALIASING==PATTERSON_ORSZAG) then
         !$omp parallel
         !$omp do
         do l=nsclf,nscll ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
            fnls(i,j,k,l)=0.0_rk
         end do; end do ; end do ; end do
         !$omp end do

         !$omp do
         do m=nu1,nu3 ; do l=nsclf,nscll ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
            fsclgrads(i,j,k,l,m)=0.0_rk
         end do; end do ; end do ; end do ; end do
         !$omp end do
         !$omp end parallel 
         call gradient(nn,fsclgrads,fu,nsclf,nscll)
         do l=nu1,nu3
            fsclgrad=>fsclgrads(:,:,:,:,l)
            call shift(nn,1_ik,fsclgrad)
            call fourier(nn,-1_ik,fsclgrad,nsclf,nscll)
         end do
         !$omp parallel do
         do l=nsclf,nscll ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
            fnls(i,j,k,l)=psu(i,j,k,nu1)*fsclgrads(i,j,k,l,nu1)+psu(i,j,k,nu2)*fsclgrads(i,j,k,l,nu2)+&
                 &psu(i,j,k,nu3)*fsclgrads(i,j,k,l,nu3)
         end do; end do ; end do ; end do
         !$omp end parallel do
         call fourier(nn,1_ik,fnls,nsclf,nscll)
      end if
      
      return

    end subroutine compute_passive_scalar

    subroutine compute_mhd
      use types
      implicit none

      call zero(nn,scratch)
      call copy(nn,u,fu)
      call fourier(nn,-1_ik,u,nfs=nu1,nfe=nu3)
      call fourier(nn,-1_ik,u,nfs=nb1,nfe=nb3)
      !compute j x b
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         ! rhs = u x b
         scratch(i,j,k,nb1) = u(i,j,k,nu2)*u(i,j,k,nb3) - u(i,j,k,nu3)*u(i,j,k,nb2)
         scratch(i,j,k,nb2) = u(i,j,k,nu3)*u(i,j,k,nb1) - u(i,j,k,nu1)*u(i,j,k,nb3)
         scratch(i,j,k,nb3) = u(i,j,k,nu1)*u(i,j,k,nb2) - u(i,j,k,nu2)*u(i,j,k,nb1)
      end do; end do ;  end do
      !$omp end parallel do


      call fourier(nn,1_ik,scratch,nfs=nb1,nfe=nb3)
      !Forward Fourier transforms
      call curl(nn,fnl,scratch,nb1)

      !Patterson-Orszag dealiasing, Magnetohydrodynamics
      if(DEALIASING==PATTERSON_ORSZAG) then
         call zero(nn,du)
         call copy(nn,psu,fu)
         call shift(nn,1_ik,psu)
         call fourier(nn,-1_ik,psu,nfs=nu1,nfe=nu3)
         call fourier(nn,-1_ik,psu,nfs=nb1,nfe=nb3)
         !$omp parallel do
         do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
            ! rhs = (u x b)
            du(i,j,k,nb1) = psu(i,j,k,nu2)*psu(i,j,k,nb3) - psu(i,j,k,nu3)*psu(i,j,k,nb2)
            du(i,j,k,nb2) = psu(i,j,k,nu3)*psu(i,j,k,nb1) - psu(i,j,k,nu1)*psu(i,j,k,nb3)
            du(i,j,k,nb3) = psu(i,j,k,nu1)*psu(i,j,k,nb2) - psu(i,j,k,nu2)*psu(i,j,k,nb1)
         end do; end do ; end do
         !$omp end parallel do
         call fourier(nn,1_ik,du,nfs=nb1,nfe=nb3)
         !Forward Fourier transforms
         call curl(nn,fnls,du,nb1)
      end if

      !Get the Lorentz force
      call lorentz_force(nn,arr_en_1,ad,fu,u)
      !Add-up to get the final terms
      !$omp parallel do
      do l=nu1,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
         fnl(i,j,k,l)=fnl(i,j,k,l)+arr_en_1(i,j,k,l)
      end do; end do ; end do ; end do
      !$omp end parallel do
      if(AMB_DIFF) then 
         !$omp parallel do
         do l=nb1,nb3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
            fnl(i,j,k,l)=fnl(i,j,k,l)+ad(i,j,k,l)
         end do; end do ; end do ; end do
         !$omp end parallel do
      end if

      return

    end subroutine compute_mhd

    subroutine compute_diffusive_terms
      use types
      implicit none
      !$omp parallel do private(ksq,ff,fvisc)
      do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1)) 
         if(isactive(i,j,k)) then
            call forcing_rhs(nn,i,j,k,fu,ff)
            ksq=wv(1_ik,i,j,k)**2+wv(2_ik,i,j,k)**2+wv(3_ik,i,j,k)**2
            fvisc(l)=-ksq*visc(l)*fu(i,j,k,l)
            !Add to get the right-hand-side
            scratch(i,j,k,l)=fnl(i,j,k,l)+fvisc(l)+ff(l)
         else
            scratch(i,j,k,l)=0.0_rk
         end if
         fnl(i,j,k,l)=scratch(i,j,k,l)
      end do; end do ; end do ; end do 
      !$omp end parallel do
      return
    end subroutine compute_diffusive_terms
    
    
  end subroutine right_hand_side

  subroutine forcing_rhs(nn,i,j,k,fu,ff)
    use types
    use data, only:trk1,trk2,trk3
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    integer(ik), intent(IN)  :: i,j,k
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN) :: fu
    complex(ck), dimension(1:nn(4)), intent(OUT) :: ff
!!!!!!!!!!!!!!!!!
    !The forcing term
!!!!!!!!!!!!!!!!!
    logical                  :: condition
    integer(ik)              :: iii,l

    if(.not.FORCED) then
       ff(:)=cmplx(0.0_rk,0_rk,ck)
       return
    end if

    iii=(i-1)/2+1
    !Forcing condition in Fourier space
    condition=(trk1(iii)/=0.and.trk2(j)/=0.and.trk3(k)/=0).and.&
         &(sqrt(real(trk1(iii)**2+trk2(j)**2+trk3(k)**2,rk))<KFORCING)

    if(condition) then
!!$       !Toy stochastic forcing
!!$       if(.not.VARIABLE_FORCING) then
!!$          if(FORCED) then
!!$             call random_number(r)
!!$             call random_number(ph)
!!$             ff(nu1)=FSCALE*r*exp(ii*2*PI*ph)
!!$
!!$             call random_number(r)
!!$             call random_number(ph)
!!$             ff(nu2)=FSCALE*r*exp(ii*2*PI*ph)
!!$
!!$             call random_number(r)
!!$             call random_number(ph)
!!$             ff(nu3)=FSCALE*r*exp(ii*2*PI*ph)
!!$          end if
!!$
!!$          !Passive scalar
!!$          if(PASSIVE_SCALAR.and.FORCED_PASSIVE_SCALAR) then
!!$             do l=nsclf,nscll
!!$                call random_number(r)
!!$                call random_number(ph)
!!$                ff(l)=FSCALE*r*exp(ii*2*PI*ph)
!!$             end do
!!$          end if
!!$
!!$          !Magnetohydrodynamics
!!$          if(MHD.and.FORCED_MHD) then
!!$             call random_number(r)
!!$             call random_number(ph)
!!$             ff(nb1)=FSCALE*r*exp(ii*2*PI*ph)
!!$
!!$             call random_number(r)
!!$             call random_number(ph)
!!$             ff(nb2)=FSCALE*r*exp(ii*2*PI*ph)
!!$
!!$             call random_number(r)
!!$             call random_number(ph)
!!$             ff(nb3)=FSCALE*r*exp(ii*2*PI*ph)
!!$          end if
!!$
!!$          !Kaneda et al. (2004) forcing
!!$       else if(VARIABLE_FORCING) then
       do l=1,nn(4)
          ff(l)=fscale(l)*cmplx(fu(i,j,k,l),fu(i+1,j,k,l),ck)
       end do
    else
       ff(:)=cmplx(0.0_rk,0.0_rk,ck)
!!$       ff(nu2)=cmplx(0.0_rk,0.0_rk,ck)
!!$       ff(nu3)=cmplx(0.0_rk,0.0_rk,ck)
!!$
!!$       if(PASSIVE_SCALAR) ff(nscl)=cmplx(0.0_rk,0.0_rk,ck)
!!$
!!$       if(MHD) then
!!$          ff(nb1)=cmplx(0.0_rk,0.0_rk,ck)
!!$          ff(nb2)=cmplx(0.0_rk,0.0_rk,ck)
!!$          ff(nb3)=cmplx(0.0_rk,0.0_rk,ck)
!!$       end if

    end if

    return

  end subroutine forcing_rhs

  subroutine lorentz_force(nn,lf,ad,fu,u)
    use types
    use mpivars
    use data,only:scratch,psu,du,nb1,nb2,nb3,nu1,nu2,nu3
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT)    :: lf
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT)    :: ad
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)     :: fu
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(INOUT)  :: u
!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate Lorentz force
!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                                :: i,j,k,l

    call zero(nn,lf)

    !Set auxiliary arrays used for dealiasing to zero
    if(DEALIASING==PATTERSON_ORSZAG) then
       call zero(nn,psu)
       call zero(nn,du)
    end if

    !Copy magnetic field
    call copy(nn,u,fu)
    !Get magnetic field in physical space
    call fourier(nn,-1_ik,u,nfs=nb1,nfe=nb3)

    !Calculate current density
    call curl(nn,scratch,fu,nb1)

    !Get current density in physical space
    call fourier(nn,-1_ik,scratch,nfs=nb1,nfe=nb3)

    !Patterson-Orszag dealiasing
    if(DEALIASING==PATTERSON_ORSZAG) then
       !Copy to auxiliary arrays
       call copy(nn,psu,u)
       call copy(nn,du,scratch)

       !Perform phase-shifts
       call shift(nn,1_ik,psu)
       call shift(nn,1_ik,du)

       !Get phase-shifted arrays in physical space
       call fourier(nn,-1_ik,psu,nb1)
       call fourier(nn,-1_ik,du,nb1)

    end if

    !Ambipolar diffusion and Hall effect
    if(AMB_DIFF.or.HALL) then

    end if

    maxj=0.0_rk
    maxlf=0.0_rk
    !$omp parallel do reduction(max:maxj,maxlf)
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       ! lf = j x b
       lf(i,j,k,nu1)=scratch(i,j,k,nb2)*u(i,j,k,nb3)-scratch(i,j,k,nb3)*u(i,j,k,nb2)
       lf(i,j,k,nu2)=scratch(i,j,k,nb3)*u(i,j,k,nb1)-scratch(i,j,k,nb1)*u(i,j,k,nb3)
       lf(i,j,k,nu3)=scratch(i,j,k,nb1)*u(i,j,k,nb2)-scratch(i,j,k,nb2)*u(i,j,k,nb1)

       maxj=max(maxj,real(sqrt(scratch(i,j,k,nb1)**2+scratch(i,j,k,nb2)**2+scratch(i,j,k,nb3)**2),rk))
       maxlf=max(maxlf,real(sqrt(lf(i,j,k,nu1)**2+lf(i,j,k,nu2)**2+lf(i,j,k,nu3)**2),rk))

       if(DEALIASING==PATTERSON_ORSZAG) then
          !Phase-shifted arrays
          ! scratch = j x b
          ad(i,j,k,nu1)=du(i,j,k,nb2)*psu(i,j,k,nb3)-du(i,j,k,nb3)*psu(i,j,k,nb2)
          ad(i,j,k,nu2)=du(i,j,k,nb3)*psu(i,j,k,nb1)-du(i,j,k,nb1)*psu(i,j,k,nb3)
          ad(i,j,k,nu3)=du(i,j,k,nb1)*psu(i,j,k,nb2)-du(i,j,k,nb2)*psu(i,j,k,nb1)
       end if


       !Ambipolar diffusion loop
       if(AMB_DIFF) then
          !Set auxiliary arrays to zero
          call zero(nn,ad)
          if(DEALIASING==PATTERSON_ORSZAG) then
             call zero(nn,psu)
          end if
          !Cubic non-linear terms
          ad(i,j,k,nu1)=AD_COEFF*(lf(i,j,k,nb2)*u(i,j,k,nb3)-lf(i,j,k,nb3)*u(i,j,k,nb2))
          ad(i,j,k,nu2)=AD_COEFF*(lf(i,j,k,nb3)*u(i,j,k,nb1)-lf(i,j,k,nb1)*u(i,j,k,nb3))
          ad(i,j,k,nu3)=AD_COEFF*(lf(i,j,k,nb1)*u(i,j,k,nb2)-lf(i,j,k,nb2)*u(i,j,k,nb1))
          !Phase-shifted terms
          if(DEALIASING==PATTERSON_ORSZAG) then
             psu(i,j,k,nu1)=AD_COEFF*(du(i,j,k,nb2)*u(i,j,k,nb3)-du(i,j,k,nb3)*u(i,j,k,nb2))
             psu(i,j,k,nu2)=AD_COEFF*(du(i,j,k,nb3)*u(i,j,k,nb1)-du(i,j,k,nb1)*u(i,j,k,nb3))
             psu(i,j,k,nu3)=AD_COEFF*(du(i,j,k,nb1)*u(i,j,k,nb2)-du(i,j,k,nb2)*u(i,j,k,nb1))
          end if
       end if


       !Main loop for the Hall term
       if(HALL) then
          do l=nu1,nu3
             ad(i,j,k,l)=ad(i,j,k,l)-HALL_COEFF*lf(i,j,k,l)
          end do
          if(DEALIASING==PATTERSON_ORSZAG) then
             do l=nu1,nu3
                psu(i,j,k,l)=psu(i,j,k,l)-HALL_COEFF*du(i,j,k,l)
             end do
          end if
       end if

    end do; end do ; end do
    !$omp end parallel do

#ifdef _MPI_
    sbuf(1)=maxj
    sbuf(2)=maxlf
    call mpi_allreduce(sbuf,rbuf,2,MPIRK,MPI_MAX,MPI_COMM_WORLD,mpierr)
    maxj=rbuf(1)
    maxlf=rbuf(2)
#endif

    !Forward Fourier transforms
    call fourier(nn,1_ik,lf,nfs=nu1,nfe=nu3)

    if(AMB_DIFF) then
       call fourier(nn,1_ik,ad,nfs=nu1,nfe=nu3)
    end if


    !Same for the phase-shifted arrays
    if(DEALIASING==PATTERSON_ORSZAG) then
       call fourier(nn,1_ik,ad,nfs=nu1,nfe=nu3)
       call shift(nn,-1_ik,ad)
       if(AMB_DIFF) then
          call fourier(nn,1_ik,psu,nu1,nu3)
          call shift(nn,-1_ik,psu)
       end if
       !Get the total (average) terms
       !$omp parallel do
       do l=nu1,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          lf(i,j,k,l) = 0.5_rk*(lf(i,j,k,l) + ad(i,j,k,l))
          if(AMB_DIFF) then
             ad(i,j,k,l)=0.5_rk*(ad(i,j,k,l)+psu(i,j,k,l))
          end if
       end do; end do ; end do ; end do
    end if

    !Truncate the Lorentz force
    call truncate(nn,lf)

    !Truncate the ambipolar diffusion terms
    if(AMB_DIFF.or.HALL) then
       call curl(nn,ad,ad,nu1)
       call truncate(nn,ad)
    end if


    return

  end subroutine lorentz_force

  subroutine mean_dissipation(nn,fu,emeans)
    use data, only: scratch,fsclgrad,fsclgrads,nu1,nb1
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)  :: fu
    real(rk), dimension(1:nn(4)), intent(out)                                  :: emeans
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate mean energy dissipation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: i,j,k,l
    real(rk) :: fac,tmp

    fac=1.0_rk/real(nn(1)*gn2*gn3,rk)

    call zero(nn,scratch)

    !Get curl of the velocity field
    call curl(nn,scratch,fu,nu1)
    if(MHD) call curl(nn,scratch,fu,nb1)

    emeans(:)=0.0_rk
    !Calculate mean-square
    call msvalue(nn,scratch,emeans)

    if(PASSIVE_SCALAR) then
       call gradient(nn,fsclgrads,fu)
       do l=nu1,nu3
          fsclgrad=>fsclgrads(:,:,:,:,l)
          call fourier(nn,-1_ik,fsclgrad)
       end do

       do l=nsclf,nscll
          tmp=0.0_rk
          !$omp parallel do reduction(+:tmp)
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             tmp = tmp + fac*(fsclgrads(i,j,k,l,nu1)**2 + fsclgrads(i,j,k,l,nu2)**2 +&
                  &fsclgrads(i,j,k,l,nu3)**2)
          end do; end do; end do 
          !$omp end parallel do
          emeans(l)=tmp
       end do
       
       call zero(nn,fsclgrads)
       
    end if

    emeans(:)=visc(:)*emeans(:)

    call zero(nn,scratch)
    
    return

  end subroutine mean_dissipation


  function mean_kinetic_helicity_dissipation(nn,fu) result(mkhdis)
    use data, only: scratch,rmsarr,rks1
    implicit none
    real(rk)                                                :: mkhdis
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)  :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate mean kinetic helicity dissipation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                             :: i,j,k
    integer(i8b)                                            :: l
    real(rk), dimension(:), allocatable                     :: tmp                     
    if(.not.allocated(tmp)) allocate(tmp(1:nn(4)))


    !Get vorticity
    call curl(nn,rmsarr,fu,nu1)
    !Get curl of vorticity
    call curl(nn,rks1,rmsarr,nu1)
    !Inverse Fourier transform
    call fourier(nn,-1_ik,rmsarr)

    call fourier(nn,-1_ik,rks1)

    !Perform inner product
    call inner_product(nn,scratch,nu1,rmsarr,rks1,nu1)
    !Forward Fourier transform
    call fourier(nn,1_ik,scratch)
    !Truncate
    call truncate(nn,scratch)
    !Calculate mean value
    call zero(nn,scratch)
    !$omp parallel do
    do l=nu2,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       scratch(i,j,k,l)=0.0_rk
    end do ; end do ; end do ; end do
    !$omp end parallel do
    tmp=mean_value(nn,scratch)

    mkhdis=2.0_rk*visc(nu1)*tmp(nu1)

    call zero(nn,scratch)
    call zero(nn,rmsarr)
    call zero(nn,rks1)
    
    return

  end function mean_kinetic_helicity_dissipation



!!$  subroutine maxima(nn,fu)
!!$    use mpivars
!!$    implicit none
!!$    integer(ik), dimension(1:4), intent(in)                   :: nn
!!$    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)  :: fu
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    !Calculate maxima of the scalar field
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$    integer(ik) :: i,j,k
!!$    
!!$    
!!$
!!$    return
!!$
!!$  end subroutine maxima

  subroutine cross_product(nn,c,a,b,nf)
    use types
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT)      :: c
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)       :: a
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)       :: b
    integer(ik), intent(in)                                     :: nf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate cross (outer) product of two vector fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: i,j,k
    
    !$omp parallel do
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       c(i,j,k,nf)=a(i,j,k,nf+1)*b(i,j,k,nf+2)-a(i,j,k,nf+2)*b(i,j,k,nf+1)
       c(i,j,k,nf+1)=a(i,j,k,nf+2)*b(i,j,k,nf)-a(i,j,k,nf)*b(i,j,k,nf+2)
       c(i,j,k,nf+2)=a(i,j,k,nf)*b(i,j,k,nf+1)-a(i,j,k,nf+1)*b(i,j,k,nf)
    end do; end do ; end do
    !$omp end parallel do

    return

  end subroutine cross_product

  subroutine curl(nn,fw,fu,nf)
    use parameters, only: dim1
    use data, only: wv,isactive
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT) :: fw
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)  :: fu
    integer(ik), intent(in)                                      :: nf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate curl of a vector field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                             :: i,j,k
    complex(ck)                                             :: tmpx, tmpy, tmpz, tmp

    !$omp parallel do private(tmpx, tmpy, tmpz, tmp)
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
       if(isactive(i,j,k)) then
          tmpx = cmplx(fu(i,j,k,nf  ),fu(i+1,j,k,nf  ),ck)
          tmpy = cmplx(fu(i,j,k,nf+1),fu(i+1,j,k,nf+1),ck)
          tmpz = cmplx(fu(i,j,k,nf+2),fu(i+1,j,k,nf+2),ck)

          tmp=ii*wv(2_ik,i,j,k)*tmpz-ii*wv(3_ik,i,j,k)*tmpy
          fw(i  ,j,k,nf)=real(tmp,rk)
          fw(i+1,j,k,nf)=aimag(tmp)

          tmp=ii*wv(3_ik,i,j,k)*tmpx-ii*wv(1_ik,i,j,k)*tmpz
          fw(i  ,j,k,nf+1)=real(tmp,rk)
          fw(i+1,j,k,nf+1)=aimag(tmp)

          tmp=ii*wv(1_ik,i,j,k)*tmpy-ii*wv(2_ik,i,j,k)*tmpx
          fw(i  ,j,k,nf+2)=real(tmp,rk)
          fw(i+1,j,k,nf+2)=aimag(tmp)
       else
          fw(i,j,k,nf:nf+2)=0.0_rk
          fw(i+1,j,k,nf:nf+2)=0.0_rk
       end if
    end do; end do ; end do
    !$omp end parallel do


    return

  end subroutine curl

  subroutine gradient(nn,fg,fu,nfs,nfe)
    use parameters, only: nsclf,nscll
    use data, only: wv,isactive
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4),nu1:nu3), intent(OUT) :: fg
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)  :: fu
    integer(ik), optional                                :: nfs, nfe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate the gradient of a vector field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                             :: i, j, k, l,nfstart, nfend
    complex(ck) :: tmp

    nfstart=nsclf
    nfend=nscll
    if(present(nfs)) nfstart=nfs
    if(present(nfe)) nfend=nfe

    call zero(nn,fg)
    
    !$omp parallel do private(tmp)
    do l=nfstart,nfend ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
       if(isactive(i,j,k)) then
          tmp=cmplx(fu(i,j,k,l),fu(i+1,j,k,l),ck)

          fg(i,j,k,l,nu1)=real(ii*wv(1_ik,i,j,k)*tmp,rk)
          fg(i+1,j,k,l,nu1)=aimag(ii*wv(1_ik,i,j,k)*tmp)

          fg(i,j,k,l,nu2)=real(ii*wv(2_ik,i,j,k)*tmp,rk)
          fg(i+1,j,k,l,nu2)=aimag(ii*wv(2_ik,i,j,k)*tmp)

          fg(i,j,k,l,nu3)=real(ii*wv(3_ik,i,j,k)*tmp,rk)
          fg(i+1,j,k,l,nu3)=aimag(ii*wv(3_ik,i,j,k)*tmp)
       else
          fg(i,j,k,l,nu1:nu3)=0.0_rk
          fg(i+1,j,k,l,nu1:nu3)=0.0_rk
       end if
    end do; end do; end do ; end do
    !$omp end parallel do


    return
  end subroutine gradient

  subroutine vector_potential(nn,fa,fu,nf)
    use types
    use data, only: wv,scratch,isactive
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT)  :: fa
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)   :: fu
    integer(ik)                                 :: nf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate the vector potential of a vector field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                              :: i,j,k
    real(rk)                                                 :: ksq

    !Get curl of the magnetic field
    call curl(nn,scratch,fu,nf)

    !$omp parallel do private(ksq)
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       ksq=wv(1_ik,i,j,k)**2+wv(2_ik,i,j,k)**2+&
               &wv(3_ik,i,j,k)**2
       if(isactive(i,j,k)) then
          fa(i,j,k,nf)=-ksq*scratch(i,j,k,nf)
          fa(i,j,k,nf+1)=-ksq*scratch(i,j,k,nf+1)
          fa(i,j,k,nf+2)=-ksq*scratch(i,j,k,nf+2)
       else
          fa(i,j,k,nf:nf+2)=0.0_rk
       end if
    end do; end do ; end do
    !$omp end parallel do

    return

  end subroutine vector_potential


  function mean_magnetic_helicity(nn,fu) result(mmh)
    use types
    use data, only: rks1,nb1
    implicit none
    real(rk)                                                :: mmh
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)               :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate mean magnetic helicity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                              :: i,j,k,l
    real(rk), dimension(1:nn(4))                                   :: tmp

    !Get vector potential of the magnetic field
    call vector_potential(nn,rks1,fu,nb1)
    !Inverse Fourier transforms
    call fourier(nn,-1_ik,rks1)
    !Perform inner product
    call inner_product(nn,rks1,nu1,rks1,&
         &fu,nb1)
    !Forward Fourier transform
    call fourier(nn,1_ik,rks1,nu1,nu1)
    !Truncate
    call truncate(nn,rks1,nu1,nu1)
    !$omp parallel do
    do l=nu2,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       rks1(i,j,k,l)=0.0_rk
    end do ; end do ; end do ; end do
    !$omp end parallel do
    !Calculate mean value
    tmp=mean_value(nn,rks1)
    !FIX!!!
    mmh=0.5_rk*tmp(nu1)

    !Set auxiliary arrays back to zero
    call zero(nn,rks1)


    return

  end function mean_magnetic_helicity


  function mean_magnetic_helicity_dissipation(nn,fu) result(mmhdis)
    use types
    use data, only: u,rmsarr,rks1,nb1,nu1,nu2,nu3
    implicit none
    real(rk)                                                :: mmhdis
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)               :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate mean magnetic helicity dissipation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk), dimension(1:nn(4))                                    :: tmp
    integer(ik) :: i,j,k,l
    
    call copy(nn,u,fu)

    !Get vorticity
    call curl(nn,rmsarr,fu,nb1)

    !Inverse Fourier transforms
    call fourier(nn,-1_ik,u,nb1)
    call fourier(nn,-1_ik,rmsarr,nb1)

    !Perform inner product
    call inner_product(nn,rks1,nu1,u,rmsarr,nb1)

    !Forward Fourier transform
    call fourier(nn,1_ik,rks1,nu1,nu1)

    !Truncate
    call truncate(nn,rks1)

    !$omp parallel do
    do l=nu2,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       rks1(i,j,k,l)=0.0_rk
    end do ; end do ; end do ;  end do
    !$omp end parallel do
    !Calculate mean value
    tmp=mean_value(nn,rks1)

    mmhdis = visc(nb1)*tmp(nu1)

    !Set auxiliary arrays back to zero
    call zero(nn,rmsarr)
    call zero(nn,rks1)

    return

  end function mean_magnetic_helicity_dissipation



  function ambipolar_diffusion_dissipation(nn,fu) result(addis)
    use types
    use mpivars
    use data, only: rmsarr,u,scratch,nb1,nb2,nb3
    implicit none
    real(rk)                                                :: addis
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)               :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate ambipolar diffusion dissipation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk)                                                  :: fac
    integer(ik) :: i,j,k
    
    fac=real(nn(1)*gn2*gn3,rk)**(-1)

    !Get magnetic field
    call copy(nn,u,fu)
    !Get current density
    call curl(nn,rmsarr,fu,nb1)

    !Inverse Fourier transforms
    call fourier(nn,-1_ik,u,nb1)
    call fourier(nn,-1_ik,rmsarr,nb1)

    call cross_product(nn,scratch,u,rmsarr,nb1)
    !$omp parallel do reduction(+:addis)
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       addis=addis+fac*(scratch(i,j,k,nb1)**2+scratch(i,j,k,nb2)**2+&
            &scratch(i,j,k,nb3)**2)
    end do; end do ; end do
    !$omp end parallel do

    !Reduce
#ifdef _MPI_
    sbuf(1)=addis
    call mpi_allreduce(sbuf,rbuf,1,MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
    addis=rbuf(1)
#endif

    addis=AD_COEFF*addis

    !Set auxiliary arrays back to zero
    call zero(nn,rmsarr)
    call zero(nn,scratch)

    return

  end function ambipolar_diffusion_dissipation




  function mean_cross_helicity_dissipation(nn,fu) result(mchdis)
    use types
    use data, only: rmsarr,rks1,nu1,nu3,nb1,nb3
    implicit none
    real(rk)                                  :: mchdis
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)              :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate mean cross helicity dissipation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk), dimension(1:nn(4)) :: tmp
    integer(ik) :: i,j,k,l

    !Get vorticity
    call curl(nn,rmsarr,fu,nu1)
    !Get current density
    call curl(nn,rks1,fu,nb1)

    !Inverse Fourier transforms
    call fourier(nn,-1_ik,rmsarr,nu1,nu3)
    call fourier(nn,-1_ik,rks1,nb1,nb3)

    !$omp parallel do
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       rmsarr(i,j,k,nb1:nb3)=rmsarr(i,j,k,nu1:nu3)
    end do; end do ; end do
    !$omp end parallel do

    !Perform inner product
    call inner_product(nn,rks1,nu1,rmsarr,rks1,nb1)

    !Forward Fourier transform
    call fourier(nn,1_ik,rks1)

    !Truncate
    call truncate(nn,rks1)

    !$omp parallel do
    do l=nu2,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       rks1(i,j,k,l)=0.0_rk
    end do; end do ; end do ;  end do
    !$omp end parallel do
    !Calculate mean value
    tmp=mean_value(nn,rks1)

    mchdis = (visc(nu1) + visc(nb1))*tmp(1)

    !Set auxiliary arrays back to zero
    call zero(nn,rmsarr)
    call zero(nn,rks1)

    return

  end function mean_cross_helicity_dissipation

  function mean_cross_helicity(nn,fu) result(mcross_hel)
    use types
    use data, only: rks1,rmsarr,nu1,nu2,nu3,nb1,nb3
    implicit none
    real(rk)                         :: mcross_hel
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)  :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate mean cross helicity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk), dimension(1:nn(4)) :: tmp 
    integer(ik) :: i,j,k,l

    !Copy to auxiliary arrays
    !$omp parallel do
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       rmsarr(i,j,k,nb1:nb3)=fu(i,j,k,nu1:nu3)
       rks1(i,j,k,nb1:nb3)=fu(i,j,k,nb1:nb3)
    end do; end do ; end do
    !$omp end parallel do


    !Inverse Fourier transforms
    call fourier(nn,-1_ik,rmsarr,nb1)
    call fourier(nn,-1_ik,rks1,nb1)

    mcross_hel=0.0_rk

    !Perform inner product
    call inner_product(nn,rks1,nu1,rmsarr,rks1,nb1)
    !Forward Fourier transform
    call fourier(nn,1_ik,rks1,nu1,nu1)
    !Truncate
    call truncate(nn,rks1,nu1,nu1)
    !$omp parallel do
    do l=nu2,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       rks1(i,j,k,l)=0.0_rk
    end do ; end do ; end do ;  end do
    !$omp end parallel do
    !Calculate mean value
    tmp=mean_value(nn,rks1)
    mcross_hel=0.5_rk*tmp(nu1)

    !Set auxiliary arrays back to zero
    call zero(nn,rmsarr)
    call zero(nn,rks1)

    return

  end function mean_cross_helicity

  function mean_kinetic_helicity(nn,fu) result(mkin_hel)
    use types
    use data, only: u,rmsarr,rks1,nu1
    implicit none
    real(rk)                                   :: mkin_hel
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)  :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate mean kinetic helicity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk), dimension(:), allocatable :: tmp
    integer(ik) :: i,j,k,l
    
    if(.not.allocated(tmp)) allocate(tmp(1:nn(4)))
    !Get velocity field
    call copy(nn,u,fu)

    !Get vorticity
    call curl(nn,rks1,fu,nu1)

    !Inverse Fourier transforms
    call fourier(nn,-1_ik,u,nu1,nu1)
    call fourier(nn,-1_ik,rks1,nu1,nu1)

    !Perform inner product
    call inner_product(nn,rmsarr,nu1,u,rks1,nu1)

    !Forward Fourier transform
    call fourier(nn,1_ik,rmsarr,nu1,nu1)

    !Truncate
    call truncate(nn,rmsarr,nu1,nu1)

    !$omp parallel do
    do l=nu2,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       rmsarr(i,j,k,l)=0.0_rk
    end do ; end do ; end do ; end do
    !$omp end parallel do
    !Calculate mean value
    tmp=mean_value(nn,rmsarr)
    mkin_hel=0.5_rk*tmp(nu1)

    !Set auxiliary arrays backto zero
    call zero(nn,rmsarr)
    call zero(nn,rks1)

    return

  end function mean_kinetic_helicity

  subroutine inner_product(nn,cross_hel,nfout,u,b,nfin)
    use types
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT) :: cross_hel
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)  :: u
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)  :: b
    integer(ik)                                :: nfout, nfin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate the inner product of two vector fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: i,j,k
    
    !$omp parallel do
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
    cross_hel(i,j,k,nfout) = u(i,j,k,nfin)*b(i,j,k,nfin)+u(i,j,k,nfin+1)*b(i,j,k,nfin+1)+&
         &u(i,j,k,nfin+2)*b(i,j,k,nfin+2)
    end do ; end do ; end do
    !$omp end parallel do

    return

  end subroutine inner_product

  subroutine dissipation(nn,e,nfout,fu)
    use data, only: wv,u,rks1, rmsarr,isactive,nu1,nu2,nu3
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT)              :: e
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)               :: fu
    integer(ik), intent(IN)                                 :: nfout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate total energy dissipation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                             :: i, j, k
    real(rk)                                                :: vtmp
    complex(ck)                                             :: tmpx, tmpy, tmpz, tmp
    !Do nothing if we're solving the Burgers equation 
    if(BURGERS) return



    !Active modes only
    !$omp parallel do private(tmpx, tmpy, tmpz, tmp, vtmp)
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
       if(isactive(i,j,k)) then
          tmpx=cmplx(fu(i,j,k,nu1),fu(i+1,j,k,nu1),ck)
          tmpy=cmplx(fu(i,j,k,nu2),fu(i+1,j,k,nu2),ck)
          tmpz=cmplx(fu(i,j,k,nu3),fu(i+1,j,k,nu3),ck)
          !du1/dx1
          tmp=ii*wv(1_ik,i,j,k)*tmpx
          u(i,j,k,nu1)=real(tmp,rk)
          u(i+1,j,k,nu1)=aimag(tmp)
          !du1/dx2
          tmp=ii*wv(2_ik,i,j,k)*tmpx
          u(i,j,k,nu2)=real(tmp,rk)
          u(i+1,j,k,nu2)=aimag(tmp)
          !du1/dx3
          tmp=ii*wv(3_ik,i,j,k)*tmpx
          u(i,j,k,nu3)=real(tmp,rk)
          u(i+1,j,k,nu3)=aimag(tmp)

          !du2/dx1
          tmp=ii*wv(1_ik,i,j,k)*tmpy
          rks1(i,j,k,nu1)=real(tmp,rk)
          rks1(i+1,j,k,nu1)=aimag(tmp)

          !du2/dx2
          tmp=ii*wv(2_ik,i,j,k)*tmpy
          rks1(i,j,k,nu2)=real(tmp,rk)
          rks1(i+1,j,k,nu2)=aimag(tmp)

          !du2/dx3
          tmp=ii*wv(3_ik,i,j,k)*tmpy
          rks1(i,j,k,nu3)=real(tmp,rk)
          rks1(i+1,j,k,nu3)=aimag(tmp)

          !du3/dx1
          tmp=ii*wv(1_ik,i,j,k)*tmpz
          rmsarr(i,j,k,nu1)=real(tmp,rk)
          rmsarr(i+1,j,k,nu1)=aimag(tmp)

          !du3/dx2
          tmp=ii*wv(2_ik,i,j,k)*tmpz
          rmsarr(i,j,k,nu2)=real(tmp,rk)
          rmsarr(i+1,j,k,nu2)=aimag(tmp)

          !du3/dx3
          tmp=ii*wv(3_ik,i,j,k)*tmpz
          rmsarr(i,j,k,nu3)=real(tmp,rk)
          rmsarr(i+1,j,k,nu3)=aimag(tmp)
       else
          u(i,j,k,nu1:nu3)=0.0_rk
          u(i+1,j,k,nu1:nu3)=0.0_rk
          rks1(i,j,k,nu1:nu3)=0.0_rk
          rks1(i+1,j,k,nu1:nu3)=0.0_rk
          rmsarr(i,j,k,nu1:nu3)=0.0_rk
          rmsarr(i+1,j,k,nu1:nu3)=0.0_rk
       end if
    end do; end do ; end do
    !$omp end parallel do



    !Inverse Fourier transforms
    call fourier(nn,-1_ik,u,nu1,nu3)
    call fourier(nn,-1_ik,rks1,nu1,nu3)
    call fourier(nn,-1_ik,rmsarr,nu1,nu3)

    if(.not.VISCOUS) then
       vtmp=1.0
    else
       vtmp=visc(nu1)
    end if

    !Calculate energy dissipation
    !$omp parallel do
    do k=1,nn(3) ;  do j=1,nn(2) ; do i=1,nn(1)
       e(i,j,k,nfout)=0.5_rk*vtmp*(4*rks1(i,j,k,nu2)**2 + 2*(rks1(i,j,k,nu3) +&
            & rmsarr(i,j,k,nu2))**2 + 4*rmsarr(i,j,k,nu3)**2 + 4*u(i,j,k,nu1)**2 + &
            &2*(rks1(i,j,k,nu1) + u(i,j,k,nu2))**2 + 2*(rmsarr(i,j,k,nu1) + &
            &u(i,j,k,nu3))**2)
    end do; end do ; end do 
    !$omp end parallel do

    !Set auxiliary arrays back to zero
    call zero(nn,u)
    call zero(nn,rmsarr)
    call zero(nn,rks1)

    return

  end subroutine dissipation



  subroutine integral_length_scale(nn,fu,L,lambda,nespec)
    use data, only: wv,nu1,nu2,nu3,nb1,nb2,nb3
    use mpivars
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    integer(ik), intent(IN)                                 :: nespec
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)               :: fu
    real(rk), intent(OUT)                                   :: L,lambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate the integral length scale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    complex(ck)                                             :: tmp1,tmp2,tmp3
    complex(ck)                                             :: tmp4,tmp5,tmp6
    real(dp), dimension(nespec+1)                          :: tmpespec
    real(dp), dimension(nespec+1)                          :: tmpespec2
    integer(ik)                                             :: i,j,k,iii,pos
    real(rk)                                                :: nrg,kk
    real(rk)                                                :: scale,dk
    real(rk)                                                :: tot,toten
    real(rk), dimension(nespec+1,2)                        :: espec
    real(rk)                                                :: tmp,diss,area
    
    !Total energy
    toten=KE
    if(MHD) toten=KE+ME

    !Scale factor
    scale=sqrt(2.0_rk/(n1*n2*n3))

    tmpespec(:)=0.0_rk
    tot=0.0_rk

    !Wave-number step
    dk=sqrt(real((nn(1)/2)**2+(gn2/2)**2+(gn3/2)**2,rk))/(real(nespec-1,rk))

    do i=1,nespec
       espec(i,1)=(i-1)*dk
    end do

    dk=(espec(nespec,1)-espec(1,1))/real(nespec-1,rk)

    !Main loop
    !$omp parallel do    private(iii,nrg,kk,pos,tmp1,tmp2,tmp3,&
    !$omp tmp4,tmp5,tmp6) reduction(+:tmpespec) reduction(+:tot) 
    do i=1,dim1(nn(1))-1,2
       iii=(i-1)/2+1
       do j=1,nn(2)
          do k=1,nn(3)
             !Calculate the energy spectrum
             !Velocity field
             tmp1=scale*cmplx(fu(i,j,k,nu1),fu(i,j,k,nu1),ck)
             tmp2=scale*cmplx(fu(i,j,k,nu2),fu(i,j,k,nu2),ck)
             tmp3=scale*cmplx(fu(i,j,k,nu3),fu(i,j,k,nu3),ck)
             !Magnetic field
             if(MHD) then
                tmp4=scale*cmplx(fu(i,j,k,nb1),fu(i+1,j,k,nb1),ck)
                tmp5=scale*cmplx(fu(i,j,k,nb2),fu(i+1,j,k,nb2),ck)
                tmp6=scale*cmplx(fu(i,j,k,nb3),fu(i+1,j,k,nb3),ck)
             else
                tmp4=cmplx(0.0_rk,0.0_rk,ck)
                tmp5=cmplx(0.0_rk,0.0_rk,ck)
                tmp6=cmplx(0.0_rk,0.0_rk,ck)
             end if
             !Energy
             nrg=0.5_rk*(abs(tmp1+tmp4)**2+abs(tmp2+tmp5)**2+abs(tmp3+&
                  &tmp6)**2)
             !Wave-vector magnitude squared
             kk=sqrt(real(wv(1_ik,i,j,k)**2+wv(2_ik,i,j,k)**2+&
                  &wv(3_ik,i,j,k)**2,rk))
             !tmpespec(:) index
             pos=int(kk/dk,ik)+1
             !Print an error message if we're out of bounds
             if(pos > nespec) pos = nespec
             if(pos < 1) pos = 1
             tmpespec(pos)=tmpespec(pos)+nrg
             tot=tot+nrg
          end do
       end do
    end do
    !$omp end parallel do

    !Reduce the array
#ifdef _MPI_
    tmpespec(nespec+1)=tot
    call mpi_allreduce(tmpespec,tmpespec2,int(nespec+1,i4b),MPIRK,MPI_SUM,&
         &MPI_COMM_WORLD,mpierr)
    tot=tmpespec2(nespec+1)
    espec(1:nespec,2)=tmpespec2(1:nespec)
#else
    espec(1:nespec,2)=tmpespec(1:nespec)
#endif

    if(tot==0.0_rk) tot=1._rk
    !Rescale
    scale=toten/tot
    !$omp parallel do
    do i=1,nespec
       espec(i,2)=scale*espec(i,2)/dk
    end do
    !$omp end parallel do

    !Integrate
    area=0.0_rk
    !$omp parallel do reduction(+:area)
    do i=1,nespec-1
       area=area+(espec(i,2)+espec(i+1,2))*0.5_rk*dk
    end do
    !$omp end parallel do

    !Reduce
#ifdef _MPI_
    sbuf(1)=area
    call mpi_allreduce(sbuf,rbuf,1,MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
    area=rbuf(1)
#endif

    !Integrate to ge the integral length scale
    L=0.0_rk
    diss=0.0_rk
    !$omp parallel do reduction(+:L,diss)
    do i=1,nespec-1
       if(i==1) then
          tmp=0.0_rk
       else
          tmp=espec(i,2)/espec(i,1)
       end if
       L=L+(tmp+espec(i+1,2)/espec(i+1,1))*0.5_rk*dk
       diss=diss+(espec(i,1)**2*espec(i,2)+&
            &espec(i+1,1)**2*espec(i+1,2))*0.5_rk*dk
    end do
    !$omp end parallel do

    !Reduce
#ifdef _MPI_
    sbuf(1)=L
    sbuf(2)=diss
    call mpi_allreduce(sbuf,rbuf,2,MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
    L=rbuf(1)
    diss=rbuf(2)
#endif
    if(area==0._rk) area=1.0_rk
    L=2.0_rk*PI*(area**(-1))*L
    if(diss==0.0_rk) diss=1.0_rk
    lambda=2.0_rk*PI*sqrt(area/diss)
    return

  end subroutine integral_length_scale


  subroutine vector_spectrum(nn,fu,nf,fname,nespec)
    use data, only:wv
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)               :: fu
    integer(ik), intent(IN)                                    :: nf
    character(*)                                               :: fname
    integer(ik), intent(IN)                                    :: nespec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate the energy spectrum of a vector field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    complex(ck)                                                :: tmp1,tmp2,tmp3
    real(rk), dimension(nespec+1)                              :: tmpespec
    integer(ik)                                                :: i,j,k,iii,pos
    real(rk)                                                   :: nrg,kk,scale
    real(rk)                                                   :: dk
    real(rk)                                                   :: tot,area
    integer(ik)                                                :: espec_file
    real(rk), dimension(nespec+1,2)                            :: espec
    real(rk)                                                   :: scale_fourier

    espec(:,:)=0.0_rk

    !Scale factors
    scale=1.0
    tmpespec(:)=0.0_rk
    tot=0.0_rk

    dk=1.0_rk

    do i=1,nespec
       espec(i,1)=(i-1)*dk
    end do

    !Calculate the energy spectrum
    !$omp parallel do private(iii,nrg,kk,pos,tmp1,tmp2,tmp3) &
    !$omp reduction(+:tmpespec) reduction(+:tot) 
    do i=1,dim1(nn(1))-1,2
       iii=(i-1)/2+1
       do j=1,nn(2)
          do k=1,nn(3)
             tmp1=scale*cmplx(fu(i,j,k,nf),fu(i+1,j,k,nf),ck)
             tmp2=scale*cmplx(fu(i,j,k,nf+1),fu(i+1,j,k,nf+1),ck)
             tmp3=scale*cmplx(fu(i,j,k,nf+2),fu(i+1,j,k,nf+2),ck)
             !Energy
             nrg=(abs(tmp1)**2+abs(tmp2)**2+abs(tmp3)**2)
             !twice the weight when kx non zero
             if (wv(1_ik,i,j,k).ne.0) nrg=nrg*2.0_rk
             !Wave-vector magnitude squared
             kk=sqrt(real(wv(1_ik,i,j,k)**2+wv(2_ik,i,j,k)**2+&
                  &wv(3_ik,i,j,k)**2,rk))
             !tempespec(:) index
             pos=min(nespec,int(kk/dk,ik)+1)
             !Print an error message if we're out of bounds
             if(pos>nespec) print '(4i4)',  pos,i,j,k
             tmpespec(pos)=tmpespec(pos)+nrg
             tot=tot+nrg
          end do
       end do
    end do
    !$omp end parallel do
    espec(1:nespec,2)=tmpespec(1:nespec)

    !get area    
    area=0.0_rk
    !$omp parallel do reduction(+:area)
    do i=1,nespec
       area=area+espec(i,2)*dk
    end do
    !$omp end parallel do

    scale_fourier=1.0_rk/(real(n1,rk)*real(n2,rk)*real(gn3,rk))
    scale_fourier=scale_fourier**2
    !Rescale
    espec(:,2)=scale_fourier*espec(:,2)



    !Write to file
    open(newunit=espec_file,file=fname)

    do i=1,nespec
       write(espec_file,'(2e30.14)') espec(i,1),espec(i,2)
    end do
    close(espec_file)



  end subroutine vector_spectrum

  subroutine scalar_spectrum(nn,fu,nf,fname,npsspec)
    use data,only:wv
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)  :: fu
    integer(ik), intent(IN)                                    :: nf
    integer(ik), intent(IN)                                    :: npsspec
    character(*), intent(IN)                                   :: fname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate the energy spectrum of a scalar field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk), dimension(npsspec+1)                             :: tmppsspec
    complex(ck)                                                :: tmp1
    integer(ik)                                                :: i,j,k,iii,pos
    real(rk)                                                   :: nrg,kk,scale
    real(rk)                                                   :: dk,tot
    real(rk)                                                   :: area
    integer(ik)                                                :: espec_file
    real(rk), dimension(npsspec+1,2)                           :: psspec
    real(rk)                                                   :: scale_fourier

    !Scale factors
    scale=1.0
    scale_fourier=1.0_rk/(real(n1,rk)*real(n2,rk)*real(n3,rk))
    scale_fourier=scale_fourier**2

    tmppsspec(:)=0.0_rk
    tot=0.0_rk

    dk=1.0_rk

    do i=1,npsspec
       psspec(i,1)=(i-1)*dk
    end do

    !$omp parallel do    private(iii,nrg,kk,pos,tmp1) &
    !$omp reduction(+:tmppsspec) reduction(+:tot)
    do i=1,dim1(nn(1))-1,2
       iii=(i-1)/2+1
       do j=1,nn(2)
          do k=1,nn(3)
             tmp1=scale*cmplx(fu(i,j,k,nf),fu(i+1,j,k,nf),ck)
             !Energy
             nrg=abs(tmp1)**2
             !Twice the weight when kx non zero
             if (wv(1_ik,i,j,k).ne.0) nrg=nrg*2
             !Wave-vector magnitude squared
             kk=sqrt(real(wv(1_ik,i,j,k)**2+wv(2_ik,i,j,k)**2&
                  &+wv(3_ik,i,j,k)**2,rk))
             !tmpespec(:) index
             pos=min(npsspec,int(kk/dk,ik)+1)
             tmppsspec(pos)=tmppsspec(pos)+nrg
             tot=tot+nrg
          end do
       end do
    end do
    !$omp end parallel do



    psspec(:,2)=tmppsspec(:)

    !Get area
    psspec(:,2)=psspec(:,2)/dk
    area=0.0_rk
    !$omp parallel do reduction(+:area)
    do i=1,npsspec
       area=area+psspec(i,2)*dk
    end do
    !$omp end parallel do

    !Rescale
    psspec(:,2)=scale_fourier*psspec(:,2)

    !Write to file
    scale=1.0
    open(newunit=espec_file,file=fname)
    do i=1,npsspec
       write(espec_file,'(2e30.14)') psspec(i,1),psspec(i,2)
    end do

    close(espec_file)

    return

  end subroutine scalar_spectrum



  function incompressibility(nn,fu,nf) result(meandivv)
    use data, only: wv,scratch,nu1
    use mpivars
    implicit none
    real(rk)                                  :: meandivv
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)              :: fu
    integer(ik), intent(in) :: nf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Check the incompressibility condition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: i,j,k
    real(rk), dimension(:), allocatable                    :: tmp
    complex(ck)                                            :: ctmp,ctmp2,ctmp3,ctmp4

    if(.not.allocated(tmp)) allocate(tmp(1:nn(4)))
    tmp = 0.0_rk
    !$omp parallel do private(ctmp,ctmp2,ctmp3,ctmp4)
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
       ctmp2=cmplx(fu(i,j,k,nf),fu(i+1,j,k,nf),ck)
       ctmp3=cmplx(fu(i,j,k,nf+1),fu(i+1,j,k,nf+1),ck)
       ctmp4=cmplx(fu(i,j,k,nf+2),fu(i+1,j,k,nf+2),ck)
       ctmp=ii*(wv(1_ik,i,j,k)*ctmp2+wv(2_ik,i,j,k)*ctmp3+&
            &wv(3_ik,i,j,k)*ctmp4)
       scratch(i,j,k,nu1)=real(ctmp,rk)
       scratch(i+1,j,k,nu1)=aimag(ctmp)
       scratch(i,j,k,nu1+1:nn(4)) = 0.0_rk
       scratch(i+1,j,k,nu1+1:nn(4)) = 0.0_rk
    end do; end do ; end do
    !$omp end parallel do

    call msvalue(nn,scratch,tmp)
    meandivv=sqrt(tmp(nu1))

    call zero(nn,scratch)
    
    return

  end function incompressibility

  subroutine truncate(nn,fu,nfs,nfe)
    use data, only: isactive
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN OUT) :: fu
    integer(ik), optional :: nfs,nfe
!!!!!!!!!!!!!!!!!!!!!!!!
    !Truncate a scalar field
!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                            :: i,j,k
    integer(ik)                                            :: nfi,nnfs,nnfe

    nnfs=1
    nnfe=nn(4)
    if(present(nfs)) nnfs=nfs
    if(present(nfe)) nnfe=nfe
    !$omp parallel do
    do nfi=nnfs,nnfe ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       if(.not.isactive(i,j,k)) fu(i,j,k,nfi)=0.0_rk
    end do; end do ; end do ; end do
    !$omp end parallel do

    return

  end subroutine truncate


  subroutine rescale(nn,fu)
    use data, only: u, nu1,nu3,nsclf,nscll,nb1,nb3
    use mpivars
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT) :: fu

!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Rescale all flow variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                                :: i,j,k,l
    real(rk), dimension(:), allocatable ::tmp
    real(rk) :: smax, smin

    if(.not.allocated(tmp)) allocate(tmp(1:nn(4)))

    call msvalue(nn,fu,tmp)
    !$omp parallel do
    do l=nu1,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       fu(i,j,k,l)=fu(i,j,k,l)/sqrt(tmp(nu1))
    end do; end do ; end do ; end do
    !$omp end parallel do
    if(PASSIVE_SCALAR) then
       !$omp parallel do
       do l=nsclf,nscll ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
          fu(i,j,k,l)= fu(i,j,k,l)/sqrt(tmp(l))
       end do; end do ; end do ; end do
       !$omp end parallel do
       
       if(RADIATION) then
          call copy(nn,u,fu)
          call fourier(nn,-1_rk,u,nfs=ntemp,nfe=ntemp)
          smax = maxval(u(:, :, :, ntemp))
          smin = minval(u(:, :, :, ntemp))
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))   
             u(i,j,k,ntemp)= 2000.0 * (u(i,j,k,ntemp) - smin) / &
                  &(smax - smin) + 300.0
          end do; end do ; end do
          !$omp end parallel do
          call copy(nn,fu,u)
          call fourier(nn,1_rk,fu,nfs=ntemp,nfe=ntemp)
       end if


    end if


    if(MHD) then
       !$omp parallel do
       do l=nb1,nb3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
          fu(i,j,k,l)=fu(i,j,k,l)/sqrt(tmp(nb1))
       end do; end do ; end do ; end do
       !$omp end parallel do
    end if

    return

  end subroutine rescale

  subroutine shift(nn,dir,fu)
    use data, only: wv,isactive,phases,iphases
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    integer(ik), intent(IN)                                    :: dir
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN OUT) :: fu
!!!!!!!!!!!!!!!!!!
    !Phase-space shift
!!!!!!!!!!!!!!!!!!
    integer(ik) :: i,j,k,l,n
    complex(ck) :: tmp, phase

    n=max(nn(1),nn(2))

    if(STORE_PHASES) then
       !$omp parallel do private(tmp,phase)
       do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ;  do i=1,dim1(nn(1))-1,2
          if(isactive(i,j,k)) then
             tmp=cmplx(fu(i,j,k,l),fu(i+1,j,k,l),ck)
             if(dir==1) then
                phase=phases(i,j,k)
             else
                phase=iphases(i,j,k)
             end if
             tmp=phase*tmp
             fu(i,j,k,l)=real(tmp,rk)
             fu(i+1,j,k,l)=aimag(tmp)
          else
             fu(i,j,k,l)=0.0_rk
             fu(i+1,j,k,l)=0.0_rk
          end if
       end do;  end do ; end do ; end do
       !$omp end parallel do
    else
       !$omp parallel do private(tmp, phase)
       do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ;  do i=1,dim1(nn(1))-1,2
          if(isactive(i,j,k)) then
             tmp=cmplx(fu(i,j,k,l),fu(i+1,j,k,l),ck)
             phase=exp(dir*ii*(PI/real(n,rk))*(wv(1_ik,i,j,k)+&
                  &wv(2_ik,i,j,k)+wv(3_ik,i,j,k)))
             tmp=phase*tmp
             fu(i,j,k,l)=real(tmp,rk)
             fu(i+1,j,k,l)=aimag(tmp)
          else
             fu(i,j,k,l)=0.0_rk
             fu(i+1,j,k,l)=0.0_rk
          end if
       end do; end do ; end do ; end do
       !$omp end parallel do
    end if

    return

  end subroutine shift

  subroutine fourier(nn,dir,u,nfs,nfe,trunc)
    use mpivars
#ifdef _MPI_
    use fft_heffte
    use mpivars
#endif
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    integer(ik), intent(IN)                                    :: dir
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN OUT)              :: u
    integer(ik), optional :: nfs, nfe
    logical, optional :: trunc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Fast Fourier transform driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical :: tr
    integer(ik) :: nnfs, nnfe

    nnfs=1
    nnfe=nn(4)
    if(present(nfs)) nnfs=nfs
    if(present(nfe)) nnfe=nfe
    tr=.true.
    if(present(trunc)) tr=trunc


    if(dir==-1.and.tr) call truncate(nn,u,nnfs,nnfe)

#ifdef _MPI_
    call fft_heffte_fourier(nn,dir,u,nnfs,nnfe)
#endif

    return

  end subroutine fourier

  subroutine cfl_condition(nn,dt)
    use types
    use data, only: u,fu,rmsarr,scratch,nu1,nu3,nb1,nb3
    use mpivars
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    real(rk), intent(OUT) :: dt
!!!!!!!!!!!!!!!!!!!!
    !Check CFL condition
!!!!!!!!!!!!!!!!!!!!
    integer(ik)            :: i,j,k
    real(rk)               :: dx,maxb2,dt_ad,mvel,vmag
    real(rk), dimension(1:3) :: v_total


    ! set auxiliary arrays to zero
    call zero(nn,rmsarr)
    call zero(nn,scratch)

    ! get velocity field in physical space
    call copy(nn,u,fu)
    call fourier(nn,-1_rk,u,nfs=nu1,nfe=nu3)

    if(MHD) then
       ! get magnetic field in physical space
       call fourier(nn,-1_rk,u,nfs=nb1,nfe=nb3)
       if(AMB_DIFF) then
          ! ambipolar diffusion
          ! rmsarr = j = nabla x b
          call curl(nn,rmsarr,fu,nb1)
          call fourier(nn,-1_ik,rmsarr,nb1)
          ! scratch = j x b
          call cross_product(nn,scratch,rmsarr,u,nb1)
       end if
       if(HALL) then
          !Hall term
          ! rmsarr = = j = nabla x b
          call curl(nn,rmsarr,fu,nb1)
          call fourier(nn,-1_ik,rmsarr,nb1)
       end if
    end if

    mvel=0.0_rk
    vmag=0.0_rk
    v_total(1:3)=0.0_rk
    !$omp parallel do private(v_total,vmag) reduction(max:mvel)
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       v_total(1:3)=u(i,j,k,nu1:nu3)
       if(MHD) then
          v_total(1:3)=v_total(1:3)+u(i,j,k,nb1:nb3)
          if(AMB_DIFF) v_total(1:3)=v_total(1:3)+AD_COEFF*rmsarr(i,j,k,nu1:nu3)
          if(HALL) v_total(1:3)=v_total(1:3)+HALL_COEFF*rmsarr(i,j,k,nb1:nb3)
       end if
       vmag=sqrt(v_total(1)**2+v_total(2)**2+v_total(3)**2)
       mvel=max(mvel,vmag)
    end do; end do ; end do 
    !$omp end parallel do

    maxb2=1.0

    !Reduce
#ifdef _MPI_
    sbuf(1)=mvel
    call mpi_allreduce(sbuf,rbuf,1,MPIRK,MPI_MAX,MPI_COMM_WORLD,mpierr)
    mvel=rbuf(1)
#endif

        
    !Set-up CFL condition
    dx=LBOX/max(nn(1),gn2,gn3)
    dt=CFL*dx/mvel
    if(MHD) then
       if(AMB_DIFF) then
          dt_ad=CFL*dx**2/(AD_COEFF*maxb2)
          dt=min(dt,dt_ad)
       end if
    end if

    !Set auxiliary arrays back to zero
    call zero(nn,rmsarr)

    return

  end subroutine cfl_condition

end module numerics
