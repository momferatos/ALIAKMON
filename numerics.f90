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
  implicit none
  !
  ! physics & numerics module
  !
contains

  subroutine interpolate_qr(nn, qr, qr_out, x)
    use mpivars, only: ljstart, ljsize, lkstart, lksize, mpirank
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    real(rks), dimension(1:nn(1),1:nn(2),1:nn(3),1:3), intent(in) :: qr
    real(rk), dimension(3), intent(out) :: qr_out
    real(rk), dimension(3), intent(in) :: x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: i0, j0, k0
    integer(ik) :: i1, j1, k1
    integer(ik) :: i2, j2, k2
    integer(ik) :: jstart, jend, kstart, kend
    real(rk) :: dx, dy, dz
    real(rk) :: xd, yd, zd
    real(rk), dimension(3) :: c000, c100, c010, c110, c001, c101, c011, c111
    real(rk), dimension(3) :: c00, c01, c10, c11
    real(rk), dimension(3) :: c0, c1
    real(rk), dimension(3) :: c

    qr_out=0.0
    
    dx=LBOX/real(nn(1)-1_ik,rk)
    dy=LBOX/real(gn2-1_ik,rk)
    dz=LBOX/real(gn3-1_ik,rk)

    i0=int(floor(x(1)/dx),ik)
    j0=int(floor(x(2)/dy),ik)
    k0=int(floor(x(3)/dz),ik)

    i1=i0+1
    j1=j0+1
    k1=k0+1

    i2=i0+2
    j2=j0+2
    k2=k0+2

    
    if(i1<1.or.i2>nn(1)) then
       qr_out=0.0_rk
       return
    end if

    jstart=ljstart
    jend=ljstart+ljsize-1
    if(j1<jstart.or.j2>jend) then
       qr_out=0.0_rk
       return
    else
       j1=j1-ljstart+1
       j2=j2-ljstart+1
    end if

    kstart=lkstart
    kend=lkstart+lksize-1
    if(k1<kstart.or.k2>kend) then
       qr_out=0.0_rk
       return
    else
       k1=k1-lkstart+1
       k2=k2-lkstart+1
    end if

    xd=modulo(x(1), dx)/dx
    yd=modulo(x(2), dy)/dy
    zd=modulo(x(3), dz)/dz 

   
    c000=qr(i1,j1,k1,1:3)
    c100=qr(i2,j1,k1,1:3)
    c010=qr(i1,j2,k1,1:3)
    c110=qr(i2,j2,k1,1:3)
    c001=qr(i1,j1,k2,1:3)
    c101=qr(i2,j1,k2,1:3)
    c011=qr(i1,j2,k2,1:3)
    c111=qr(i2,j2,k2,1:3)

    c00=c000*(1-xd)+c100*xd
    c01=c001*(1-xd)+c101*xd
    c10=c010*(1-xd)+c110*xd
    c11=c011*(1-xd)+c111*xd

    c0=c00*(1-yd)+c10*yd
    c1=c01*(1-yd)+c11*yd

    c=c0*(1-zd)+c1*zd

    qr_out=c

    return
  end subroutine interpolate_qr

  subroutine integrate_qr(nn,qr,time)
    use data, only: nsects, s, ga
    use mpi
    use mpivars
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    real(rks), dimension(1:nn(1),1:nn(2),1:nn(3),1:3), intent(in) :: qr
    real(rk), intent(in) :: time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk) :: radius, halfbox, quarterbox, qr_proj, red, vol
    real(rk), dimension(3) :: center, point, qr_vec, shat
    integer(ik) :: ns
    real(rk) :: eps=1.0e-3
    real(rk) :: val
    real(rk) :: dx, dy, dz
    real(rk) :: tmp_send(2), tmp_recv(2)

    !$acc update self(qr(1:nn(1),1:nn(2),1:nn(3),1:3),&
    !$acc& ga(1:nn(1),1:nn(2),1:nn(3)))
    
    halfbox=LBOX/2.0_rk
    quarterbox=LBOX/4.0_rk
    radius=quarterbox
    center=halfbox
    vol=LBOX**3

    dx=LBOX/real(nn(1),rk)
    dy=LBOX/real(gn2,rk)
    dz=LBOX/real(gn3,rk)

    val=0.0_rk
    do ns=1,nsects
       shat=s(ns,1:3)/dot_product(s(ns,1:3),s(ns,1:3))
       point(1:3)=center(1:3)+radius*s(ns,1:3)
       call interpolate_qr(nn,qr,qr_vec,point)
       qr_proj=dot_product(qr_vec,s(ns,1:3))
       val=val+qr_proj
    end do

    red=(sum(ga)*dx*dy*dz)/CLIGHT
    tmp_send(1)=val; tmp_send(2)=red
    call mpi_reduce(tmp_send,tmp_recv,2_i4b,MPIRK,&
         &MPI_SUM,MPIROOT,MPI_COMM_WORLD,mpierr)
    val=tmp_recv(1)
    red=tmp_recv(2)
    if(mpirank==MPIROOT) then
       write(432,*) time, val*radius**2, red
       flush(432)
    end if
        
    return
    
  end subroutine integrate_qr
  

  subroutine shift(nn,dir,fu,nfs,nfe,idx)
    use data, only: isactive, phases, iphases
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    integer(ik), intent(IN)                                    :: dir
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(INOUT) :: fu
    integer(ik), optional :: nfs, nfe
    integer(ik), optional :: idx
!!!!!!!!!!!!!!!!!!
    !Phase-space shift
!!!!!!!!!!!!!!!!!!
    integer(ik) :: i,j,k,l,n
    complex(ck) :: tmp
    integer(ik) :: nnfs, nnfe
    integer(ik) :: iidx
    real(rk)                                     :: hdx
    real(rk), dimension(3)                       :: ddx
    complex(ck)                                  :: phase, iphase
    nnfs=1
    if(present(nfs)) nnfs=nfs

    nnfe=nn(4)
    if(present(nfe)) nnfe=nfe

    iidx=1
    if(present(idx)) iidx=idx

    n=max(nn(1),nn(2))

    hdx = PI / real(n, rk)

    !$omp parallel do private(tmp,phase)
    do l=nnfs,nnfe ; do k=1,nn(3) ; do j=1,nn(2) ;  do i=1,dim1(nn(1))-1,2
       if(isactive(i,j,k)) then
          tmp=cmplx(fu(i,j,k,l),fu(i+1,j,k,l),ck)
          if(dir==1) then
             phase=phases(i,j,k,iidx)
          else
             phase=iphases(i,j,k,iidx)
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


    return

  end subroutine shift
  
  subroutine msvalue(nn,fu,msv)
    use mpivars
    use data, only: rmsarr,nu1,nu2,nu3,nb1,nb2,nb3, copy, zero
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
    call mpi_allreduce(MPI_IN_PLACE,msv,int(nn(4)),MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

    tmp=msv(nu1)+msv(nu2)+msv(nu3)
    msv(nu1:nu3)=tmp
    

    if(MHD) then
       tmp=msv(nb1)+msv(nb2)+msv(nb3)
       msv(nb1:nb3)=tmp
    end if

    call zero(nn,rmsarr)

    return

  end subroutine msvalue

  function mean_value(nn,fu) result(meanval)
    use mpivars
    use data, only: rmsarr,nu1,nu3,nsclf,nscll,nb1,nb3,copy,zero
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rk), dimension(1:nn(4))                            :: meanval
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)                   :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !calculate mean-square value for vector array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                                 :: i,j,k,l
    real(rk) :: fmean(nn(4))
    real(rk)                                                    :: fac,tmp
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
    call mpi_allreduce(MPI_IN_PLACE,meanval,int(nn(4)),MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

    return

  end function mean_value

  subroutine check_free_slip_bcs(nn,fu)
    use data, only: u, copy, zero, fsclgrads, fsclgrad
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(inout)               :: fu
    integer(ik) :: j0, jpi
    real(rk) :: v0, vpi, d0, dpi

    call copy(nn,u,fu)
    call fourier(nn,-1_ik,u)

    j0=1
    jpi=nn(2)/2+1
    
    v0=maxval(abs(u(1:nn(1),j0,1:nn(3),nu2)))
    vpi=maxval(abs(u(1:nn(1),jpi,1:nn(3),nu2)))
    print '(a,2(a,e10.3))', 'BCs for v: ', ' y=0: ', v0, '| y=PI: ', vpi

    call gradient(nn,fsclgrads,fu,nu1,ntemp)
    fsclgrad=>fsclgrads(:,:,:,:,2)
    call fourier(nn,-1_ik,fsclgrad)

    d0=maxval(abs(fsclgrad(1:nn(1),j0,1:nn(3),nu1)))
    dpi=maxval(abs(fsclgrad(1:nn(1),jpi,1:nn(3),nu1)))
    print '(a,2(a,e10.3))', 'BCs for du/dy: ', ' y=0: ', d0, '| y=PI: ', dpi

    d0=maxval(abs(fsclgrad(1:nn(1),j0,1:nn(3),nu3)))
    dpi=maxval(abs(fsclgrad(1:nn(1),jpi,1:nn(3),nu3)))
    print '(a,2(a,e10.3))', 'BCs for dw/dy: ', ' y=0: ', d0, '| y=PI: ', dpi
    
    d0=maxval(abs(fsclgrad(1:nn(1),j0,1:nn(3),ntemp)))
    dpi=maxval(abs(fsclgrad(1:nn(1),jpi,1:nn(3),ntemp)))
    print '(a,2(a,e10.3))', 'BCs for dT/dy: ', ' y=0: ', d0, '| y=PI: ', dpi
    call zero(nn,fsclgrad);

    return

  end subroutine check_free_slip_bcs
  
  subroutine apply_free_slip_bcs(nn, fu)
    use data, only: scratch,copy,tr_wv_idx
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(inout)               :: fu
    integer(ik) :: i,j,k,l,jplus,jminus
    real(rk) ::  tmp1, tmp2

    call copy(nn,scratch,fu)

    do l=1,nn(4) ; do k=1,nn(3) ; do j=0,nn(2)/2 ; do i=1,dim1(nn(1))
       jplus=tr_wv_idx(j,nn(2))
       jminus=tr_wv_idx(-j,nn(2))
       if(l==nu2) then
          tmp1 = 0.5_rk*(scratch(i,jplus,k,l)-scratch(i,jminus,k,l))
          tmp2 = 0.5_rk*(scratch(i,jminus,k,l)-scratch(i,jplus,k,l))
          fu(i,jplus,k,l) = tmp1
          fu(i,jminus,k,l) = tmp2
       else
          tmp1 = 0.5_rk*(scratch(i,jplus,k,l)+scratch(i,jminus,k,l))
          tmp2 = 0.5_rk*(scratch(i,jminus,k,l)+scratch(i,jplus,k,l))
          fu(i,jplus,k,l) = tmp1
          fu(i,jminus,k,l) = tmp2
       end if
    end do; end do ; end do ; end do

    return

  end subroutine apply_free_slip_bcs
  
  subroutine project(nn,fu,press)
    use mpi
    use mpivars
    use data, only: wv,nu1,nb1,isactive,scratch,zero
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(INOUT)               :: fu
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3)), intent(OUT)                  :: press
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
                   if(RADIATION) press(i,j,k)=DENS*p
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

    if(RADIATION) then
       scratch(:,:,:,1)=press(:,:,:)
       call fourier(nn,-1_ik,scratch,nfs=1_ik,nfe=1_ik)
       press(:,:,:)=scratch(:,:,:,1)
       call zero(nn,scratch)
    end if
    
    return

  end subroutine project



  subroutine timestep(nn,fu,dt)
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(INOUT)  :: fu
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

    if(BOUNDCOND==FREESLIP) call apply_free_slip_bcs(nn,fu)
    
    return

  end subroutine timestep



  subroutine runge_kutta2(nn,fu,dt)
    use data, only: wv,rhs,rks1,isactive
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(INOUT)               :: fu
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
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(INOUT)               :: fu
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
    use data, only: isactive,u,du,psu,scratch,fnls,nu1,nu2,nu3,nb1,nb2,nb3,&
         &nsclf,nscll,ad,&
         &fsclgrads,fsclgrad,wv,arr_en_1,fdivqr,copy,zero,scratch2,press
    use mpivars
    use fvdom, only: calcia
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT) :: fnl
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN) :: fu 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Non-linear terms in rotational form
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                             :: i, j, k,l,m
    real(rk) :: ksq

    call compute_hydro

    if(PASSIVE_SCALAR) call compute_passive_scalar

    if(BOUSSINESQ) call compute_boussinesq
    
    if(MHD) then
       call compute_mhd
       !Get the Lorentz force
       call zero(nn,scratch2)
       call lorentz_force(nn,scratch2,ad,fu)
    end if
    

    if(RADIATION.and.RADIATION_COUPLING) call compute_radiation
    
    ! handle Patterson-Orszag deliasing
    if(DEALIASING==PATTERSON_ORSZAG) then
       call shift(nn,-1_ik,fnls)
       !$omp parallel do
       do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
          if(isactive(i,j,k)) then
             fnl(i,j,k,l)=0.5_rk*(fnl(i,j,k,l)+fnls(i,j,k,l))+&
                  &scratch2(i,j,k,l)
          end if
       end do; end do ; end do ; end do
       !$omp end parallel do
    end if

    if(RADIATION.and.RADIATION_COUPLING) then
       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
          if(isactive(i,j,k)) then
             fnl(i,j,k,ntemp)=fnl(i,j,k,ntemp)+(CP)**(-1)*fdivqr(i,j,k)
          end if
       end do ; end do ; end do
       !$omp end parallel do
    end if
    
    ! truncate non-linear terms to remove aliasing errors
    call truncate(nn,fnl)

    
    
    
    
    ! compute diffusive terms
    call compute_diffusive_terms

    ! project to zero divergence
    if(.not.BURGERS) call project(nn,fnl,press)
    
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
         fnl(i,j,k,nu1)= u(i,j,k,nu2)*scratch(i,j,k,nu3) - &
              &u(i,j,k,nu3)*scratch(i,j,k,nu2)
         fnl(i,j,k,nu2)= u(i,j,k,nu3)*scratch(i,j,k,nu1) - &
              &u(i,j,k,nu1)*scratch(i,j,k,nu3)
         fnl(i,j,k,nu3)= u(i,j,k,nu1)*scratch(i,j,k,nu2) - &
              &u(i,j,k,nu2)*scratch(i,j,k,nu1)
      end do; end do ; end do   
      !$omp end parallel do

#ifdef _MPI_
      call mpi_allreduce(MPI_IN_PLACE,MAXVEL,1,MPIRK,MPI_MAX,MPI_COMM_WORLD,mpierr)
      call mpi_allreduce(MPI_IN_PLACE,MAXVORT,1,MPIRK,MPI_MAX,MPI_COMM_WORLD,mpierr)
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

      !Set passive scalar components
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
         fnl(i,j,k,l)=u(i,j,k,nu1)*fsclgrads(i,j,k,l,nu1)+&
              &u(i,j,k,nu2)*fsclgrads(i,j,k,l,nu2)+&
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
            fnls(i,j,k,l)=psu(i,j,k,nu1)*fsclgrads(i,j,k,l,nu1)+&
                 &psu(i,j,k,nu2)*fsclgrads(i,j,k,l,nu2)+&
                 &psu(i,j,k,nu3)*fsclgrads(i,j,k,l,nu3)
         end do; end do ; end do ; end do
         !$omp end parallel do
         call fourier(nn,1_ik,fnls,nsclf,nscll)
      end if



      return

    end subroutine compute_passive_scalar

    subroutine compute_boussinesq
      real(rk) :: ksq
      !$omp parallel do private(ksq)
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         if(isactive(i,j,k)) then
            ksq=wv(1_ik,i,j,k)**2+wv(2_ik,i,j,k)**2+wv(3_ik,i,j,k)**2
            if(ksq/=0.0) then
               fnl(i,j,k,nu2)=fnl(i,j,k,nu2)-AEXP*GGRAV*u(i,j,k,ntemp)
            else
               fnl(i,j,k,nu2)=fnl(i,j,k,nu2)+AEXP*GGRAV*0.5_rk*TEMP0
            end if
         end if
      end do; end do ; end do
      !$omp end parallel do
      
    end subroutine compute_boussinesq
    
    subroutine compute_radiation
      use data, only: qr, fqr, fdivqr, fdivqr_tmp, ia, iba, copy, temp,&
           &scratch2
      implicit none
      real(rk) :: scale, tmp, tt, cv, dens, vappress, pp, dens_air, y
      integer(ik), dimension(4) :: nnn
      integer :: icode
      integer(ik) :: nidx
      real(8), external :: cvtp
      real(8), external :: cptp
      real(8), external :: psatt
      real(8), external :: dtp

      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         fdivqr(i,j,k) = 0.0_rk
         fdivqr_tmp(i,j,k) = 0.0_rk
      end do; end do ;  end do
      !$omp end parallel do
      call zero(nn,scratch2)
      call copy(nn,scratch2,fu,nfs=ntemp,nfe=ntemp)
      call fourier(nn,-1_ik,scratch2,nfs=ntemp,nfe=ntemp)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         temp(i,j,k)=scratch2(i,j,k,ntemp)
      end do; end do ;  end do
      !$omp end parallel do

      call calcia(nn, temp)

      nnn(1:3) = nn(1:3)
      nnn(4) = 3_ik
      call copy(nnn, fqr, qr, nfs=1_ik, nfe=3_ik)
      call fourier(nnn, 1_ik, fqr, nfs=1_ik, nfe=3_ik)
      call divergence(nnn, fdivqr, fqr)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nidx=1
      call zero(nn,scratch2)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         fdivqr_tmp(i,j,k) = 0.0_rk
      end do; end do ;  end do
      !$omp end parallel do

      call zero(nn,scratch2)
      call copy(nn,scratch2,fu, nfs=ntemp, nfe=ntemp)
      call shift(nn,1_ik, scratch2, nfs=ntemp, nfe=ntemp, idx=nidx)
      call fourier(nn,-1_ik,scratch2, nfs=ntemp, nfe=ntemp)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         temp(i,j,k)=scratch2(i,j,k,ntemp)
      end do; end do ;  end do
      !$omp end parallel do

      call calcia(nn, temp)

      nnn(1:3) = nn(1:3)
      nnn(4) = 3_ik
      call copy(nnn, fqr, qr, nfs=1_ik, nfe=3_ik)
      call fourier(nnn, 1_ik, fqr, nfs=1_ik, nfe=3_ik)
      call divergence(nnn, fdivqr_tmp, fqr)

      call zero(nn,scratch2)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         scratch2(i,j,k,ntemp) = fdivqr_tmp(i,j,k)
      end do; end do ;  end do
      !$omp end parallel do

      call shift(nn,-1_ik,scratch2,nfs=ntemp, nfe=ntemp,idx=nidx)

      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         fdivqr(i,j,k) = fdivqr(i,j,k) + scratch2(i,j,k,ntemp)
      end do; end do ;  end do
      !$omp end parallel do

      call zero(nn,scratch2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nidx=1
      call zero(nn,scratch2)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         fdivqr_tmp(i,j,k) = 0.0_rk
      end do; end do ;  end do
      !$omp end parallel do

      call zero(nn,scratch2)
      call copy(nn,scratch2,fu, nfs=ntemp, nfe=ntemp)
      call shift(nn,1_ik, scratch2, nfs=ntemp, nfe=ntemp, idx=nidx)
      call fourier(nn,-1_ik,scratch2, nfs=ntemp, nfe=ntemp)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         temp(i,j,k)=scratch2(i,j,k,ntemp)
      end do; end do ;  end do
      !$omp end parallel do

      call calcia(nn, temp)

      nnn(1:3) = nn(1:3)
      nnn(4) = 3_ik
      call copy(nnn, fqr, qr, nfs=1_ik, nfe=3_ik)
      call fourier(nnn, 1_ik, fqr, nfs=1_ik, nfe=3_ik)
      call divergence(nnn, fdivqr_tmp, fqr)

      call zero(nn,scratch2)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         scratch2(i,j,k,ntemp) = fdivqr_tmp(i,j,k)
      end do; end do ;  end do
      !$omp end parallel do

      call shift(nn,-1_ik,scratch2,nfs=ntemp, nfe=ntemp,idx=nidx)

      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         fdivqr(i,j,k) = fdivqr(i,j,k) + scratch2(i,j,k,ntemp)
      end do; end do ;  end do
      !$omp end parallel do

      call zero(nn,scratch2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nidx=2
      call zero(nn,scratch2)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         fdivqr_tmp(i,j,k) = 0.0_rk
      end do; end do ;  end do
      !$omp end parallel do

      call zero(nn,scratch2)
      call copy(nn,scratch2,fu, nfs=ntemp, nfe=ntemp)
      call shift(nn,1_ik, scratch2, nfs=ntemp, nfe=ntemp, idx=nidx)
      call fourier(nn,-1_ik,scratch2, nfs=ntemp, nfe=ntemp)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         temp(i,j,k)=scratch2(i,j,k,ntemp)
      end do; end do ;  end do
      !$omp end parallel do

      call calcia(nn, temp)

      nnn(1:3) = nn(1:3)
      nnn(4) = 3_ik
      call copy(nnn, fqr, qr, nfs=1_ik, nfe=3_ik)
      call fourier(nnn, 1_ik, fqr, nfs=1_ik, nfe=3_ik)
      call divergence(nnn, fdivqr_tmp, fqr)

      call zero(nn,scratch2)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         scratch2(i,j,k,ntemp) = fdivqr_tmp(i,j,k)
      end do; end do ;  end do
      !$omp end parallel do

      call shift(nn,-1_ik,scratch2,nfs=ntemp, nfe=ntemp,idx=nidx)

      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         fdivqr(i,j,k) = fdivqr(i,j,k) + scratch2(i,j,k,ntemp)
      end do; end do ;  end do
      !$omp end parallel do

      call zero(nn,scratch2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nidx=3
      call zero(nn,scratch2)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         fdivqr_tmp(i,j,k) = 0.0_rk
      end do; end do ;  end do
      !$omp end parallel do

      call zero(nn,scratch2)
      call copy(nn,scratch2,fu, nfs=ntemp, nfe=ntemp)
      call shift(nn,1_ik, scratch2, nfs=ntemp, nfe=ntemp, idx=nidx)
      call fourier(nn,-1_ik,scratch2, nfs=ntemp, nfe=ntemp)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         temp(i,j,k)=scratch2(i,j,k,ntemp)
      end do; end do ;  end do
      !$omp end parallel do

      call calcia(nn, temp)

      nnn(1:3) = nn(1:3)
      nnn(4) = 3_ik
      call copy(nnn, fqr, qr, nfs=1_ik, nfe=3_ik)
      call fourier(nnn, 1_ik, fqr, nfs=1_ik, nfe=3_ik)
      call divergence(nnn, fdivqr_tmp, fqr)

      call zero(nn,scratch2)
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         scratch2(i,j,k,ntemp) = fdivqr_tmp(i,j,k)
      end do; end do ;  end do
      !$omp end parallel do

      call shift(nn,-1_ik,scratch2,nfs=ntemp, nfe=ntemp,idx=nidx)

      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         fdivqr(i,j,k) = fdivqr(i,j,k) + scratch2(i,j,k,ntemp)
      end do; end do ;  end do
      !$omp end parallel do

      call zero(nn,scratch2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         fdivqr(i,j,k) = 0.5_rk * fdivqr(i,j,k)
      end do; end do ;  end do
      !$omp end parallel do



      !call fourier(nn, 1_ik, icp, nfs=1_ik, nfe=1_ik)

      return

    end subroutine compute_radiation

    subroutine compute_mhd
      use types
      implicit none

      call zero(nn,scratch)
      call copy(nn,u,fu)
      call fourier(nn,-1_ik,u,nfs=nu1,nfe=nu3)
      call fourier(nn,-1_ik,u,nfs=nb1,nfe=nb3)
      !compute u x b
      !$omp parallel do
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
         ! rhs = u x b
         scratch(i,j,k,nb1) = u(i,j,k,nu2)*u(i,j,k,nb3) - &
              &u(i,j,k,nu3)*u(i,j,k,nb2)
         scratch(i,j,k,nb2) = u(i,j,k,nu3)*u(i,j,k,nb1) - &
              &u(i,j,k,nu1)*u(i,j,k,nb3)
         scratch(i,j,k,nb3) = u(i,j,k,nu1)*u(i,j,k,nb2) - &
              &u(i,j,k,nu2)*u(i,j,k,nb1)
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
            du(i,j,k,nb1) = psu(i,j,k,nu2)*psu(i,j,k,nb3) - &
                 &psu(i,j,k,nu3)*psu(i,j,k,nb2)
            du(i,j,k,nb2) = psu(i,j,k,nu3)*psu(i,j,k,nb1) - &
                 &psu(i,j,k,nu1)*psu(i,j,k,nb3)
            du(i,j,k,nb3) = psu(i,j,k,nu1)*psu(i,j,k,nb2) - &
                 &psu(i,j,k,nu2)*psu(i,j,k,nb1)
         end do; end do ; end do
         !$omp end parallel do
         call fourier(nn,1_ik,du,nfs=nb1,nfe=nb3)
         !Forward Fourier transforms
         call curl(nn,fnls,du,nb1)
      end if


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
      complex(ck) :: ff(nn(4))
      real(rk) :: ksq, fvisc_val
      !$omp parallel do private(ksq,ff,fvisc_val,l)
      do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
         if(isactive(i,j,k)) then
            call forcing_rhs(nn,i,j,k,fu,ff)
            ksq=wv(1_ik,i,j,k)**2+wv(2_ik,i,j,k)**2+wv(3_ik,i,j,k)**2
            do l=1,nn(4)
               fvisc_val=-ksq*visc(l)*fu(i,j,k,l)
               fnl(i,j,k,l)=fnl(i,j,k,l)+fvisc_val+ff(l)
            end do
         else
            fnl(i,j,k,:)=0.0_rks
         end if
      end do; end do; end do
      !$omp end parallel do
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
  

    iii=(i-1)/2+1
    !Forcing condition in Fourier space
    condition=(trk1(iii)/=0.and.trk2(j)/=0.and.trk3(k)/=0).and.&
         &(sqrt(real(trk1(iii)**2+trk2(j)**2+trk3(k)**2,rk))<KFORCING)

    if(condition) then
       do l=1,nn(4)
          ff(l)=fscale(l)*cmplx(fu(i,j,k,l),fu(i+1,j,k,l),ck)
       end do
    else
       ff(:)=cmplx(0.0_rk,0.0_rk,ck)

    end if

    
    if(.not.FORCED) then
       ff(nu1:nu3)=cmplx(0.0_rk,0_rk,ck)
    end if

    if(PASSIVE_SCALAR.and..not.FORCED_PASSIVE_SCALAR) then
       ff(nsclf:nscll)=cmplx(0.0_rk,0_rk,ck)
    end if

    if(MHD.and..not.FORCED_MHD) then
       ff(nb1:nb3)=cmplx(0.0_rk,0_rk,ck)
    end if
    
    return

  end subroutine forcing_rhs

  subroutine lorentz_force(nn,lf,ad,fu)
    use types
    use mpivars
    use data,only:u,scratch,psu,du,nb1,nb2,nb3,nu1,nu2,nu3,copy,zero
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT)    :: lf
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT)    :: ad
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)     :: fu
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

    !Calculate current density
    call curl(nn,scratch,fu,nb1)


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

    !Get current density in physical space
    call fourier(nn,-1_ik,scratch,nfs=nb1,nfe=nb3)
    !Get magnetic field in physical space
    call fourier(nn,-1_ik,u,nfs=nb1,nfe=nb3)

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
    call mpi_allreduce(MPI_IN_PLACE,maxj,1,MPIRK,MPI_MAX,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,maxlf,1,MPIRK,MPI_MAX,MPI_COMM_WORLD,mpierr)
#endif

    !Forward Fourier transforms
    call fourier(nn,1_ik,lf,nfs=nu1,nfe=nu3)

    if(AMB_DIFF) then
       call fourier(nn,1_ik,ad,nfs=nb1,nfe=nb3)
    end if


    !Same for the phase-shifted arrays
    if(DEALIASING==PATTERSON_ORSZAG) then
       call fourier(nn,1_ik,ad,nfs=nu1,nfe=nu3)
       call shift(nn,-1_ik,ad,nfs=nu1,nfe=nu3)
       if(AMB_DIFF) then
          call fourier(nn,1_ik,psu,nu1,nu3)
       end if
       !Get the total (average) terms
       !$omp parallel do
       do l=nu1,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          lf(i,j,k,l) = 0.5*(lf(i,j,k,l) + ad(i,j,k,l))
          if(AMB_DIFF) then
             ad(i,j,k,l)=0.5_rk*(ad(i,j,k,l)+psu(i,j,k,l))
          end if
       end do; end do ; end do ; end do
    end if


    !Truncate the ambipolar diffusion terms
    if(AMB_DIFF.or.HALL) then
       call curl(nn,ad,ad,nu1)
       call truncate(nn,ad)
    end if


    return

  end subroutine lorentz_force

  subroutine mean_dissipation(nn,fu,emeans)
    use data, only: scratch,fsclgrad,fsclgrads,nu1,nu2,nu3,nb1,zero
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

    emeans(:)=visc*emeans(:)

    call zero(nn,scratch)
    
    return

  end subroutine mean_dissipation


  function mean_kinetic_helicity_dissipation(nn,fu) result(mkhdis)
    use data, only: scratch,scratch2,rks1,nu1,nu3,zero
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
    call curl(nn,scratch2,fu,nu1)
    !Get curl of vorticity
    call curl(nn,rks1,scratch2,nu1)
    !Inverse Fourier transform
    call fourier(nn,-1_ik,scratch2)

    call fourier(nn,-1_ik,rks1)

    !Perform inner product
    call inner_product(nn,scratch,nu1,scratch2,rks1,nu1)
    !Forward Fourier transform
    call fourier(nn,1_ik,scratch)
    !Truncate
    call truncate(nn,scratch)
    !Calculate mean value

    tmp=mean_value(nn,scratch)

    mkhdis=2.0_rk*visc(nu1)*tmp(nu1)

    call zero(nn,scratch)
    call zero(nn,scratch2)
    call zero(nn,rks1)
    
    return

  end function mean_kinetic_helicity_dissipation




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

  subroutine divergence(nn,fdiv,fu)
    use parameters, only: dim1
    use data, only: wv,isactive
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3)), intent(OUT) :: fdiv
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:3), intent(IN)  :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate curl of a vector field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                             :: i,j,k
    complex(ck)                                             :: tmpx, tmpy, tmpz, tmp
    
    !$omp parallel do private(tmpx, tmpy, tmpz, tmp)
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
       if(isactive(i,j,k)) then
          tmpx = cmplx(fu(i,j,k,1),fu(i+1,j,k,1),ck)
          tmpy = cmplx(fu(i,j,k,2),fu(i+1,j,k,2),ck)
          tmpz = cmplx(fu(i,j,k,3),fu(i+1,j,k,3),ck)
          tmp=ii*(wv(1_ik, i, j, k)*tmpx + &
               &wv(2_ik, i, j, k)*tmpy +&
               &wv(3_ik, i, j, k)*tmpz)
          fdiv(i, j, k) = real(tmp, rk)
          fdiv(i+1, j, k) = aimag(tmp)
       else
          fdiv(i,j,k)=0.0_rk
       end if
    end do; end do ; end do
    !$omp end parallel do


    return

  end subroutine divergence

  subroutine gradient(nn,fg,fu,nfs,nfe)
    use parameters
    use data, only: wv,isactive, nu1, nu2, nu3, nsclf,nscll, zero
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4),1:3), intent(OUT) :: fg
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

          fg(i,j,k,l,1)=real(ii*wv(1_ik,i,j,k)*tmp,rk)
          fg(i+1,j,k,l,1)=aimag(ii*wv(1_ik,i,j,k)*tmp)

          fg(i,j,k,l,2)=real(ii*wv(2_ik,i,j,k)*tmp,rk)
          fg(i+1,j,k,l,2)=aimag(ii*wv(2_ik,i,j,k)*tmp)

          fg(i,j,k,l,3)=real(ii*wv(3_ik,i,j,k)*tmp,rk)
          fg(i+1,j,k,l,3)=aimag(ii*wv(3_ik,i,j,k)*tmp)
       else
          fg(i,j,k,l,1:3)=0.0_rk
          fg(i+1,j,k,l,1:3)=0.0_rk
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
       if(isactive(i,j,k).and.ksq>small) then
          fa(i,j,k,nf)=ksq**(-1)*scratch(i,j,k,nf)
          fa(i,j,k,nf+1)=ksq**(-1)*scratch(i,j,k,nf+1)
          fa(i,j,k,nf+2)=ksq**(-1)*scratch(i,j,k,nf+2)
       else
          fa(i,j,k,nf:nf+2)=0.0_rk
       end if
    end do; end do ; end do
    !$omp end parallel do

    return

  end subroutine vector_potential


  function mean_magnetic_helicity(nn,fu) result(mmh)
    use types
    use data, only: scratch,scratch2,nb1,nu1,nu2,nu3, zero, u, copy
    implicit none
    real(rk)                                                :: mmh
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)               :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate mean magnetic helicity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                              :: i,j,k,l
    real(rk), dimension(1:nn(4))                             :: tmp
    real(rk)                                                 :: fac
    !Get vector potential of the magnetic field
    call vector_potential(nn,scratch2,fu,nb1)
    !Inverse Fourier transforms
    call fourier(nn,-1_ik,scratch2,nb1,nb3)
    call copy(nn,u,fu)
    call fourier(nn,-1_ik,u,nb1,nb3)
    !Perform inner product
    call inner_product(nn,scratch,nu1,scratch2,&
         &u,nb1)
    !Forward Fourier transform
    call fourier(nn,1_ik,scratch,nu1,nu1)
    !Truncate
    call truncate(nn,scratch,nu1,nu1)

    !$omp parallel do
    do l=nu2,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       scratch(i,j,k,l)=0.0_rk
    end do ; end do ; end do ;  end do
    !$omp end parallel do
    !Calculate mean value
    tmp=mean_value(nn,scratch)

    mmh=tmp(nu1)

    !Set auxiliary arrays back to zero
    call zero(nn,u)
    call zero(nn,scratch)
    call zero(nn,scratch2)


    return

  end function mean_magnetic_helicity


  function mean_magnetic_helicity_dissipation(nn,fu) result(mmhdis)
    use types
    use data, only: u,rmsarr,scratch,nb1,nu1,nu2,nu3,copy,zero
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
    call inner_product(nn,scratch,nu1,u,rmsarr,nb1)

    !Forward Fourier transform
    call fourier(nn,1_ik,scratch,nu1,nu1)

    !Truncate
    call truncate(nn,scratch)

    !$omp parallel do
    do l=nu2,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       scratch(i,j,k,l)=0.0_rk
    end do ; end do ; end do ;  end do
    !$omp end parallel do
    !Calculate mean value
    tmp=mean_value(nn,scratch)

    mmhdis = visc(nb1)*tmp(nu1)

    !Set auxiliary arrays back to zero
    call zero(nn,rmsarr)
    call zero(nn,scratch)

    return

  end function mean_magnetic_helicity_dissipation



  function ambipolar_diffusion_dissipation(nn,fu) result(addis)
    use types
    use mpivars
    use data, only: rmsarr,u,scratch,nb1,nb2,nb3,copy,zero
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
    call mpi_allreduce(MPI_IN_PLACE,addis,1,MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

    addis=AD_COEFF*addis

    !Set auxiliary arrays back to zero
    call zero(nn,rmsarr)
    call zero(nn,scratch)

    return

  end function ambipolar_diffusion_dissipation




  function mean_cross_helicity_dissipation(nn,fu) result(mchdis)
    use types
    use data, only: scratch,scratch2,nu1,nu2,nu3,nb1,nb3,zero
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
    call curl(nn,scratch,fu,nu1)
    !Get current density
    call curl(nn,scratch2,fu,nb1)

    !Inverse Fourier transforms
    call fourier(nn,-1_ik,scratch,nu1,nu3)
    call fourier(nn,-1_ik,scratch2,nb1,nb3)

    !$omp parallel do
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       scratch(i,j,k,nb1:nb3)=scratch(i,j,k,nu1:nu3)
    end do; end do ; end do
    !$omp end parallel do

    !Perform inner product
    call inner_product(nn,scratch2,nu1,scratch,scratch2,nb1)

    !Forward Fourier transform
    call fourier(nn,1_ik,scratch2)

    !Truncate
    call truncate(nn,scratch2)

    !$omp parallel do
    do l=nu2,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       scratch2(i,j,k,l)=0.0_rk
    end do; end do ; end do ;  end do
    !$omp end parallel do
    !Calculate mean value
    tmp=mean_value(nn,scratch2)

    mchdis = (visc(nu1) + visc(nb1))*tmp(1)

    !Set auxiliary arrays back to zero
    call zero(nn,scratch)
    call zero(nn,scratch2)

    return

  end function mean_cross_helicity_dissipation

  function mean_cross_helicity(nn,fu) result(mcross_hel)
    use types
    use data, only: scratch2,scratch,nu1,nu2,nu3,nb1,nb3,zero
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
       scratch(i,j,k,nb1:nb3)=fu(i,j,k,nu1:nu3)
       scratch2(i,j,k,nb1:nb3)=fu(i,j,k,nb1:nb3)
    end do; end do ; end do
    !$omp end parallel do


    !Inverse Fourier transforms
    call fourier(nn,-1_ik,scratch,nb1,nb3)
    call fourier(nn,-1_ik,scratch2,nb1,nb3)

    mcross_hel=0.0_rk

    !Perform inner product
    call inner_product(nn,scratch2,nu1,scratch,scratch2,nb1)
    !Forward Fourier transform
    call fourier(nn,1_ik,scratch2,nu1,nu1)
    !Truncate
    call truncate(nn,scratch2,nu1,nu1)
    !$omp parallel do
    do l=nu2,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       scratch2(i,j,k,l)=0.0_rk
    end do ; end do ; end do ;  end do
    !$omp end parallel do
    !Calculate mean value
    tmp=mean_value(nn,scratch2)
    mcross_hel=0.5_rk*tmp(nu1)

    !Set auxiliary arrays back to zero
    call zero(nn,scratch)
    call zero(nn,scratch2)

    return

  end function mean_cross_helicity

  function mean_kinetic_helicity(nn,fu) result(mkin_hel)
    use types
    use data, only: u,scratch,scratch2,nu1,nu2,copy,zero
    implicit none
    real(rk)                                   :: mkin_hel
    integer(ik), dimension(1:4), intent(in)                   :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)  :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Calculate mean kinetic helicity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk), dimension(:), allocatable :: tmp
    integer(ik) :: i,j,k,l
    real(rk) :: fac

    !Set auxiliary arrays to zero
    call zero(nn,scratch)
    call zero(nn,scratch2)

    if(.not.allocated(tmp)) allocate(tmp(1:nn(4)))
    !Get velocity field
    call copy(nn,u,fu)

    !Get vorticity
    call curl(nn,scratch2,fu,nu1)

    !Inverse Fourier transforms
    call fourier(nn,-1_ik,u,nu1,nu3)
    call fourier(nn,-1_ik,scratch2,nu1,nu3)

    !Perform inner product
    call inner_product(nn,scratch,nu1,u,scratch2,nu1)
    !Forward Fourier transform
    call fourier(nn,1_ik,scratch,nu1,nu1)

    !Truncate
    call truncate(nn,scratch,nu1,nu1)

    !Calculate mean value
    tmp=mean_value(nn,scratch)

    mkin_hel=0.5_rk*tmp(nu1)
    
    !Set auxiliary arrays backto zero
    deallocate(tmp)
    call zero(nn,scratch)
    call zero(nn,scratch2)
    call zero(nn,u)

    return

  end function mean_kinetic_helicity

  subroutine inner_product(nn,cross_hel,nfout,u,b,nfin)
    use types
    use data, only: zero
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
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       cross_hel(i,j,k,nfout) = u(i,j,k,nfin)*b(i,j,k,nfin)+&
            &u(i,j,k,nfin+1)*b(i,j,k,nfin+1)+&
            &u(i,j,k,nfin+2)*b(i,j,k,nfin+2)
    end do; end do ; end do
    !$omp end parallel do

    return

  end subroutine inner_product

  subroutine dissipation(nn,e,nfout,fu)
    use data, only: wv,u,rks1, rmsarr,isactive,nu1,nu2,nu3,zero
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
             tmp1=scale*cmplx(fu(i,j,k,nu1),fu(i+1,j,k,nu1),ck)
             tmp2=scale*cmplx(fu(i,j,k,nu2),fu(i+1,j,k,nu2),ck)
             tmp3=scale*cmplx(fu(i,j,k,nu3),fu(i+1,j,k,nu3),ck)
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
    call mpi_allreduce(MPI_IN_PLACE,area,1,MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
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
    call mpi_allreduce(MPI_IN_PLACE,L,1,MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,diss,1,MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
    if(area==0._rk) area=1.0_rk
    L=(area**(-1))*L
    if(diss==0.0_rk) diss=1.0_rk
    lambda=sqrt(area/diss)
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

  subroutine radiation_spectra(nn,nfilespec,nespec)
    use data, only: ga, qr
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    integer(ik), intent(IN)                                   :: nfilespec
    integer(ik), intent(IN)                                   :: nespec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Energy spectra of the radiative flux q (vector) and the incident
    !radiation G (scalar). Both live in physical space, so transform a
    !local copy to Fourier space and reuse the generic spectrum routines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik)                                :: i,j,k,l
    integer(ik), dimension(1:4)                :: nnn
    real(rks), dimension(:,:,:,:), allocatable :: fwork
    character(256)                             :: fname

    !radiative heat flux q : vector spectrum
    nnn = nn
    nnn(4) = 3_ik
    allocate(fwork(1:dim1(nn(1)),1:nn(2),1:nn(3),1:3))
    fwork = 0.0_rks
    !$omp parallel do
    do l=1,3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       fwork(i,j,k,l) = qr(i,j,k,l)
    end do; end do ; end do ; end do
    !$omp end parallel do
    call fourier(nnn,1_ik,fwork,nfs=1_ik,nfe=3_ik)
    write(fname,'(a,i0.5,a)') 'qrspec.', nfilespec, '.dat'
    call vector_spectrum(nnn,fwork,1_ik,fname,nespec)
    deallocate(fwork)

    !incident radiation G : scalar spectrum
    nnn = nn
    nnn(4) = 1_ik
    allocate(fwork(1:dim1(nn(1)),1:nn(2),1:nn(3),1:1))
    fwork = 0.0_rks
    !$omp parallel do
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       fwork(i,j,k,1) = ga(i,j,k)
    end do; end do ; end do
    !$omp end parallel do
    call fourier(nnn,1_ik,fwork,nfs=1_ik,nfe=1_ik)
    write(fname,'(a,i0.5,a)') 'gspec.', nfilespec, '.dat'
    call scalar_spectrum(nnn,fwork,1_ik,fname,nespec)
    deallocate(fwork)

    return

  end subroutine radiation_spectra



  function incompressibility(nn,fu,nf) result(meandivv)
    use data, only: wv,scratch,nu1,zero
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
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(INOUT) :: fu
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
    use data, only: u, nu1,nu3,nsclf,nscll,nb1,nb3,copy,zero
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

    if(tmp(nu1)/=0.0) then
       !$omp parallel do
       do l=nu1,nu3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
          fu(i,j,k,l)=fu(i,j,k,l)/sqrt(tmp(nu1))
       end do; end do ; end do ; end do
       !$omp end parallel do
    end if

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
          call mpi_allreduce(MPI_IN_PLACE,smax,1,MPIRK,MPI_MAX,MPI_COMM_WORLD,mpierr)
          smin = minval(u(:, :, :, ntemp))
          call mpi_allreduce(MPI_IN_PLACE,smin,1,MPIRK,MPI_MIN,MPI_COMM_WORLD,mpierr)
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))   
             u(i,j,k,ntemp)= (TEMPMAX - TEMPMIN) * (u(i,j,k,ntemp) - smin) / &
                  &(smax - smin) + TEMPMIN
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

  subroutine fourier(nn,dir,u,nfs,nfe,trunc)
    use mpivars
#ifdef _MPI_
    use fft_heffte
    use mpivars
#endif
    implicit none
    integer(ik), dimension(1:4), intent(in)                   :: nn
    integer(ik), intent(IN)                                    :: dir
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(INOUT)              :: u
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
    use data, only: u,fu,rmsarr,scratch,nu1,nu3,nb1,nb3,copy,zero
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
    call mpi_allreduce(MPI_IN_PLACE,mvel,1,MPIRK,MPI_MAX,MPI_COMM_WORLD,mpierr)
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

    if(dt>0.05) dt=0.05
    
    !Set auxiliary arrays back to zero
    call zero(nn,rmsarr)

    return

  end subroutine cfl_condition


end module numerics

