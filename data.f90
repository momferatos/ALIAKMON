!!$     ___ __                                       
!!$ (  / _ \\ \       /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!$Copyright (c) 2009-2020 Georgios Momferatos
module data
  use types
  use parameters
  implicit none
  ! data and initialization
  ! primary arrays
  ! fields in real space
  real(rks), dimension(:,:,:,:), allocatable    :: u
  ! fields in Fourier space
  real(rks), dimension(:,:,:,:), allocatable    :: fu
  ! Lorentz force
  real(rks), dimension(:,:,:,:), allocatable    :: lf
  ! ambipolar diffusion terms
  real(rks), dimension(:,:,:,:), allocatable    :: ad
  ! arrays used for phase shifting
  complex(cks), dimension(:,:,:), allocatable :: phases,iphases   
  ! auxiliary arrays used as temporary storage during dealiasing
  real(rks), dimension(:,:,:,:), allocatable    :: du
  real(rks), dimension(:,:,:,:), allocatable    :: psu
  ! array for the right-hand-side of the equations of motion
  real(rks), dimension(:,:,:,:), allocatable    :: rhs
  ! auxiliary arrays used as temporary storage during Runge-Kutta integration
  real(rks), dimension(:,:,:,:), allocatable    :: rks1, rks2
  real(rks), dimension(:,:,:,:), allocatable    :: rmsarr
  ! auxiliary arrays for the energy conservation test
  real(rks), dimension(:,:,:,:), allocatable    :: arr_en_1
  real(rks), dimension(:,:,:,:), allocatable    :: arr_en_2
  real(rks), dimension(:,:,:,:), allocatable    :: arr_en_3
  real(rks), dimension(:,:,:,:), allocatable    :: scratch
  real(rks), dimension(:,:,:,:), allocatable    :: fnls
  ! auxiliary arrays for the passive scalar gradients
  real(rks), dimension(:,:,:,:,:), allocatable, target  :: fsclgrads
  real(rks), dimension(:,:,:,:), pointer        :: fsclgrad
  ! auxiliary array for slice HDF5 output
  real(rks), dimension(:, :, :), allocatable    :: slice
  ! wavevector arrays
  integer(ik),  dimension(:), allocatable     :: k1, k2, k3,trk1, trk2, trk3,&
       & gk2, trgk2, gk3,trgk3
  ! maximum wavenumber
  real(rk)                                    :: k_max
  ! used for truncation
  logical, dimension(:,:,:),allocatable       :: isactive
  integer(i2b), dimension(:,:), allocatable   :: modes, trunc                   
  integer(ik)                                 :: nespec 
  ! variables used for particles
  real(sp), dimension(:,:), allocatable       :: x,x0,vp
  real(sp), dimension(:), allocatable         :: accel
  integer(i1b), dimension(:), allocatable     :: plbl
  integer(ik)                                 :: np=1000
  integer(ik)                                 :: ncg=4
  integer(ik)                                 :: nrho
  real(rks), dimension(:,:,:), allocatable    :: rho
  real(rks), dimension(:,:), allocatable      :: rho2d
  integer(ik), dimension(1:4)                   :: nn

  interface zero
     module procedure zero3d
     module procedure zero4d
  end interface zero
  
contains

  
  subroutine zero3d(nn,array)
    use types
    implicit none
    integer(ik), dimension(4), intent(in) :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)) :: array
    !
    ! fills array with zeros
    !
    integer(ik) :: i,j,k,l
    
    !$omp parallel do
    do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       array(i,j,k,l)=0.0_rks
    end do; end do ; end do ; end do
    !$omp end parallel do
    
    return
  end subroutine zero3d

  subroutine zero4d(nn,array)
    use types
    implicit none
    integer(ik), dimension(4), intent(in) :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4),1:3), intent(out) :: array
    !
    ! fills array with zeros
    !
    integer(ik) :: i,j,k,l,m
    
    !$omp parallel do
    do m=1,3 ; do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1)) 
       array(i,j,k,l,m)=0.0_rks
    end do; end do ; end do ; end do ; end do
    !$omp end parallel do
    
    return
  end subroutine zero4d

  subroutine copy(nn,dest,source,ilim)
    use types
    implicit none
    integer(ik), dimension(4), intent(in) :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(out) :: dest
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(in) :: source
    integer(ik), optional :: ilim
    !
    ! copies array source to array dest
    !
    integer(ik) :: i,j,k,l,iilim

    iilim=dim1(nn(1))
    if(present(ilim)) iilim=ilim
    !$omp parallel do
    do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,iilim
       dest(i,j,k,l)=source(i,j,k,l)
    end do ; end do ; end do ; end do
    !$omp end parallel do

    return
    
  end subroutine copy
  
  subroutine alloc_init
    use mpivars
    implicit none
    !
    ! allocates and initializes all arrays
    ! 
    integer(ik) :: nmax,i,j,k,n
    real(rk)    :: kk
    logical     :: trunc_crit

    nmax=max(n1,n2,n3)

    ! set number of fields
    ! velocity is certain
    nu1=1
    nu2=2
    nu3=3
    nfields=3
    ! arbitrary number of passive scalars
    if(PASSIVE_SCALAR) then
       nscl=numscls
       nsclf=nu3+1
       nscll=nu3+nscl
       nfields=nu3+nscl
    end if
    ! magnetic field
    if(MHD) then
       nb1=nfields+1
       nb2=nfields+2
       nb3=nfields+3
       nfields=nfields+3
    end if

    ! set vector of array dimensions
    nn(1)=n1
    nn(2)=n2
    nn(3)=n3
    nn(4)=nfields


    ! viscosities, scalar diffusivities, magnetic diffusivities
    allocate(visc(1:nn(4)))
    visc(:)=0.0_rk

    ! forcing scales
    allocate(fscale(1:nn(4)))
    fscale(:)=0.0_rk

    ! passive scalar dissipation rates
    allocate(emeanscl(1:nn(4)))
    emeanscl(:)=0.0_rk

    ! passive scalar variances
    allocate(sclvarprev(1:nn(4)))
    sclvarprev(:)=1.0

    ! allocate truncation mask
    ! true if the mode is active, false if it is inactive
    allocate(isactive(1:dim1(nn(1)),1:nn(2),1:nn(3)))
    isactive(:,:,:) = .false.

    ! allocate secondary arrays
    ! used as scratch
    allocate(scratch(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,scratch)

    ! used for comupting mean-square values in subroutine msvalue
    allocate(rmsarr(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,rmsarr)

    ! used for invariants validation only 
    if(MHD) then
       allocate(arr_en_1(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
       call zero(nn,arr_en_1)
    end if

    if(PASSIVE_SCALAR) then
       allocate(arr_en_2(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
       call zero(nn,arr_en_2)

       allocate(arr_en_3(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
       call zero(nn,arr_en_3)
    end if

    if(PASSIVE_SCALAR) then
       ! gradient of the passive scalars in fourier space
       allocate(fsclgrads(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4),nu1:nu3))
       call zero(nn,fsclgrads)
    end if

    if(DEALIASING==PATTERSON_ORSZAG) then
       ! phase-shifted non-linear terms
       allocate(fnls(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
       call zero(nn,fnls)
    end if

    !allocate particle arrays
    if(PARTICLES) then
       np=NPART
       allocate(x0(3,np))
       allocate(accel(np))
       allocate(x(3,np))
       allocate(rho(1:dim1(n1),1:n2,1:n3))
       allocate(rho2d(1:n1,1:n2))
       if(INERTIAL) then
          allocate(vp(3,np))
       end if
       if(MIXING) then
          allocate(plbl(1:np))
       end if
    end if

    !allocate primary arrays
    ! fields in physical space
    allocate(u(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,u)
    ! fields in fourier space
    allocate(fu(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,fu)
    ! right-hand side of the equations of motion
    allocate(rhs(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,rhs)
    ! scratch used for time integration
    allocate(rks1(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,rks1)

    !only for ambipolar diffusion or hall effect
    if(MHD.and.AMB_DIFF.or.HALL) then
       ! ambipolar diffusion or Hall effect terms
       allocate(ad(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nfields))
       call zero(nn,ad)
    end if

    !phases array for phase-shifting
    if(DEALIASING/=NONE.and.STORE_PHASES) then
       ! forward phase factors
       allocate(phases(dim1(n1),n2,n3))
       phases(:,:,:) = 0.0_rk
       ! inverse phase factors
       allocate(iphases(dim1(n1),n2,n3))
       iphases(:,:,:) = 0.0_rk
    end if
    !only for Patterson-Orszag dealiasing
    if(DEALIASING==PATTERSON_ORSZAG) then
       ! phase-shifted fields in physical space
       allocate(psu(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nfields))
       call zero(nn,psu)
       ! scratch of other phase-shifted terms
       allocate(du(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nfields))
       call zero(nn,du)
    end if

    !only for Runge-Kutta integration
    if(INTEGRATION_METHOD==MRUNGE_KUTTA4) then
       ! scratch used for 4th order Runge-Kutta
       allocate(rks2(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
       call zero(nn,rks2)
    end if

    nespec=max(n1,n2,n3)  

    !allocate wave-vector arrays
    ! x
    allocate(k1(1:dim1(n1)))
    k1(:) = 0.0_rk
    ! x, used for truncation
    allocate(trk1(1:dim1(n1)))
    trk1(:) = 0.0_rk


    ! if we're using MPI, z direction is spread across processes,
    ! each process has a slice of z-width lksize
#ifdef _MPI_
    ! y
    allocate(k2(1:ljsize))
    ! y, used for truncation
    allocate(trk2(1:ljsize))
    ! z
    allocate(k3(1:lksize))
    ! z, used for truncation
    allocate(trk3(1:lksize))
#else
    ! y
    allocate(k2(1:n3))
    ! y, used for truncation
    allocate(trk2(1:n3))
    ! z
    allocate(k3(1:n3))
    ! z, used for truncation
    allocate(trk3(1:n3))
#endif

    k3(:) = 0.0_rk
    trk3(:) = 0.0_rk
    ! if we're using MPI, we also need the global y-wavevector array
    allocate(gk2(1:gn2))
    gk2(:) = 0.0_rk
    ! and its counterpart for truncation
    allocate(trgk2(1:gn2))
    trgk2(:) = 0.0_rk
    ! if we're using MPI, we also need the global z-wavevector array
    allocate(gk3(1:gn3))
    gk3(:) = 0.0_rk
    ! and its counterpart for truncation
    allocate(trgk3(1:gn3))
    trgk3(:) = 0.0_rk
    !initialize wave-vector arrays
    call wave_vectors(k1,k2,k3,trk1,trk2,trk3,gk2,trgk2,gk3,trgk3)

    ! initialize the mask of active modes
    n = maxval(nn(1:3))
    k_max = 0.0 
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1
       ! truncation criterion
       select case(TRUNCATION)
          ! two-thirds rule
       case(TWO_THIRDS)
          trunc_crit = twothirds_tr(i, j, k, n / 3)
          ! spherical
       case(SPHERICAL)
          trunc_crit = spherical_tr(i,j,k,TRFAC*n)
          ! polyhedral
       case(POLYHEDRAL)
          trunc_crit = polyhedral_tr(i,j,k)
       end select

       ! truncate
       if(trunc_crit) then
          ! real part of the mode
          isactive(i,j,k)=.false.
          ! imaginary part 
          isactive(i+1,j,k)=.false.
       else 
          ! wavevector modulus
          kk = sqrt(real(abs(wv(1_ik,i,j,k)**2+wv(2_ik,i,j,k)**2&
               &+wv(3_ik,i,j,k)**2)))
          ! maximum wavevector
          k_max=max(k_max,kk)
          ! real part of the mode
          isactive(i,j,k)=.true.
          ! imaginary part 
          isactive(i+1,j,k)=.true.  
       end if
    end do; end do ; end do

    return

  end subroutine alloc_init

  subroutine alloc_init_post_proc
    use mpivars
    implicit none
    !
    ! allocates and initializes all arrays
    ! 
    integer(ik) :: nmax,i,j,k,n
    real(rk)    :: kk
    logical     :: trunc_crit

    nmax=max(n1,n2,n3)

    ! set number of fields
    ! velocity is certain
    nu1=1
    nu2=2
    nu3=3
    nfields=3
    ! arbitrary number of passive scalars
    if(PASSIVE_SCALAR) then
       nscl=numscls
       nsclf=nu3+1
       nscll=nu3+nscl
       nfields=nu3+nscl
    end if
    ! magnetic field
    if(MHD) then
       nb1=nfields+1
       nb2=nfields+2
       nb3=nfields+3
       nfields=nfields+3
    end if

    ! set vector of array dimensions
    nn(1)=n1
    nn(2)=n2
    nn(3)=n3
    nn(4)=nfields


    ! viscosities, scalar diffusivities, magnetic diffusivities
    allocate(visc(1:nn(4)))
    visc(:)=0.0_rk

    ! forcing scales
    allocate(fscale(1:nn(4)))
    fscale(:)=0.0_rk

    ! passive scalar dissipation rates
    allocate(emeanscl(1:nn(4)))
    emeanscl(:)=0.0_rk

    ! passive scalar variances
    allocate(sclvarprev(1:nn(4)))
    sclvarprev(:)=1.0

    ! allocate truncation mask
    ! true if the mode is active, false if it is inactive
    allocate(isactive(1:dim1(nn(1)),1:nn(2),1:nn(3)))
    isactive(:,:,:) = .false.

    ! allocate secondary arrays
    ! used as scratch
    allocate(scratch(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,scratch)

    allocate(rks1(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,rks1)

    allocate(rmsarr(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,rmsarr)

    !allocate primary arrays
    ! fields in physical space
    allocate(u(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,u)
    ! fields in fourier space
    allocate(fu(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,fu)

    nespec=max(n1,n2,n3)  

    !allocate wave-vector arrays
    ! x
    allocate(k1(1:dim1(n1)))
    k1(:) = 0.0_rk
    ! x, used for truncation
    allocate(trk1(1:dim1(n1)))
    trk1(:) = 0.0_rk

    
    ! if we're using MPI, z direction is spread across processes,
    ! each process has a slice of z-width lksize
#ifdef _MPI_
    ! y
    allocate(k2(1:ljsize))
    ! y, used for truncation
    allocate(trk2(1:ljsize))
    ! z
    allocate(k3(1:lksize))
    ! z, used for truncation
    allocate(trk3(1:lksize))
#else
    ! y
    allocate(k2(1:n3))
    ! y, used for truncation
    allocate(trk2(1:n3))
    ! z
    allocate(k3(1:n3))
    ! z, used for truncation
    allocate(trk3(1:n3))
#endif

    k3(:) = 0.0_rk
    trk3(:) = 0.0_rk
    ! if we're using MPI, we also need the global y-wavevector array
    allocate(gk2(1:gn2))
    gk2(:) = 0.0_rk
    ! and its counterpart for truncation
    allocate(trgk2(1:gn2))
    trgk2(:) = 0.0_rk
    ! if we're using MPI, we also need the global z-wavevector array
    allocate(gk3(1:gn3))
    gk3(:) = 0.0_rk
    ! and its counterpart for truncation
    allocate(trgk3(1:gn3))
    trgk3(:) = 0.0_rk
    
    !initialize wave-vector arrays
    call wave_vectors(k1,k2,k3,trk1,trk2,trk3,gk2, trgk2, gk3,trgk3)

    ! initialize the mask of active modes
    n = maxval(nn(1:3))
    k_max = 0.0 
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1
       ! truncation criterion
       select case(TRUNCATION)
          ! two-thirds rule
       case(TWO_THIRDS)
          trunc_crit = twothirds_tr(i, j, k, n / 3)
          ! spherical
       case(SPHERICAL)
          trunc_crit = spherical_tr(i,j,k,TRFAC*n)
          ! polyhedral
       case(POLYHEDRAL)
          trunc_crit = polyhedral_tr(i,j,k)
       end select

       ! truncate
       if(trunc_crit) then
          ! real part of the mode
          isactive(i,j,k)=.false.
          ! imaginary part 
          isactive(i+1,j,k)=.false.
       else 
          ! wavevector modulus
          kk = sqrt(real(abs(wv(1_ik,i,j,k)**2+wv(2_ik,i,j,k)**2&
               &+wv(3_ik,i,j,k)**2)))
          ! maximum wavevector
          k_max=max(k_max,kk)
          ! real part of the mode
          isactive(i,j,k)=.true.
          ! imaginary part 
          isactive(i+1,j,k)=.true.  
       end if
    end do; end do ; end do

    return

  end subroutine alloc_init_post_proc
  pure function twothirds_tr(i,j,k,kmax)
    implicit none
    logical                 :: twothirds_tr
    integer(ik), intent(IN) :: i,j,k,kmax
    !
    ! two-thirds rule truncation
    !
    integer(ik)             :: iii,m

    iii=(i-1)/2+1

    m=0
    if(abs(trk1(iii))>kmax) then
       m=m+1
    end if

    if(abs(trk2(j))>kmax) then
       m=m+1
    end if

    if(abs(trk3(k))>kmax) then
       m=m+1
    end if


    if(m==1.or.m==2.or.m==3) then
       twothirds_tr=.true.
    else if(m==0) then
       twothirds_tr=.false.
    end if


    return

  end function twothirds_tr

  pure function polyhedral_tr(i,j,k)
    implicit none
    logical                 :: polyhedral_tr
    integer(ik), intent(IN) :: i,j,k
    !
    ! polyhedral truncation
    !
    integer(ik)             :: iii,n,m
    logical                 :: cond

    n=max(n1,n2,gn3)
    iii=(i-1)/2+1

    m=0
    if(abs(trk1(iii))>=n/2.or.abs(trk2(j))>=n/2.or.abs(trk3(k))>=n/2) then
       polyhedral_tr=.true.
    else 
       cond=(abs(trk1(iii)+trk2(j))>=(2*n)/3).or.&
            &(abs(trk2(j)+trk3(k))>=(2*n)/3).or.&
            &(abs(trk1(iii)+trk3(k))>=(2*n)/3).or.&
            &(abs(trk1(iii)-trk2(j))>=(2*n)/3).or.&
            &(abs(trk2(j)-trk3(k))>=(2*n)/3).or.&
            &(abs(trk1(iii)-trk3(k))>=(2*n)/3)
       if(cond) then
          polyhedral_tr=.true.
       else
          polyhedral_tr=.false.
       end if
    end if

    return

  end function polyhedral_tr

  pure function spherical_tr(i,j,k,kmax)
    implicit none
    logical                 :: spherical_tr
    integer(ik), intent(IN) :: i,j,k
    real(rk), intent(IN)    :: kmax
    !
    ! spherical truncation with radius kmax
    !
    integer(ik)             :: iii
    iii=(i-1)/2+1

    if(sqrt(real(trk1(iii)**2+trk2(j)**2+trk3(k)**2,rk))>=kmax) then
       spherical_tr=.true.
    else
       spherical_tr=.false.
    end if

    return

  end function spherical_tr

  pure elemental function wv(dim,i,j,k)
    implicit none
    integer(ik), intent(IN) :: dim,i,j,k
    !
    ! returns wavevector component
    !
    integer(ik)             :: wv

    ! ATTENTION: FFTW inverts the last two dimensions or all arrays in
    ! Fourier space
#ifdef _MPI_
    if(dim==1) then
       wv=k1((i-1)/2+1)
    else if(dim==2) then
       wv=k3(k)
    else if(dim==3) then
       wv=k2(j)
    end if
#else
    if(dim==1) then
       wv=k1((i-1)/2+1)
    else if(dim==2) then
       wv=k2(j)
    else if(dim==3) then
       wv=k3(k)
    end if
#endif

    return

  end function wv

  subroutine wave_vectors(k1,k2,k3,trk1,trk2,trk3,gk2,trgk2,gk3,trgk3)
    use mpivars
    implicit none
    integer(ik), dimension(1:), intent(OUT) :: k1,k2,k3,trk1,trk2,trk3,&
         &gk2, trgk2, gk3,trgk3
    !
    ! sets up wavevectors in standard FFTW order
    !
    integer(ik)                             :: i, j
    integer(ik)                             :: n1,n2,n3


    n1 = nn(1)
    n2 = nn(2)
    n3 = nn(3)

    ! x-direction
    k1(1:dim1(n1))=0_ik
    do i=2,n1/2
       k1(i)=(i-1)
    end do

    k1(n1/2+1)=(n1/2)
    do i=1,n1/2-1
       j=-n1/2+i
       k1(i+n1/2+1)=j
    end do
    trk1(1:dim1(n1))=k1(1:dim1(n1))
    trk1(n1/2+1)=(n1/2)

!!$    ! y-direction
!!$    k2(1)=0_ik
!!$    do i=2,n2/2
!!$       k2(i)=(i-1)
!!$    end do
!!$
!!$    k2(n2/2+1)=(n2/2)
!!$    do i=1,n2/2-1
!!$       j=-n2/2+i
!!$       k2(i+n2/2+1)=j
!!$    end do
!!$    trk2(:)=k2(:)
!!$    trk2(n2/2+1)=(n2/2)
    ! y-direction is divided across MPI processes
    ! global wave-vector array
    gk2(1)=0_ik
    do i=2,gn2/2
       gk2(i)=(i-1)
    end do

    gk2(gn2/2+1)=(gn2/2)
    do i=1,gn2/2-1
       j=-gn2/2+i
       gk2(i+gn2/2+1)=j
    end do
    trgk2(:)=gk2(:)
    trgk2(gn2/2+1)=(gn2/2)
#ifdef _MPI_
    do i=1,ljsize
       ! local wave-vector array
       k2(i)=gk2(i+ljstart-1)
       trk2(i)=trgk2(i+ljstart-1)
    end do
#else
    do i=1,n2
       k2(i)=gk2(i)
       trk2(i)=trgk2(i)
    end do
#endif
    
    ! z-direction is divided across MPI processes
    ! global wave-vector array
    gk3(1)=0_ik
    do i=2,gn3/2
       gk3(i)=(i-1)
    end do

    gk3(gn3/2+1)=(gn3/2)
    do i=1,gn3/2-1
       j=-gn3/2+i
       gk3(i+gn3/2+1)=j
    end do
    trgk3(:)=gk3(:)
    trgk3(gn3/2+1)=(gn3/2)
#ifdef _MPI_
    do i=1,lksize
       ! local wave-vector array
       k3(i)=gk3(i+lkstart-1)
       trk3(i)=trgk3(i+lkstart-1)
    end do
#else
    do i=1,n3
       k3(i)=gk3(i)
       trk3(i)=trgk3(i)
    end do
#endif  


    ! set maximum wavenumber
    if(TRUNCATION==SPHERICAL&
         &.or.TRUNCATION==POLYHEDRAL) then
       kmax=TRFAC*max(n1,gn2,gn3)
    else
       kmax=max(n1/3,n2/3,n3/3)
       if(MHD.and.AMB_DIFF) kmax=max(n1/3,n2/3,n3/3)
    end if

    return

  end subroutine wave_vectors

  subroutine make_phases_array(nn,phases,iphases) 
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    complex(cks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3)), intent(OUT) :: phases
    complex(cks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3)), intent(OUT) :: iphases
    !
    ! prepares phases for Patterson & Orszag (1972) dealiasing
    !
    integer(ik)                                  :: i,j,k,n

    n=maxval(nn(1:3))

    do k=1,nn(3)
       do j=1,nn(2)
          do i=1,dim1(nn(1))
             ! forward phase
             phases(i,j,k)=exp(ii*(PI/real(n,rk))*(wv(1_ik,i,j,k)+&
                  &wv(2_ik,i,j,k)+wv(3_ik,i,j,k)))
             ! inverse phase
             if(abs(phases(i,j,k)) > small) then
                iphases(i,j,k)=1.0_rks/phases(i,j,k)
             else
                iphases(i,j,k)=0.0_rk
             end if
          end do
       end do
    end do

    return

  end subroutine make_phases_array

end module data
