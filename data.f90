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
  ! pressure
  real(rks), dimension(:,:,:),   allocatable    :: press
  !$acc declare create(press)
  ! fields in real space
  real(rks), dimension(:,:,:,:), allocatable    :: u
  ! fields in Fourier space
  real(rks), dimension(:,:,:,:), allocatable    :: fu
  ! Lorentz force
  real(rks), dimension(:,:,:,:), allocatable    :: lf
  ! ambipolar diffusion terms
  real(rks), dimension(:,:,:,:), allocatable    :: ad
  ! arrays used for phase shifting
  complex(cks), dimension(:,:,:,:), allocatable :: phases, iphases
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
  real(rks), dimension(:,:,:,:), allocatable    :: scratch2
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
  !$acc declare create(nn)

  real(rks), dimension(:, :, :, :), allocatable, target ::  ia    ! radiation intensity
  !$acc declare create(ia)
  real(rks), dimension(:, :, :, :), allocatable, target ::  iba   ! radiation intensity at previous iteration
  !$acc declare create(iba)
  real(rks), dimension(:, :, :), allocatable, target ::  temp   ! temperature buffer
  !$acc declare create(temp)
  real(rks), dimension(:, :, :), allocatable, target ::  ga   ! incident radiation
  !$acc declare create(ga)
  real(rks), dimension(:, :, :, :), allocatable, target ::  qr   ! incident radiation
  !$acc declare create(qr)
  real(rks), dimension(:, :, :, :), allocatable, target ::  fqr   ! incident radiation
  real(rks), dimension(:, :, :), allocatable, target ::  fdivqr
  real(rks), dimension(:, :, :), allocatable, target ::  fdivqr_tmp
  real(rk), dimension(:, :), allocatable :: s      ! vector pointing at the direction of the center of the angular
  !$acc declare create(s)
  real(rk), dimension(:), allocatable :: omeg                                 ! solid angle per angular finite volume
  !$acc declare create(omeg)
  integer(ik) :: nnphi
  ! sign vectors used to check the direction of the shat vector
  real(rk), dimension(8, 3) :: sgn
  !$acc declare create(sgn(1:8, 1:3))
  integer(ik), dimension(8) :: istart, iend, istep ! starts, ends and steps of the i-sweep (step is +/- 1)
  !$acc declare create(istart(1:8), istep(1:8), iend(1:8))
  integer(ik), dimension(8) :: jstart, jstep, jend ! starts, ends and steps of the j-sweep (step is +/- 1)
  !$acc declare create(jstart(1:8), jstep(1:8), jend(1:8))
  integer(ik), dimension(8) :: kstart, kstep, kend ! starts, ends and steps of the k-sweep (step is +/- 1)
  !$acc declare create(kstart(1:8), kstep(1:8), kend(1:8))
  real(rks), dimension(:, :, :), allocatable :: sendbuf, recvbuf
  real(rks), dimension(:, :, :), allocatable :: ghostleft
  !$acc declare create(ghostleft)
  real(rks), dimension(:, :, :), allocatable :: ghostright
  !$acc declare create(ghostright)
  real(rks), dimension(:, :, :), allocatable :: left
  !$acc declare create(left)
  real(rks), dimension(:, :, :), allocatable :: right
  !$acc declare create(right)
  logical, dimension(:, :), allocatable :: is_wq
  !$acc declare create(is_wq)
  interface zero
     module procedure zero3d
     module procedure zero4d
  end interface zero

contains


  subroutine allocate_fvdom
    use types
    use parameters
    implicit none
    integer(ik) :: nphi
!!!!!!!!!!!!!!!!!!!!
    ! allocates memory !
!!!!!!!!!!!!!!!!!!!!

    ! find the closest integer nphi for which nsects = nphi * (nphi + 2)
    if(EQSECTS == 2) then
       nphi = int(sqrt(real(nsects, 8)), 8)
       nphi = 2*int(nphi / 2., 8)
       nsects = nphi * (nphi + 2)
       nnphi = nphi
    end if

    !allocate memory
    if(.not.allocated(ia)) allocate(ia(1:nn(1), 1:nn(2), 0:nn(3)+1, 1:nsects))
    !$acc enter data create(ia(1:nn(1), 1:nn(2), 0:nn(3) + 1, 1:nsects))
    if(.not.allocated(iba)) allocate(iba(1:nn(1), 1:nn(2), 0:nn(3)+1, &
         &1:nsects))
    !$acc enter data create(iba(1:nn(1), 1:nn(2), 0:nn(3)+1, 1:nsects))
    if(.not.allocated(temp)) allocate(temp(1:nn(1), 1:nn(2), 1:nn(3)))
    !$acc enter data create(temp(1:nn(1), 1:nn(2), 1:nn(3)))
    if(.not.allocated(ga)) allocate(ga(1:nn(1), 1:nn(2), 1:nn(3)))
!!$acc enter data create(ga(1:nn(1), 1:nn(2), 1:nn(3)))
    if(.not.allocated(qr)) allocate(qr(1:nn(1), 1:nn(2), 1:nn(3), 1:3))
!!$acc enter data create(qr(1:nn(1), 1:nn(2), 1:nn(3), 1:3))
    if(.not.allocated(fqr)) allocate(fqr(1:dim1(nn(1)), 1:nn(2), 1:nn(3), 1:3))
    if(.not.allocated(fdivqr)) allocate(fdivqr(1:dim1(nn(1)), 1:nn(2),&
         & 1:nn(3)))
    if(.not.allocated(fdivqr_tmp)) allocate(fdivqr_tmp(1:dim1(nn(1)), 1:nn(2),&
         & 1:nn(3)))
    if(.not.allocated(s)) allocate(s(1:nsects, 1:3))
    !$acc enter data create(s(1:nsects, 1:3))
    if(.not.allocated(omeg)) allocate(omeg(1:nsects))
    !$acc enter data create(omeg(1:nsects))

    !$acc enter data create(istart(1:8), iend(1:8),  jstart(1:8), jend(1:8),  kstart(1:8), kend(1:8), &
    !$acc& istep(1:8), jstep(1:8), kstep(1:8), sgn(1:8, 1:3))

    if(.not.allocated(sendbuf)) allocate(sendbuf(1:n1, 1:n2, 1:nsects))
    if(.not.allocated(recvbuf)) allocate(recvbuf(1:n1, 1:n2, 1:nsects))

    if(.not.allocated(ghostleft)) allocate(ghostleft(1:n1, 1:n2, 1:nsects))
    !$acc enter data create(ghostleft(1:n1, 1:n2, 1:nsects))
    if(.not.allocated(ghostright)) allocate(ghostright(1:n1, 1:n2, 1:nsects))
    !$acc enter data create(ghostright(1:n1, 1:n2, 1:nsects))

    if(.not.allocated(left)) allocate(left(1:n1, 1:n2, 1:nsects))
    !$acc enter data create(left(1:n1, 1:n2, 1:nsects))
    if(.not.allocated(right)) allocate(right(1:n1, 1:n2, 1:nsects))
    !$acc enter data create(right(1:n1, 1:n2, 1:nsects))

    if(.not.allocated(is_wq)) allocate(is_wq(1:8, 1:nsects))
    !$acc enter data create(is_wq(1:8, 1:nsects))

    return

  end subroutine allocate_fvdom

  subroutine deallocate_fvdom

    implicit none
!!!!!!!!!!!!!!!!!!!!!!
    ! deallocates memory !
!!!!!!!!!!!!!!!!!!!!!!

    if(allocated(ia)) deallocate(ia)
    !$acc exit data delete(ia)
    if(allocated(iba)) deallocate(iba)
    !$acc exit data delete(iba)
    if(allocated(temp)) deallocate(temp)
    !$acc exit data delete(temp)
    if(allocated(ga)) deallocate(ga)
!!$acc exit data delete(ga)
    if(allocated(qr)) deallocate(qr)
    if(allocated(fqr)) deallocate(fqr)
    if(allocated(fdivqr)) deallocate(fdivqr)
    if(allocated(fdivqr_tmp)) deallocate(fdivqr_tmp)
!!$acc exit data delete(qr)
    if(allocated(s)) deallocate(s)
    !$acc exit data delete(s)
    if(allocated(omeg)) deallocate(omeg)
    !$acc exit data delete(omeg)

    !$acc exit data delete(istart(1:8), iend(1:8),  jstart(1:8), jend(1:8),  kstart(1:8), kend(1:8), &
    !$acc& istep(1:8), jstep(1:8), kstep(1:8), sgn(1:8, 1:3))

    if(allocated(sendbuf)) deallocate(sendbuf)
    if(allocated(recvbuf)) deallocate(recvbuf)

    if(allocated(ghostleft)) deallocate(ghostleft)
    !$acc exit data delete(ghostleft)
    if(allocated(ghostright)) deallocate(ghostright)
    !$acc exit data delete(ghostright)

    if(allocated(left)) deallocate(left)
    !$acc exit data delete(left)
    if(allocated(right)) deallocate(right)
    !$acc exit data delete(right)

    if(allocated(is_wq)) deallocate(is_wq)
    !$acc exit data delete(is_wq)

    return

  end subroutine deallocate_fvdom

  subroutine init_fvdom
    use parameters
    use mpivars, only: ljstart, ljsize, lkstart, lksize
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initializes all relevant module variables and arrays !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: m, mp, ns, i, j, k, l
    integer(ik) :: isl, iel, jsl, jel, ksl, kel
    real(rk), dimension(3) :: shat
    integer(ik) :: sd

    !nsects = 80
    
    ! initialize radiative intensities to zero
    !$omp parallel do 
    do l=1,nsects ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       ia(i, j, k, l) = 0.0_rk
    end do; end do ; end do ; end do
    !$omp end parallel do
    !$acc update device(ia(1:n1, 1:n2, 0:nn(3)+1, 1:nsects))

    !$omp parallel do 
    do l=1,nsects ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       iba(i, j, k, l) = 0.0_rk
    end do; end do ; end do ; end do
    !$omp end parallel do
    !$acc update device(iba(1:n1, 1:n2, 1:nn(3), 1:nsects))

    !$omp parallel do 
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       ga(i, j, k) = 0.0_rk
    end do; end do ; end do
    !$omp end parallel do

    !$omp parallel do 
    do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       fdivqr(i, j, k) = 0.0_rk
       fdivqr_tmp(i, j, k) = 0.0_rk
    end do; end do ; end do
    !$omp end parallel do


    !$omp parallel do 
    do l=1,3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       qr(i, j, k, l) = 0.0_rk
    end do; end do ; end do; end do
    !$omp end parallel do

    !$omp parallel do 
    do l=1,3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       fqr(i, j, k, l) = 0.0_rk
    end do; end do ; end do; end do
    !$omp end parallel do

    !$omp parallel do
    do ns=1,nsects
       do j=1,n2
          do i=1,n1
             sendbuf(i, j, ns) = 0.0_rks
             recvbuf(i, j, ns) = 0.0_rks
             ghostleft(i, j, ns) = 0.0_rks
             ghostright(i, j, ns) = 0.0_rks
             left(i, j, ns) = 0.0_rks
             right(i, j, ns) = 0.0_rks
          end do
       end do
    end do
    !$omp end parallel do

    ! perform subdivision of the unit sphere in equal-area angular control volumes        
    call sectors
    !$acc update device(s(1:nsects,1:3), omeg(1:nsects))

    isl = 1
    iel = nn(1)

    istart(1:8) = (/isl,  isl,  isl,  isl,  iel, iel, iel, iel/)  ! starts of the i-sweep
    iend(1:8) = (/iel, iel, iel, iel, isl,  isl,  isl,  isl /)  ! ends of the i-sweep
    istep(1:8) = (/1,  1,  1,  1, -1, -1, -1, -1 /)  ! steps of the i-sweep

    jsl = 1
    jel = nn(2)

    jstart(1:8) = (/jsl,  jsl,  jel, jel, jsl,  jsl,  jel, jel/)  ! starts of the j-sweep
    jend(1:8) = (/jel, jel, jsl,  jsl,  jel, jel, jsl,  jsl /)  ! ends of the j-sweep
    jstep(1:8) = (/1,  1, -1, -1,  1,  1, -1, -1 /)  ! steps of the j-sweep

    ksl = 1
    kel = nn(3)

    kstart(1:8) = (/ksl,  kel, ksl,  kel, ksl,  kel, ksl,  kel/)  ! starts of the k-sweep
    kend(1:8) = (/kel, ksl,  kel, ksl,  kel, ksl,  kel, ksl /)  ! ends of the k-sweep
    kstep(1:8) = (/1, -1,  1, -1,  1, -1,  1, -1 /)  ! steps of the k-sweep

    sgn(1, 1:3) = (/  1,  1,  1 /)
    sgn(2, 1:3) = (/  1,  1, -1 /)
    sgn(3, 1:3) = (/  1, -1,  1 /)
    sgn(4, 1:3) = (/  1, -1, -1 /)
    sgn(5, 1:3) = (/ -1,  1,  1 /)
    sgn(6, 1:3) = (/ -1,  1, -1 /)
    sgn(7, 1:3) = (/ -1, -1,  1 /)
    sgn(8, 1:3) = (/ -1, -1, -1 /)
    
    !$acc update device(istart(1:8), iend(1:8),  jstart(1:8), jend(1:8),  kstart(1:8), kend(1:8), &
    !$acc& istep(1:8), jstep(1:8), kstep(1:8), sgn(1:8, 1:3))

    do ns=1,nsects
       do sd=1,8
          ! direction cosines of s
          shat = s(ns, 1:3) / sqrt(dot_product(s(ns, 1:3), s(ns, 1:3))) 
          is_wq(sd, ns) = all(sgn(sd, 1:3) * shat(1:3) >= 0.0)
       end do
    end do
    !$acc update device(is_wq(1:8, 1:nsects))

    return

  end subroutine init_fvdom

  subroutine sectors
    implicit none

    if(EQSECTS == 1) then
       call sectors_equal ! exactly equal angular control volumes
    else if(EQSECTS == 0) then
       call sectors_unequal ! unequal angular control volumes
    else if(EQSECTS == 2) then
       call sectors_almost_equal ! almost equal angular control volumes
    end if

  end subroutine sectors

  subroutine sectors_unequal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! subdivision of the unit sphere into nsects = nphi * ntheta sectors
    ! of unequal area
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    real(rk) :: theta, phi, dphi, dtheta, phi0, theta0,&
         &phi1, phi2, theta1, theta2 
    integer(ik) :: nphi, ntheta, ns, np, nt

    ! find integer number of subdivisions across phi, theta
    nphi = int(sqrt(real(nsects, 8)), 8)
    ntheta = int(sqrt(real(nsects, 8)), 8)

    ! new total number of directions
    nsects = nphi * ntheta

    ! phi, theta increments
    dphi = PI / nphi
    dtheta = 2 * PI / ntheta

    ! phi, theta starting points
    phi0 = 0.0
    theta0 = 0.0

    ns  = 0
    do np=1,nphi
       phi1 = phi0 + (np - 1) * dphi ! phi startpoint of sector
       phi2 = phi0 + np * dphi ! phi endpoint of sector
       do nt=1,ntheta
          theta1 = theta0 + (nt - 1) * dtheta ! theta startpoint of sector
          theta2 = theta0 + nt * dtheta ! theta endpoint of sector
          ns = ns  + 1 ! next sector

          ! total sector area
          omeg(ns) = (theta2 - theta1) * (cos(phi1) - cos(phi2)) 

          ! sector vector components (vector has magnitude equal to sector area) 
          s(ns, 1) = 0.5 * (phi1 - phi2 - cos(phi1) * sin(phi1) + &
               &cos(phi2) * sin(phi2)) * (sin(theta1) - sin(theta2))

          s(ns, 2) = 0.25 * (cos(theta1) - cos(theta2)) * (-2.0 * phi1 + &
               & 2.0 * phi2 + sin(2 * phi1) - sin(2 * phi2))

          s(ns, 3) = -0.25 * (theta1 - theta2) * (cos(2.0 * phi1) - &
               &cos(2.0 * phi2))

       end do
    end do

    nsects = ns
    print *, ns
    ! write vectors to file for checks
    open(789, file = 's.dat', form = 'formatted')
    do ns=1,nsects
       write(789, '(3e16.8)') s(ns, :) 
    end do
    close(789, status = 'keep')

    return

  end subroutine sectors_unequal

  subroutine sectors_almost_equal
    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! subdivision of the unit sphere into nsects = nphi * (nphi + 2) sectors
    ! of almost equal area
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk) :: theta, phi, dphi, dtheta, phi0, theta0,&
         &phi1, phi2, theta1, theta2 
    integer(ik) :: nphi, ntheta, ns, np, nt, nphi_middle, nn

    nphi = nnphi ! number of phi subdivisions
    nphi_middle = nphi / 2 ! middle phi subdivision
    dphi = PI / nphi ! phi increment

    ! phi, theta starting points
    phi0 = 0.0
    theta0 = 0.0

    ns  = 0
    ntheta = 4 ! starting value of theta subdivisions
    do np=1,nphi
       phi1 = phi0 + (np - 1) * dphi ! phi startpoint of sector 
       phi2 = phi0 + np * dphi ! phi endpoint of sector
       dtheta = 2 * PI / ntheta ! theta increment
       do nt=1,ntheta
          theta1 = theta0 + (nt - 1) * dtheta ! theta startpoint of sector 
          theta2 = theta0 + nt * dtheta ! theta endpoint of sector
          ns = ns + 1 ! next sector

          ! sector area
          omeg(ns) = (theta2 - theta1) * (cos(phi1) - cos(phi2)) 

          ! sector vector components (vector has magnitude equal to sector area)             
          s(ns, 1) = 0.5 * (phi1 - phi2 - cos(phi1) * sin(phi1) + &
               &cos(phi2) * sin(phi2)) * (sin(theta1) - sin(theta2))

          s(ns, 2) = 0.25 * (cos(theta1) - cos(theta2)) * (-2.0 * phi1 + &
               & 2.0 * phi2 + sin(2 * phi1) - sin(2 * phi2))

          s(ns, 3) = -0.25 * (theta1 - theta2) * (cos(2 * phi1) - &
               &cos(2 * phi2))

       end do

       ! update number of theta subdivisions
       if(np < nphi_middle) then
          ntheta = ntheta + 4
       else if(np > nphi_middle) then
          ntheta = ntheta - 4
       end if

    end do

    !write vectors to file for checks
    open(789, file = 's.dat', form = 'formatted')
    do ns=1,nsects
       write(789, '(3e16.8)') s(ns, :) 
    end do
    close(789, status = 'keep')

    return

  end subroutine sectors_almost_equal


  subroutine sectors_equal

    implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! perform subdivision of the unit sphere in equal-area angular control volumes !
    ! using the method of:                                                         !
    !                                                                              !
    ! Leopardi, Paul. "A partition of the unit sphere into regions of equal area   !
    ! and small diameter." Electronic Transactions on Numerical Analysis 25.12     !
    ! (2006): 309-327.                                                             !
    !                                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk) :: omega, V_R, phi_c, delta_i, n_i, delta_f
    real(rk) :: delta_theta, theta1, theta2, theta, delta_phi
    real(rk) :: phi1, phi2, phi
    integer(ik) :: n, ns
    real(rk), dimension(:), allocatable :: theta_f, yy, a
    integer(ik), dimension(:), allocatable :: m
    integer(ik) :: i, j, np
    ! total sphere area
    omega = 4.0 * PI

    ! partition (sector) area
    V_R = omega / real(nsects, 8)

    ! colattitude of north polar spherical cap
    phi_c = 2.0 * asin(sqrt(1.0_8 / real(nsects, 8)))

    ! ideal collar angle
    delta_i = sqrt(V_R)

    ! ideal number of collars
    n_i = (PI - 2.0 * phi_c) / delta_i

    ! n: actual number of collars
    if(nsects == 2) then
       n = 0
    else
       n = max(1, floor(n_i + 0.5))
    end if

    allocate(theta_f(1:n + 1))

    ! ideal number of regions for each collar
    delta_f = n_i / n * delta_i
    do i=1,n+1
       theta_f(i) = phi_c + (i - 1) * delta_f
    end do

    allocate(yy(1:n))

    do i=1,n
       yy(i) = (V(theta_f(i+1)) - V(theta_f(i))) / V_R
    end do

    ! m(1:n): actual number of regions for each collar
    allocate(m(1:n))
    allocate(a(0:n))
    a = 0.0
    do i=1,n
       m(i) = floor(yy(i) + a(i-1) + 0.5)
       a(i) = a(i-1) + yy(i) - m(i)
    end do

    deallocate(theta_f)
    deallocate(yy)
    deallocate(a)

    ! calculate direction vectors
    s(1, 1:3) = V_R * (/0.0, 0.0, 1.0/) ! north pole vector
    s(2, 1:3) = V_R * (/0.0, 0.0, -1.0/) ! south pole vector
    ! for the rest of the direction vectors, subdivide spherical
    ! coordinates theta and phi
    delta_phi = (PI - 2.0 * phi_c) / real(n, 8) ! theta increment
    ns = 3
    do i=1,n
       phi1 = phi_c + (i - 1) * delta_phi ! theta startpoint of sector
       phi2 = phi_c + i * delta_phi ! theta endpoint of sector
       delta_theta = 2.0 * PI / real(m(i), 8) ! phi increment
       do j=1,m(i)
          theta1 = (j - 1) * delta_theta ! phi startpoint of sector
          theta2 = j * delta_theta ! phi endpoint of sector

          ! sector vector components (vector has magnitude equal to sector area)             
          s(ns, 1) = 0.5 * (phi1 - phi2 - cos(phi1) * sin(phi1) + &
               &cos(phi2) * sin(phi2)) * (sin(theta1) - sin(theta2))

          s(ns, 2) = 0.25 * (cos(theta1) - cos(theta2)) * (-2.0 * phi1 + &
               & 2.0 * phi2 + sin(2 * phi1) - sin(2 * phi2))

          s(ns, 3) = -0.25 * (theta1 - theta2) * (cos(2 * phi1) - &
               &cos(2 * phi2))

          ns = ns + 1 ! next sector
       end do
    end do

    omeg(:) = V_R ! solid angles of sectors are equal

    ! wtite vectors to file for checks
    open(789, file = 's.dat', form = 'formatted')
    do ns=1,nsects
       write(789, '(3e16.8)') s(ns, :)
    end do
    close(789, status = 'keep')


    deallocate(m)

    return

  end subroutine sectors_equal

  function V(theta) 

    implicit none
    real(rk) :: V
    real(rk) :: theta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! area of spherical cap of angle theta !
    ! function (2.9) in (Leonardi, 2006)   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    V = 4.0 * PI * sin(0.5 * theta) * sin(0.5 * theta)

    return

  end function V

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

  subroutine copy(nn,dest,source,nfs,nfe,ilim)
    use types
    implicit none
    integer(ik), dimension(4), intent(in) :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(out) :: dest
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(in) :: source
    integer(ik), optional :: nfs,nfe,ilim
    !
    ! copies array source to array dest
    !
    integer(ik) :: i,j,k,l,iilim,nnfs,nnfe

    nnfs=1
    if(present(nfs)) nnfs=nfs
    nnfe=nn(4)
    if(present(nfe)) nnfe=nfe
    iilim=dim1(nn(1))
    if(present(ilim)) iilim=ilim
    !$omp parallel do
    do l=nnfs,nnfe ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,iilim
       dest(i,j,k,l)=source(i,j,k,l)
    end do; end do ; end do ; end do
    !$omp end parallel do

    return

  end subroutine copy

  subroutine alloc_init
    use mpivars
    implicit none
    !
    ! allocates and initializes all arrays
    ! 
    integer(ik) :: nmax,i,j,k,l,n
    real(rk)    :: kk
    logical     :: trunc_crit

    nmax=max(n1,n2,n3)

    ! set number of fields
    ! velocity is certain
    nu1=1
    nu2=2
    nu3=3
    nfields=3

    if(.not.PASSIVE_SCALAR) numscls = 0_ik

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

    ntemp = nsclf

    ! set vector of array dimensions
    nn(1)=n1
    nn(2)=n2
    nn(3)=n3
    nn(4)=nfields
    !$acc update device(nn(1:4), TEMPMIN, TEMPMAX)
    if(nsects == 0) nsects = nn(1)
    !$acc update device(nsects)

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

    allocate(scratch2(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    allocate(press(1:dim1(nn(1)),1:nn(2),1:nn(3)))
    !$acc enter data create(press(1:dim1(nn(1)),1:nn(2),1:nn(3)))
    press(:,:,:)=0.0_rk
    !$acc update device(press(1:dim1(nn(1)),1:nn(2),1:nn(3)))
    if(RADIATION) then
       call allocate_fvdom
       call init_fvdom
    end if

    ! allocate secondary arrays
    ! used as scratch
    allocate(scratch(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,scratch)

    ! used for comupting mean-square values in subroutine msvalue
    allocate(rmsarr(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
    call zero(nn,rmsarr)

    ! used for invariants validation only 
    if(MHD.or.VALID) then
       allocate(arr_en_1(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)))
       call zero(nn,arr_en_1)
    end if

    if(PASSIVE_SCALAR.or.VALID) then
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
    if(MHD) then
       ! ambipolar diffusion or Hall effect terms
       allocate(ad(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nfields))
       call zero(nn,ad)
    end if

    !phases array for phase-shifting
    if(DEALIASING/=NONE.or.(RADIATION.and.RADIATION_COUPLING)) then
       ! forward and inverse phase factors
       allocate(phases(1:dim1(nn(1)),1:nn(2),1:nn(3),1:3))
       allocate(iphases(1:dim1(nn(1)),1:nn(2),1:nn(3),1:3))
       do l=1,3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1
          phases(i,j,k,l) = 0.0_rk
          iphases(i,j,k,l) = 0.0_rk
       end do;  end do ;  end do ;  end do
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

    !phases array for phase-shifting
    if(DEALIASING/=NONE.or.(RADIATION.and.RADIATION_COUPLING)) then
       ! forward and inverse phase factors
       call make_phases_array
    end if

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


    if(dim==1) then
       wv=k1((i-1)/2+1)
    else if(dim==2) then
       wv=k2(j)
    else if(dim==3) then
       wv=k3(k)
    end if


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

  function tr_wv_idx(idx,n) result(tr_idx)
    integer(ik), intent(in) :: tr_idx,n
    integer(ik) :: idx
    
    if(idx>=0) then
       tr_idx=idx+1
    else if(idx<0) then
       tr_idx=n-abs(idx)
    end if
        
    return
  end function tr_wv_idx
  
  subroutine make_phases_array
    implicit none
    !
    ! prepares phases for Patterson & Orszag (1972) dealiasing
    !
    integer(ik)                                  :: i,j,k,l,n
    real(rk)                                     :: hdx
    real(rk), dimension(3)                       :: ddx, wvvec
    complex(ck) :: phase, iphase
    n=maxval(nn(1:3))
    hdx = PI / real(n, rk)
    ddx(1) = hdx
    ddx(2) = 2.0_rk * hdx
    ddx(3) = 0.5_rk * hdx

    do l=1,3
       do k=1,nn(3)
          do j=1,nn(2)
             do i=1,dim1(nn(1))
                wvvec(1) = wv(1_ik,i,j,k)
                wvvec(2) = wv(2_ik,i,j,k)
                wvvec(3) = wv(3_ik,i,j,k)
                ! forward phase
                phase=exp(ii*ddx(l)*sum(wvvec))
                ! inverse phase
                if(abs(phase) > small) then
                   iphase=1.0_rks/phase
                else
                   iphases(i,j,k,l)=0.0_rk
                end if
                phases(i,j,k,l)=phase
                iphases(i,j,k,l)=iphase
             end do
          end do
       end do
    end do

    return

  end subroutine make_phases_array



end module data
