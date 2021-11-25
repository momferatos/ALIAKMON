q!!$     ___ __                                       
!!$ (  / _ \\ \        /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!$Copyright (c) 2009-2020 Georgios Momferatos
module types
  use iso_fortran_env
  implicit none
  !
  ! Integer, real and complex types
  ! 
  ! Integer types
  integer, parameter     :: i1b=1 ! integer 1 bit
  integer, parameter     :: i2b=2 ! integer 2 bits
  integer, parameter     :: i4b=4 ! integer 4 bits
  integer, parameter     :: i8b=8 ! integer 8 bits
  integer, parameter     :: ik=i8b! main integer type
  ! Real types
  integer, parameter     :: sp=4  ! single precision real
  integer, parameter     :: dp=8  ! double precision real
  integer, parameter     :: qp=real128 ! quadruple precision real
  integer, parameter     :: rk=dp ! main real type
  ! Real type for arrays
#ifdef _DOUBLE_
  integer, parameter     :: rks=dp
#else
  integer, parameter     :: rks=sp
#endif
  ! Complex types
  integer, parameter     :: csp=4  ! single precision complex
  integer, parameter     :: cdp=8  ! double precision real
  integer, parameter     :: ck=cdp ! main complex type
  integer, parameter     :: cks=rks! complex type for arrays
  ! Real types for I/O
  integer(ik), parameter :: inrk=sp
  integer(ik), parameter :: outrk=dp
end module types

module parameters
  use types
  use iso_c_binding
  implicit none
  ! 
  ! Main parameter module
  ! 
  ! number of fields
  integer(c_int), bind(C)     :: cudaerror, numdevices, numdevice
  integer(ik)                                   :: nfields
  integer(ik)                                   :: nu1
  integer(ik)                                   :: nu2
  integer(ik)                                   :: nu3
  integer(ik)                                   :: nscl
  integer(ik)                                   :: nsclf
  integer(ik)                                   :: nscll
  integer(ik)                                   :: nb1
  integer(ik)                                   :: nb2
  integer(ik)                                   :: nb3
  integer(ik)                                   :: ntemp
  ! For the energy conservation tests
  real(rk)                                      :: teinitial
  real(rk)                                      :: mchinitial
  real(rk)                                      :: mmhinitial
  real(rk)                                      :: mkhinitial
  real(rk)                                      :: rkvdis
  real(rk)                                      :: rkodis
  real(rk)                                      :: rkaddis
  real(rk)                                      :: rkten
  real(rk)                                      :: rkeinj
  real(rk)                                      :: rkmkh
  real(rk)                                      :: rkmkhdis
  ! Output files
  integer(ik)                                   :: dispeak_dat=111
  integer(ik)                                   :: hydro_dat=222
  integer(ik)                                   :: passive_dat=333
  integer(ik)                                   :: ank_stdout=444
  integer(ik)                                   :: maxima_dat=555
  integer(ik)                                   :: magnetic_dat=666
  integer(ik)                                   :: distest_dat=777
  real(rk), parameter    :: small=1.0e-16
  ! PI
  real(rk), parameter    :: PI=3.14159265358979323846264338327950&
       &2884197169399375105820974944592307816406286208998628034825&
       &342117068_rk
  real(rk), parameter    :: CP=4.0e3_rk
  real(rk), parameter    :: PRESS=1.01325e5
  real(rk)               :: TEMPMIN
  !$acc declare create(TEMPMIN)
  real(rk)               :: TEMPMAX
  !$acc declare create(TEMPMAX)
  ! factor for mean-square value calculation in Fourier space
  real(rk)               :: MSFAC
  ! Factor for spherical truncation
  real(rk)               :: TRFAC
  ! Length of the simulation box
  real(rk), parameter    :: LBOX=2.0_rk*PI
  ! Integral length scale
  real(rk)               :: ILS
  ! Taylor microscale, root-mean-square velocity, Taylor microscale Reynolds
  ! number
  real(rk)               :: lambda,rmsu,REl
  ! Target root-mean-square velocity
  real(rk)               :: RMSUTAR
  ! Target kinetic energy
  real(rk), parameter    :: KENTAR=0.5
  ! Target magnetic energy
  real(rk), parameter    :: MENTAR=0.5
  ! Eddy turnover time
  real(rk)               :: ETT
  ! Parameter for the estimation of target Reynolds number
  real(rk), parameter    :: D=0.4_rk
  ! Imaginary unity
  complex(ck), parameter :: ii=cmplx(0.0_rk,1.0_rk,ck)
  ! Factor for rescaling the viscosity
  ! Particle relaxation time
  real(rk)               :: TP
  ! Kinetic energy
  real(rk)               :: KE
  ! Magnetic energy
  real(rk)               :: ME
  ! Passive scalar variance
  real(rk)               :: PSCV
  ! Maximum wave-number
  real(rk)               :: kmax
  ! Number of OpenMP threads
  integer                :: nt=1
  ! Local dimensions of the domain
  integer(ik)            :: n1,n2,n3
  ! Global y dimension across MPI processes
  integer(ik)            :: gn2
  ! Global z dimension across MPI processes
  integer(ik)            :: gn3
  ! Grid Subdivision
  ! Slabs
  integer(ik), parameter :: SLABS=0
  ! Pencils
  integer(ik), parameter :: PENCILS=1
  integer(ik)            :: FFT_DECOMPOSITION=SLABS 
  ! Maximum dimension
  integer(ik)            :: nmax
  integer(ik), parameter :: nmaxsqrt=128
  ! Number of active and truncated modes
  integer(i8b)           :: nmodes, ntrunc
  ! Mean energy dissipation, mean passive scalar dissipation, mean magnetic
  ! energy dissipation, mean cross-helicity, mean magnetic helicity,
  ! maximum of the vorticity
  real(rk)               :: emean, emeanb, mcross_hel, mmh,wmax
  real(rk), dimension(:), allocatable :: emeanscl
  real(rk), dimension(:), allocatable :: sclvarprev
  ! Mean Lorentz force, mean kinetic helicity, mean magnetic helicity dissipation
  ! ,mean cross helicity dissipation, mean kinetic helicity dissipation
  real(rk)               :: mlor,mkin_hel,mmhdis,mchdis,mkhdis,gdedt
  ! Total dissipations
  real(rk)               :: totdis,totdisprev
  real(rk)               :: totdisprev2
  ! Ambipolar diffusion dissipation
  real(rk)               :: addis
  ! Kaneda et al. (2004) forcing
  logical                :: VARIABLE_FORCING=.true.
  ! For phase-shifting
  ! Constants for choice of initial conditions
  integer(ik), parameter :: ZERO_INITCOND=0,STOCHASTIC_INITCOND_FLAT=1
  integer(ik), parameter :: STOCHASTIC_INITCOND_WITH_SPECTRUM=2
  integer(ik), parameter :: ORSZAG_TANG_VORTEX=3, ABC=4, TAYLOR_GREEN_VORTEX=5, RADSPHERE=6
  ! Constants for choice of time integration method
  integer(ik), parameter :: MEULER=0, MRUNGE_KUTTA2=1,MRUNGE_KUTTA4=2
  ! Constants for choice of truncation method
  integer(ik), parameter :: TWO_THIRDS=0, SPHERICAL=1,POLYHEDRAL=2
  ! Constants for choice of dealiasing method
  integer(ik), parameter :: NONE=0,PATTERSON_ORSZAG=1
  ! Constants for choice of particle  initial conditions 
  integer(ik), parameter :: HOMOGENEOUS=0,SPHERE=1,SHEET=2
  ! Random number scaling factor
  real(rk)               :: RSCALE
  ! Maximum time
  real(rk)               :: TMAX
  ! Stop at the dissipation peak?
  logical                :: STOP_AT_DISSPEAK=.true.
  ! Maximum timesteps
  integer(ik)            :: TIMESTEPS
  ! Courant-Friedrichs-Lewy (CFL) condition
  real(rk)               :: CFL
  ! Timestepx
  real(rk)               :: dt
  ! Timestep for the particles
  real(rk)               :: ldt
  ! Logical variabkle for viscosity
  logical                :: VISCOUS
  ! Reynolds number
  real(rk)               :: RE
  ! Maximum wavenumber * Kolmogorov microscale
  real(rk)               :: KMAXETA
  ! Viscosity
  real(rk), dimension(:), allocatable  :: VISC
  ! Logical variabkle for Burgers equation
  logical                :: BURGERS
  ! Logical variabkle for forcing
  logical                :: FORCED
  ! Integer variable for choice of forcing
  integer(ik)            :: FORCING
  ! Constants for choice of forcing
  integer(ik), parameter :: STOCHASTIC=0,KANEDA=1
  ! Forcing wavenumber
  real(rk)               :: KFORCING=2
  ! Scale factor of the forcing term
  real(rk), dimension(:), allocatable    :: fscale
  ! Same for MHD
  real(rk)               :: FSCALEMHD
  ! Same for passive scalar
  logical                :: PASSIVE_SCALAR
  ! Number of passive scalars
  integer(ik)            :: numscls=0
  ! Logical variable for passive scalar heating
  logical                :: HEATING
  ! Logical variable for diffusivity
  logical                :: DIFFUSIVE
  ! Prandtl number
  real(rk)               :: PR=1.0_rk
  ! Logical variable for passive scalar forcing
  logical                :: FORCED_PASSIVE_SCALAR=.true.
  ! Scale factor of the forcing term of the passive scalar
  real(rk)               :: FSCALESCL=1.0_rk
  ! Logical variable for magnetohydrodynamics
  logical                :: MHD
  ! Beta MHD factor: <u**2>/<b**>
  real(rk)               :: BETA
  ! Magnetic Prandtl number
  real(rk)               :: MAGNETIC_PR
  ! Logical variable for resistivity
  logical                :: RESISTIVE
  ! Logical variable for ambipolar diffusion
  logical                :: AMB_DIFF
  ! Ambipolar diffusion coefficient
  real(rk)               :: AD_COEFF
  ! Logical variable for Hall effect
  logical                :: HALL
  ! Hall coefficient
  real(rk)               :: HALL_COEFF
  ! Logical variable for MHD forcing
  logical                :: FORCED_MHD
  ! Logical variable for radiation
  logical                :: RADIATION
  logical                :: RADIATION_COUPLING
  integer(ik)            :: EQSECTS
  integer(ik)            :: nsects
  !$acc declare create(nsects)
  integer(ik)            :: NITERDO
  real(rk)               :: FVTOL
  real(rk), parameter    :: STEFB=5.67037321e-8_rk
  ! Logical variable for particles
  logical                :: PARTICLES
  ! Default value for particle initial condition
  integer                :: PART_INITCOND
  ! Number of particles
  integer(ik)            :: NPART
  ! Number of particle timesteps per flow timestep
  integer(ik)            :: NPTS=3
  ! Logical variable for inertial particles
  logical                :: INERTIAL
  ! Logical variable for particle periodicity
  logical                :: PERIODIC_PARTICLES
  ! Logical variable for passive scalar mixing
  logical                :: MIXING
  ! Stokes number for inertial particles
  real(rk)               :: STK
  ! Output lagrangian hisotry
  logical                :: LAGRANGIAN_HISTORY
  ! Default integrationn method
  integer(ik)            :: INTEGRATION_METHOD
  ! Deafault truncation
  integer(ik)            :: TRUNCATION
  ! Default dealiasing
  integer(ik)            :: DEALIASING
  ! Crank-Nicholson method for dissipative terms
  logical                :: CRANK_NICHOLSON
  ! Default initial condition
  integer                :: INITCOND
  ! Wave-number for the initial conditions
  real(rk)               :: KINITCOND
  ! Seed random number generator
  logical                :: SEEDRANDOM
  ! Number of output files
  logical                :: OUTPUTFILES
  ! Read an input field as initial condition
  logical                :: INPUT_FIELD
  ! Filename of the input field
  character(128)         :: INPUT_FIELD_FILENAME
  ! File number to start file output
  integer(ik)            :: NFILESTART=0
  ! Frame rate of HDF5 file output
  real(rk)               :: hdf5frate
  ! Frame rate of VTK slice output
  real(rk)               :: slicefrate
  ! HDF5 gzip Compression level
  integer                :: COMPRESSION_LEVEL=0
  ! Maximum vorticity, maximum current density, maximum Lorentz force
  real(rk)               :: MAXVORT, MAXJ, MAXLF, MAXVEL, MAXB
  ! Maximum passive scalar gradient
  real(rk)               :: MAXGRADTHETA=0.0_rk
  ! 
  ! Variables for post-processing mode
  ! 
  logical                :: OUTPUT_W
  logical                :: OUTPUT_DISS
  logical                :: OUTPUT_SCL_DISS
  logical                :: OUTPUT_J

contains
  pure elemental function dim1(n)
    integer(ik) :: dim1
    integer(ik), intent(IN) :: n
    ! 
    ! Used to handle x dimension of the FFTW arrays
    ! 

    dim1=int(2_ik*(floor(real(n,rk)/2.0_rk)+1),ik)

    return
  end function dim1

  function round(x)
    integer(ik) :: round
    real(rk) :: x
    ! 
    ! Round a real number
    ! 
    integer(ik) :: floor
    real(rk) :: frac

    floor=int(x,ik)
    frac=x-floor
    if(frac>0.5_rk) then
       round=floor+1
    else
       round=floor
    end if

    return

  end function round

end module parameters

module mpivars
  use iso_c_binding
  use types
  use parameters, only: dim1
  ! 
  ! MPI variable and constants
  ! 
#ifdef _MPI_
  use mpi
#endif
  implicit none
  ! MPI rank
  integer(i4b), bind(C) :: mpirank=0
  ! MPI root process
  integer(i4b), parameter :: MPIROOT=0
  ! Number of MPI processes
  integer(i4b), bind(C)  :: mpisize=0
  ! Size of MPI slice across y-dimension of FFTW array
  integer(ik) :: ljsize = 0
  ! X-dimension index where FFTW array starts
  integer(ik) :: ljstart = 1
  ! Size of MPI slice across z-dimension of FFTW array
  integer(ik) :: lksize=0
  !$acc declare create(lksize)
  ! X-dimension index where FFTW array starts
  integer(ik) :: lkstart=1
  !$acc declare create(lkstart)
#ifdef _MPI_
  integer(i4b) :: mpirequest
  integer(i4b) :: mpistat(MPI_STATUS_SIZE)
  integer(i4b) :: mpierr
  ! Global value for maximum wavenumber 
  integer(ik) :: glkmax,kstep
  ! Constants for MPI data types
  integer :: MPIRK, MPI2RK, MPISP,MPICKS, MPIIK,MPIRKS
  ! Buffers used in reductions and Broadcasts
  real(dp), dimension(1024) :: sbuf,rbuf
  ! Maximum integer
  integer, parameter :: MPIMAXINT=2**16
contains
  subroutine set_mpi_types
    implicit none
    ! 
    ! Set-up MPI data types
    ! 

    ! Integer
    MPIIK=MPI_INTEGER

    ! Real
    if(rk==sp) then
       MPIRK=MPI_REAL
       MPI2RK=MPI_2REAL
    else if(rk==dp) then
       MPIRK=MPI_DOUBLE_PRECISION
       MPI2RK=MPI_2DOUBLE_PRECISION
    end if

    if(rks==sp) then
       MPIRKS=MPI_REAL
    else if(rk==dp) then
       MPIRKS=MPI_DOUBLE_PRECISION
    end if

    MPISP=MPI_REAL

    ! Complex
    if(cks==sp) then
       MPICKS=MPI_COMPLEX
    else if(cks==dp) then
       MPICKS=MPI_DOUBLE_COMPLEX
    end if

    return

  end subroutine set_mpi_types

  subroutine initialize_mpi(nt)
    implicit none
    integer(i4b), intent(in) :: nt
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(i4b) :: provided
    ! 
    ! Initialize MPI environment
    ! 

    ! Set data types
    call set_mpi_types

    ! Initialize
    if(nt==1) then
       ! Single-threaded
       call mpi_init(mpierr)
    else
       ! Multi-threaded
       call mpi_init_thread(MPI_THREAD_MULTIPLE,provided,mpierr)
       if(provided < MPI_THREAD_MULTIPLE) then
          print *, 'error: MPI failed to provide multi-thread support.'
          stop
       end if
    end if
    ! Get number of processes
    call mpi_comm_size(MPI_COMM_WORLD,mpisize,mpierr)
    ! Get process rank
    call mpi_comm_rank(MPI_COMM_WORLD, mpirank,mpierr)
    
    return
    
  end subroutine initialize_mpi

  subroutine finalize_mpi
    implicit none
    ! 
    ! Clean-up MPI environment
    ! 

    ! Finalize
    call mpi_finalize(mpierr)

    return
    
  end subroutine finalize_mpi
  
#endif
end module mpivars
