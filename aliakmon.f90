!!$     ___ __                                       
!!$ (  / _ \\ \       /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!!$   Copyright (c) 2009-2020 George Momferatos
program aliakmon
  use types
  use parameters
  use data
  use initial_conditions 
  use numerics
  use input_output
  use validation
  use mpivars
  use hdf5_aliakmon
#ifdef _MPI_
  use mpi
  use fft_heffte
#endif
#ifdef _OPENMP_
  use omp_lib
#endif
#ifdef _CUFFT_
  use cudafor
#endif
  use iso_c_binding
  implicit none
  ! 
  ! Main program
  ! 
  integer(ik)                         :: i, j, k
  integer(ik)                         :: days,hours,mins,secs
  integer(ik)                         :: nfilespec,nfilestrfun
  real(rk)                            :: eta, tau, time,tstart
  real(rk)                            :: t1,t2,maxu
  real(rk)                            :: REtmp,et,lambdatar,emeantar
  character(256)                      :: outfile
  integer(i4b)                        :: RAND_SIZE
  integer(i4b), dimension(33)         :: rand_seed
  integer(i4b), dimension(8)          :: date_time
  real(rk), dimension(:), allocatable :: nrgs
  integer(ik)                         :: nhdf5file, nvortfile
  integer                             :: nodenamelen
  character(len=256)                  :: nodename
#ifdef _CUFFT_
  type(cudadeviceprop)                :: prop
#endif
  real(rk) :: maxc
  character(len=128)                  :: strlocrank
  integer(ik)                         :: locrank

  interface
     subroutine print_gpus() bind(C)
     end subroutine print_gpus
  end interface

  interface
     subroutine heffte_set_num_device(numd) bind(c)
       use iso_c_binding
       implicit none
       integer(c_int), value :: numd
     end subroutine heffte_set_num_device
  end interface
  
  tstart=0.0_rk
  dt=1.0e-9_rk
  RAND_SIZE=1
  MSFAC=1.0_rk
  TRFAC=sqrt(2.0_rk)/3.0_rk
  RMSUTAR=sqrt(1.0_rk/3.0_rk)
  nfilestrfun=0
  nhdf5file=0
  nvortfile=0


#ifdef _MPI_
  ! Initilize MPI enviroment
#ifdef _OPENMP_
  !$omp parallel
  nt=omp_get_num_threads()
  !$omp end parallel
  !make#else
  nt=1
#endif
  call initialize_mpi(nt)
#else
  mpirank=MPIROOT
  mpisize=1
#endif

  !call print_gpus()

  nodenamelen=256
  call mpi_get_processor_name(nodename, nodenamelen, mpierr)
#ifdef _CUFFT_
  call mpi_barrier(MPI_COMM_WORLD, mpierr)
  cudaerror = cudagetdevicecount(numdevices)
  call get_environment_variable('OMPI_COMM_WORLD_LOCAL_RANK', strlocrank)
  read(strlocrank, *) locrank
  if(numdevices > 1) then
     numdevice = locrank
  else
     numdevice = 0
  end if
  cudaerror = cudasetdevice(numdevice)
  cudaerror = cudagetdevice(numdevice)
  call heffte_set_num_device(numdevice)
!!$  print '(a,i0,2a,3(a,i0),a)', 'MPI Process ', mpirank, ' running @ node ', &
!!$       &trim(nodename), ' with localrank ', locrank, ' uses GPU # ', &
!!$       &numdevice + 1, ' out of ', &
!!$       &numdevices, ' available on the node.'
  call mpi_barrier(MPI_COMM_WORLD, mpierr)
#endif
  ! Read input file
  if(mpirank==MPIROOT) print *, 'Reading input file...'
  call read_namelist_file
  lkstart=0
  ljstart=0
#ifdef _MPI_
  ! Allocate FFT structures

  call fft_heffte_alloc(n1,n2,n3,ljsize,ljstart,lksize,lkstart)

  n2=ljsize
  n3=lksize

#else
#ifdef _OPENMP_
  call fft_omp_alloc(n1,n2,n3)
#endif
#endif

  if(command_argument_count() > 0) then
     call alloc_init_post_proc
     call post_processing_mode
#ifdef _MPI_
     call finalize_mpi
#endif
     stop
  end if

  ! Allocate and initialize all arrays


  call alloc_init



  if(mpirank == mpiroot) print *, 'Finished memory allocation and &
       &initialization.'

  !call test_fft

  MSFAC=1.0_rk/real(nn(1)*gn2*gn3,rk)

  ! Calculate viscosity
  if(RE==0.0_rk) then
     eta=KMAXETA/kmax
     REtmp=((D**(1.0_rk/4.0_rk))*eta)**(-4.0_rk/3.0_rk)
     RE=6.0_rk*sqrt(REtmp)
  else
     REtmp=(RE/6.0_rk)**2
     eta=D**(-1.0_rk/4.0_rk)*REtmp**(-3.0_rk/4.0_rk)
  end if

  visc(:)=0.0_rk
  visc(:)=15.0_rk**(1.0_rk/4.0_rk)*(eta*RMSUTAR)/(sqrt(RE))

  emeantar=visc(nu1)**3/eta**4
  lambdatar=sqrt(15.0_rk*visc(nu1)*RMSUTAR**2/emeantar)

  ! Calculate particle timescale
  if(PARTICLES.and.INERTIAL) then
     TP=STK*(lambdatar/RMSUTAR)
     if(mpirank==0) print *, 'Particle Relaxation Time: Tp=', TP
  end if

  ! FIX this: to multiple scalars
  if(.not.VISCOUS) visc(nu1:nu3)=0.0_rk
  if(PASSIVE_SCALAR) then
     if(DIFFUSIVE) then
        visc(nsclf:nscll)=visc(nsclf:nscll)/PR
     else
        visc(nsclf:nscll)=0.0_rk
     end if
  end if

  if(MHD) then
     if(RESISTIVE) then
        visc(nb1:nb3)=visc(nb1:nb3)/MAGNETIC_PR
     else
        visc(nb1:nb3)=0.0_rk
     end if
  end if

  emean=D
  emeanscl=1.0_rk
  tau=(REtmp)**(-1.0_rk/2.0_rk)

  ! Get number of OpenMP threads
  nt=1
#ifdef _OPENMP_
  !$omp parallel
  nt=omp_get_num_threads()
  !$omp end parallel
#endif


  ! Print header
  if(mpirank==MPIROOT) then
     print '(a)', "---------------------------------------------------------   "
     print '(a)', "|      ___ __                                           |   "
     print '(a)', "|  (  / _ \\\\ \\       /                                  |"
     print '(a)', "|    | |_| |\\ \\  _  __  ___  ___   _  __   __  _  __    | "
     print '(a)', "|    |  _  | > \\| |/  \/ / |/ / | | |/ / _ \\ \\| |/ /   &
          &  |"
     print '(a)', "|    | | | |/ ^ \\ ( ()  <|   <| |_| | |_/ \\_| | / /     | "
     print '(a)', "|    |_| |_/_/ \\_\\_)__/\\_\\_|\\_\\ ._,_|\\___^___/|__/ &
          &     |"
     print '(a)', "|                             |_|                       |" 
     print '(a)', "-------------------------------------------------------- "
     print '(a)', "|ALIAKMON Spectral Code for Fluid Turbulence Simulations|"
     print '(a)', "|                                                       |"  
     print '(a)', "|       (C) 2009-2020 Georgios Momferatos               |"
     print '(a)', "---------------------------------------------------------"
     write(*,*)  ' '

     write(*,'(a)') '###########################################&
          &####################################'
     write(*,'(a,i5,a,i5,a)')  '#',mpisize, ' MPI processes with', nt, &
          &' OpenMP threads each. ###########################'
     write(*,'(a)') '###########################################&
          &####################################'

     write(*,*)  ' '

     write(*,'(a)') '-------------------------------------------&
          &------------------------------------'
     write(*,'(a)') '+++++++++++++++++++++++++++++++++++++++++++&
          &++++++++++++++++++++++++++++++++++++'
     write(*,&
          &'(a,i7,3(a,f10.3),a,e10.3)')  '| N= ', max(n1,n2,n3),'| kmax=', &
          &k_max,'| REl=', RE, '| RE=', REtmp, '| VISC: ', visc(nu1)
     write(*,&
          & '(2(a,f7.2,a,e10.3))')  '| kmax*eta: ',kmax*eta,'| eta: ',eta,&
          & ' |kmax*lambda=', kmax*lambdatar ,'| lambda:', lambdatar 
     write(*,'(a)') '-------------------------------------------&
          &------------------------------------'
     write(*,'(a)') '+++++++++++++++++++++++++++++++++++++++++++&
          &++++++++++++++++++++++++++++++++++++'

     write(*,*)  ' '
  end if

  ! seed random number generator

  rand_seed(:)=1_i4b
  if(SEEDRANDOM) then
     call date_and_time(values=date_time)
     rand_seed(1:8)=date_time(:)+mpirank
  end if

  call random_seed(put=rand_seed)

  if(INPUT_FIELD) nvortfile=NFILESTART
  ! read input field or set initial conditions
  if(INPUT_FIELD) then
     if(mpirank==MPIROOT) print *, 'Reading restart file'
     call zero(nn,fu)
     call read_hdf5_file(nn,u,trim(INPUT_FIELD_FILENAME))

     call copy(nn,fu,u)

     call fourier(nn,1_ik,fu)

     if(mpirank==MPIROOT) print *, 'Restart file OK.'
  else
     if(mpirank == mpiroot) print *, 'Initial conditions...'
     call set_initial_conditions(nn,u,fu)
     if(mpirank == mpiroot) print *, 'Initial conditions are set.'
     maxc=incompressibility(nn,fu,nu1)
     if(mpirank == mpiroot) print '(a,e10.3)', 'max(comp) = ', &
          &maxc
     !stop
     nhdf5file=0_ik
     !call output_files(0_ik)
     nhdf5file=nhdf5file+1_ik

     nvortfile=0
     call output_slices(nvortfile, time)
     nvortfile=nvortfile+1_ik

  end if

  allocate(nrgs(1:nn(4)))
  ! Kinetic energy
  call msvalue(nn,fu,nrgs)
  KE=nrgs(nu1)


  ! Magnetic energy
  ME=0.0_rk
  if(MHD) then
     ME=nrgs(nb1)
  else
     ME=0.0_rk
  end if
  teinitial=0.5_rk*(ME+KE)
  rkten=teinitial

  ! Mean kinetic helicity
  mkhinitial=mean_kinetic_helicity(nn,fu)
  rkmkh=mkhinitial
  ! MHD helicities
  if(MHD) then
     mchinitial=mean_cross_helicity(nn,fu)
     mmhinitial=mean_magnetic_helicity(nn,fu)
  end if

  ! set up .dat files for output
  if(mpirank==MPIROOT) then
     call setup_output_dat_files
  end if

  ! Start counting time
  call timing(t1)
  time=0.0_rk
  k=1
  totdis=0.0_rk
  totdisprev=0.0_rk
  totdisprev2=0.0_rk
  if(mpirank == mpiroot) print *, 'Entering main loop...'
  ! main time loop
  timeloop:do while((TIMESTEPS==0.and.time-tstart<=TMAX).or.&
       &(TIMESTEPS.ne.0.and.k<=TIMESTEPS))
     
     ! Print progress to stdout
     call print_progress(k,time,&
          &tstart,t1)

     ! Calculate dissipation rate
     totdis=emean
     if(MHD) then
        totdis=totdis+emeanb
        if(AMB_DIFF) then
           addis=ambipolar_diffusion_dissipation(nn,fu)
           ! Uncomment the following line to include AD dissipation
           !totdis=totdis+addis
        end if
     end if

     ! Broadcast dissipation rates
#ifdef _MPI_
     sbuf(1)=totdis
     sbuf(2)=totdisprev
     sbuf(3)=totdisprev2
     call mpi_bcast(sbuf,3,MPIRK,MPIROOT,MPI_COMM_WORLD,mpierr)
     totdis=sbuf(1)
     totdisprev=sbuf(2)
     totdisprev2=sbuf(3)
#endif

     ! Check for temporal total dissipation peak
     if(totdis-totdisprev<0.0_rk.and.k>3.and.&
          &totdisprev-totdisprev2>0.0_rk.and..not.FORCED) then
        ! Output files
        call reached_dissipation_peak
        if(STOP_AT_DISSPEAK) exit timeloop
     end if
     totdisprev2=totdisprev
     totdisprev=totdis

     call cfl_condition(nn,dt)
     ! Advance in time
     call timestep(nn,fu,dt)
     
     !call energy_test(nn,u,fu,rhs)     
     ! Write maxima to file
     write(maxima_dat,'(6e17.8)') time,maxu,maxb,MAXVORT, MAXJ, MAXLF
     if(PARTICLES) then
        ldt=dt/NPTS
        do i=1,NPTS
           !call ltimestep(np,x,vp,ldt)
        end do
     end if

     time=time+dt     
     k=k+1

     if(int(floor(time * slicefrate), ik) == nvortfile) then
        call output_slices(nvortfile, time)
        call output_spectra(0_ik)
        nvortfile = nvortfile + 1
     end if

     if(int(floor(time * hdf5frate), ik) == &
          &nhdf5file .and. hdf5frate /= 0.0_rk) then
        call output_files(nhdf5file)
        nhdf5file = nhdf5file + 1
     end if
          
  end do timeloop

  ! Print estimated total time information
  if(mpirank==MPIROOT) then
     call timing(t2) 
     et=(t2-t1)
     days=int(et/(60.*60.*24.))
     et=et-days*(60.*60.*24.)
     hours=int(et/(60.*60.))
     et=et-hours*(60.*60.)
     mins=int(et/60.)
     secs=et-mins*60.

     print '(4(a,i4),2(a,f10.3))', 'Time: ', days,':',hours,':',mins,':', &
          &secs, ' Time/timestep: ', ((t2-t1))/(k-1), ' timesteps/hour: ',&
          & k/((t2-t1))*60*60


     close(maxima_dat,status='keep')
  end if
  call output_files(888888_ik)
  call output_spectra(888888_ik)

#ifdef _MPI_
  call fft_heffte_dealloc
  call finalize_mpi
#endif

  stop 

contains

  subroutine reached_dissipation_peak
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Outputs files on the dissipation peak
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(mpirank==MPIROOT) print *, 'Dissipation Peak'
    if(.not.FORCED.and..not.FORCED_MHD) then
       TMAX=time+5.0*ETT
    end if

    nfilestrfun=nfilestrfun+1
    if(mpirank==MPIROOT) print *, 'Writing files'

    !Write output file
    call output_files(999999_ik)

    !Output information about the dissipation peak
    open(newunit=dispeak_dat,file='dispeak.dat',form='formatted')
    write(dispeak_dat,'(11e17.8)') time,KE,mkin_hel,emean,ils,lambda,eta,&
         &rmsu*ils/visc(nu1),REl,mkhdis,fscale(nu1)
    if(MHD) write(987,'(9e17.8)') time,emeanb,mcross_hel,ME,mmh,&
         &AD_COEFF*mlor,mmhdis,mchdis,fscale(nb1)
    close(dispeak_dat,status='keep')


    call output_spectra(99999_ik)
    nfilespec=nfilespec+1

    return

  end subroutine reached_dissipation_peak

  subroutine output_spectra(nfilespec)
    implicit none
    integer(ik) :: nfilespec
!!!!!!!!!!!!!!!!!!!!!!!
    !Outputs energy spectra
!!!!!!!!!!!!!!!!!!!!!!!

    write(outfile,'(a,i0.5,a)') 'espec.', nfilespec,'.dat'
    call vector_spectrum(nn,fu,nu1,outfile,nespec)
    if(PASSIVE_SCALAR) then
       write(outfile,'(a,i0.5,a)') 'psspec.', nfilespec,'.dat'
       call scalar_spectrum(nn,fu,nscl,outfile,nespec)
    end if
    if(MHD) then
       write(outfile,'(a,i0.5,a)') 'bespec.', nfilespec,'.dat'
       call vector_spectrum(nn,fu,nb1,outfile,nespec)
    end if

    return

  end subroutine output_spectra


  subroutine output_files(num)
    use hdf5_aliakmon
    use numerics, only: calcqr
    implicit none
    integer(ik), intent(IN) :: num
    !Output fields in files
    
    if(RADIATION) then
       call calcqr
    end if
    
    call copy(nn,u,fu)

    call fourier(nn,-1_ik,u)
    call write_hdf5_file(nn,u,time,num)

    call zero(nn,u)

    return

  end subroutine output_files

  subroutine setup_output_dat_files
    if(INPUT_FIELD) then
       open(newunit=maxima_dat,file='maxima.dat',form='formatted',&
            &action='write',&
            &status='old',position='append')
    else
       open(maxima_dat,file='maxima.dat',form='formatted',action='write')
    end if
    write(maxima_dat,'(6a17)') 'Time|','max(u)|','max(b)|','max(w)|', &
         &'max(j)|', 'max(j x b)|'
    write(maxima_dat,'(6a17)') 't|','maxu|','maxb|','maxw|', &
         &'maxj|', 'maxjxb|'
    if(INPUT_FIELD) then
       open(newunit=hydro_dat,file='hydro.dat',form='formatted',action='write',&
            &status='old',position='append')
    else
       open(newunit=hydro_dat,file='hydro.dat',form='formatted',action='write')
       write(hydro_dat,'(11a37)') 'Time|','Kinetic energy|',&
            &'Mean kinetic helicity|','Mean energy dissipation|',&
            &'Integral length scale|','Taylor microscale|',&
            &'Kolmogorov microscale|','Reynolds number|',&
            &'Taylor microscale Reynolds number|',&
            &'Mean kinetic helicity dissipation|','Forcing scale'
       write(hydro_dat,'(11a)') 't|','ke|', 'mkh|','e|',&
            &'ils|','lambda|', 'eta|','Re|', 'Rel|', 'mkhdiss|','fscale'
    end if
    if(INPUT_FIELD.and.PASSIVE_SCALAR) then
       open(newunit=passive_dat,file='passive.dat',form='formatted',&
            &action='write',&
            &status='old',position='append')
    else if(PASSIVE_SCALAR) then
       open(newunit=passive_dat,file='passive.dat',form='formatted',&
            &action='write')
       write(passive_dat,'(11a41)') 'Time|','Mean scalar value|',&
            &'Scalar variance|','Maximum scalar value|',&
            &'Minimum scalar value|','Mean scalar dissipation|',&
            &'Maximum scalar gradient|',&
            &'Passive scalar microscale Peclet number|',&
            &'Passive scalar microscale|',&
            &'Obukhov-Corsin microscale|',&
            &'Scalar forcing scale|'
       write(passive_dat,'(11a37)') 't|','mscl|', 'sclvar|','maxscl|',&
            & 'minscl|','scldiss|', 'maxsclgrad|', 'Pel|', 'scllambda|',&
            &'etaoc|', 'sclfscale|'
    end if
    if(INPUT_FIELD.and.MHD) then
       open(magnetic_dat,file='magnetic.dat',form='formatted',action='write',&
            &status='old',position='append')
    else
       open(newunit=magnetic_dat,file='magnetic.dat',form='formatted',&
            &action='write')
       write(magnetic_dat,'(7a37)') 'time|','Magnetic energy|',&
            &'Ohmic dissipation|','Mean magnetic helicity|',&
            &'Mean magnetic helicity dissipation|',&
            &'Mean cross-helicity|','Mean cross-helicity dissipation|'
       write(magnetic_dat,'(7a37)') 't|','me|','odiss|','mmh|','mmhdiss|',&
            &'mch|','mchdiss|'
    end if
    if(INPUT_FIELD) then
       open(newunit=distest_dat,file='distest.dat',form='formatted',&
            &action='write',&
            &status='old',position='append')
    else
       open(newunit=distest_dat,file='distest.dat',form='formatted',&
            &action='write')
    end if

    return

  end subroutine setup_output_dat_files

  subroutine post_processing_mode
    use hdf5_aliakmon,only:read_hdf5_file,write_hdf5_file
    implicit none
    integer                         :: arg_num
    character(len=256)              :: arg
    character(len=256)              :: filename_in
    character(len=256)              :: filename_out

    if(command_argument_count()==1) then

       if(mpirank == MPIROOT) then
          print '(a)',  ''
          print '(a)',  'usage: aliakmon [OPTIONS] filename.h5'
          print '(a)',  ''
          print '(a)',  'Without further options, read aliakmon.nml and &
               &compute.'
          print '(a)',  ''
          print '(a)',  'post-processing options:'
          print '(a)',  ''
          print '(a)',  '  -v, --vorticity'
          print '(a)',  '  -e, --dissipation'
          print '(a)',  '  -s, --scalar_dissipation'
          print '(a)',  '  -j, --current'
          print '(a)',  '  -a, --all'
       end if

       return

    end if

    !parse command line arguments
    do arg_num=1,command_argument_count()
       call get_command_argument(arg_num, arg)
       if(arg(1:1)=='-') then
          select case(arg)
          case('-v','--vorticity')
             OUTPUT_W=.true.
          case('-e','--dissipation')
             OUTPUT_DISS=.true.
          case('-s','--scalar-dissipation')
             OUTPUT_SCL_DISS=.true.
          case('-j','--current')
             OUTPUT_J=.true.
          case('-a','--all')
             OUTPUT_W=.true.
             OUTPUT_DISS=.true.
             OUTPUT_SCL_DISS=.true.
             OUTPUT_J=.true.
          case default
             if(mpirank==mpiroot) then
                print '(2a)', 'Unrecognized command-line option: ', trim(arg)
                print '(a)',  'usage: aliakmon [OPTIONS] filename.h5'
                print '(a)',  ''
                print '(a)',  'Without further options, read aliakmon.nml and &
                     &compute.'
                print '(a)',  ''
                print '(a)',  'post-processing options:'
                print '(a)',  ''
                print '(a)',  '  -v, --vorticity'
                print '(a)',  '  -e, --dissipation'
                print '(a)',  '  -s, --scalar_dissipation'
                print '(a)',  '  -j, --current'
                print '(a)',  '  -a, --all'
             end if
          end select
       else
          !input HDF5 file name
          filename_in=trim(arg)
       end if
    end do

    !read input HDF5 file
    call read_hdf5_file(nn,u,trim(filename_in))

    !copy velocity
    call copy(nn, fu, u)
    call fourier(nn, 1_ik, fu)
    !calculate vorticity
    if(OUTPUT_W) then
       call curl(nn, u, fu, nu1)
       call fourier(nn, -1_ik, u, nu1, nu3)
    end if

    if(OUTPUT_DISS) then
       ! Calculate viscosity
       if(RE==0.0_rk) then
          eta=KMAXETA/kmax
          REtmp=((D**(1.0_rk/4.0_rk))*eta)**(-4.0_rk/3.0_rk)
          RE=6.0_rk*sqrt(REtmp)
       else
          REtmp=(RE/6.0_rk)**2
          eta=D**(-1.0_rk/4.0_rk)*REtmp**(-3.0_rk/4.0_rk)
       end if

       visc(:)=0.0_rk
       visc(:)=15.0_rk**(1.0_rk/4.0_rk)*(eta*RMSUTAR)/(sqrt(RE))

       call dissipation(nn, scratch, 1_ik, fu)

    end if
!!$
!!$    if(PASSIVE_SCALAR.and.OUTPUT_SCL_DISS) then
!!$       call gradient(nn,u,fscl,.true.)
!!$       call fourier(nn,-1_ik,u1)
!!$       call fourier(nn,-1_ik,u2)
!!$       call fourier(nn,-1_ik,u3)
!!$
!!$       rhs5=u1**2+u2**2+u3**2
!!$
!!$    end if
!!$
!!$    !calculate vorticity
!!$    if(MHD.and.OUTPUT_J) then
!!$       call curl(nn,u,fb,.true.)
!!$       call fourier(nn,-1_ik,u1)
!!$       call fourier(nn,-1_ik,u2)
!!$       call fourier(nn,-1_ik,u3)   
!!$    end if

    !create output HDF5 filename
    write(filename_out,'(2a)') 'pp_',trim(filename_in)

    !write output HDF5 file
    call write_hdf5_file(nn,u,scratch,trim(filename_out),time)

    return

  end subroutine post_processing_mode

  subroutine output_slices(nfile, time)
    use mpivars
    use hdf5_aliakmon
    implicit none
    integer(ik), intent(IN) :: nfile
    real(rk), intent(IN)    :: time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=1024)     :: fname, comment
    integer(ik) :: l, m, nj, nf
    character(len=64), dimension(:), allocatable :: datanames
    
    nfields = 2
    if(PASSIVE_SCALAR) nfields = nfields + 2 * numscls
    if(MHD) nfields = nfields + 1

    allocate(datanames(1:nfields))
    if(.not.allocated(slice)) allocate(slice(1:nfields,1:nn(2),1:nn(3)))

    datanames(nu1)='w'
    datanames(nu2)='e'
    if(PASSIVE_SCALAR) then
       do l=3,3+numscls-1
          write(datanames(l),'(a,i3.3)') 'scl_', l - 3 + 1
       end do
       do l=3+numscls,3+numscls+numscls-1
          m = l - (3 + numscls) + 1 
          write(datanames(l),'(a,i3.3)') 'sg_', m
       end do
    end if

    if(MHD) then
       datanames(nfields) = 'j'
    end if



    call dissipation(nn,scratch,nu1,fu)
    call curl(nn,rmsarr,fu,nu1)
    call fourier(nn,-1_ik,rmsarr,nfs=nu1,nfe=nu3)
    call copy(nn,u,fu)
    call fourier(nn,-1_ik,u)
    if(PASSIVE_SCALAR) then
       call gradient(nn,fsclgrads,fu,nsclf,nscll)
       do l=1,3
          fsclgrad=>fsclgrads(:,:,:,:,l)
          call fourier(nn,-1_ik,fsclgrad)
       end do
    end if
    write(fname,'(a,i5.5,a)') 'slice-', nfile, '.h5'

    if(MHD) then
       call curl(nn,rks1,fu,nb1)
       call fourier(nn,-1_ik,rks1,nfs=nb1,nfe=nb3)
    end if
    
    !$omp parallel do
    do l=1,nfields ; do k=1,nn(3) ; do j=1,nn(2) 
       slice(l,j,k) = 0.0
    end do; end do ; end do
    !$omp end parallel do
    
    i=1
    !$omp parallel do
    do k=1,nn(3) ; do j=1,nn(2) 
       slice(1,j,k) = sqrt(rmsarr(i,j,k,nu1)**2 + &
            &rmsarr(i,j,k,nu2)**2 + rmsarr(i,j,k,nu3)**2)
       slice(2,j,k) = scratch(i,j,k,nu1)
    end do; end do
    !$omp end parallel do
    if(PASSIVE_SCALAR) then
       do l=3,3+numscls-1
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) 
             slice(l,j,k) = u(i,j,k,nsclf+l-3)
          end do; end do
          !$omp end parallel do
       end do
       do l=3+numscls,3+numscls+numscls-1
          m = l - (3 + numscls) + nsclf
          !$omp parallel do 
          do k=1,nn(3) ; do j=1,nn(2)
             slice(l,j,k) = sqrt(fsclgrads(i,j,k,m,1)**2+&
                  &fsclgrads(i,j,k,m,1)**2+fsclgrads(i,j,k,m,1)**2)
          end do; end do
          !$omp end parallel do
       end do
    end if

    if(PASSIVE_SCALAR) then
       nj = 3+numscls+numscls
    else
       nj = 3
    end if

    if(MHD) then
       !$omp parallel do 
       do k=1,nn(3) ; do j=1,nn(2)
          slice(nj,j,k)=sqrt(rks1(i,j,k,nb1)**2+rks1(i,j,k,nb2)**2+&
               &rks1(i,j,k,nb3)**2)
       end do; end do
       !$omp end parallel do
    end if

    call write_hdf5_file(nn,nfields,slice,datanames,time,nfile)

    call zero(nn,rmsarr)
    call zero(nn,scratch)
    if(PASSIVE_SCALAR) call zero(nn,fsclgrads)
    if(MHD) call zero(nn,rks1)

    deallocate(datanames)

    return
  end subroutine output_slices

  subroutine test_fft
    real(rk) :: maxerr
    integer(ik) :: i, l

    do l=1,5
       call random_number(u(1:nn(1),1:nn(2),1:nn(3),1:nn(4)))
       fu=u

       call fourier(nn,1_ik,fu)
       call fourier(nn,-1_ik,fu,trunc=.false.)


       maxerr=maxval(abs(u(1:nn(1),1:nn(2),1:nn(3),1:nn(4))-&
            &fu(1:nn(1),1:nn(2),1:nn(3),1:nn(4))))
#ifdef _MPI_
       sbuf(1)=maxerr
       call mpi_allreduce(sbuf,rbuf,1,MPIRK,MPI_MAX,MPI_COMM_WORLD,mpierr)
       maxerr=rbuf(1)
#endif
       if(mpirank == mpiroot) print '(a,e10.3)', 'FFT max error:', maxerr
    end do

    call finalize_mpi

    stop

  end subroutine test_fft

end program aliakmon




