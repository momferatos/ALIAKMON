!!$     ___ __                                       
!!$ (  / _ \\ \        /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!$Copyright (c) 2009-2020 Georgios Momferatos
module input_output
  use types
  use parameters
  implicit none
!!!!!!!!!!!
  !I/O Module
!!!!!!!!!!!

contains

  subroutine read_namelist_file
    implicit none
    integer(ik)                  :: n
    integer(ik)                  :: aliakmon_nml
    
    !Read the input namelist file
    namelist /general/           n,TIMESTEPS, TMAX, CFL, RE, KMAXETA,&
         &STOP_AT_DISSPEAK
    namelist /hydro/             VISCOUS,BURGERS
    namelist /force/           FORCED, VARIABLE_FORCING,KFORCING
    namelist /passivescalar/     PASSIVE_SCALAR, DIFFUSIVE, PR, HEATING,&
         &FORCED_PASSIVE_SCALAR,NUMSCLS
    namelist /radiation/         RADIATION, EQSECTS, TEMPMIN, TEMPMAX, &
         &NSECTS, FVTOL, NITERDO
    namelist /magnetohydro/      MHD, MAGNETIC_PR,RESISTIVE,FORCED_MHD,BETA,&
         &AMB_DIFF, AD_COEFF,HALL, HALL_COEFF
    namelist /particle/          PARTICLES, PART_INITCOND,NPART,NPTS,INERTIAL,&
         &PERIODIC_PARTICLES,LAGRANGIAN_HISTORY,STK
    namelist /numerics/          INTEGRATION_METHOD,TRUNCATION,DEALIASING,&
         &CRANK_NICHOLSON, FFT_DECOMPOSITION
    namelist /initialconditions/ INITCOND,KINITCOND,SEEDRANDOM
    namelist /inputoutput/       NOUTPUTFILES,INPUT_FIELD,INPUT_FIELD_FILENAME,&
         &NFILESTART,hdf5frate,slicefrate,COMPRESSION_LEVEL


    !open namelist input file
    open(newunit=aliakmon_nml,file='aliakmon.nml',action='read',err=100)

    read(aliakmon_nml,nml=general)
    read(aliakmon_nml,nml=hydro)
    read(aliakmon_nml,nml=force)
    read(aliakmon_nml,nml=passivescalar)
    read(aliakmon_nml,nml=radiation)
    read(aliakmon_nml,nml=magnetohydro)
    read(aliakmon_nml,nml=particle)
    read(aliakmon_nml,nml=numerics)
    read(aliakmon_nml,nml=initialconditions)
    read(aliakmon_nml,nml=inputoutput)

    n1=n
    n2=n
    n3=n
    gn2=n
    gn3=n
    
    if(.not.VISCOUS) FORCED=.false.
    if(.not.DIFFUSIVE) FORCED_PASSIVE_SCALAR=.false.
    if(.not.VISCOUS.or..not.DIFFUSIVE) CRANK_NICHOLSON=.false.

    close(aliakmon_nml, status='keep')

    return

100 print *, 'read_input_file:: Cant open file aliakmon.nml.'
    stop
  end subroutine read_namelist_file

  subroutine print_progress(ntimestep,t,&
       &tstart,t1)
    use data, only: nespec,nu1,nu3,nb1,nb3,fu,nn
    use parameters, only: dt, emean,kmax, sclvarprev
    use numerics
    use mpivars
    implicit none
    integer(ik), intent(IN)                                   :: ntimestep
    real(rk), intent(IN)                                      :: t,tstart
    real(rk), intent(INOUT)                                   :: t1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Print progress of the simulation on stdout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk)                                                  :: maxdivv
    real(rk)                                                  :: maxdivb
    real(rk), dimension(:),allocatable                        :: sclvar, lambdascl, PEl, nrgs, dissips
    real(rk)                                                  :: rmsb
    real(rk)                                                  :: t2
    real(rk)                                                  :: percent
    real(rk)                                                  :: etl
    real(rk)                                                  :: et,eta,etaoc
    real(rk)                                                  :: tol
    real(rk)                                                  :: tntimestep
    integer(ik)                                               :: edays, ehours
    integer(rk)                                               :: emins, esecs
    integer(ik)                                               :: ldays,lhours
    integer(ik)                                               :: lmins, lsecs
    real(rk), save                                            :: reprev=0.0
    real(rk), save                                            :: kenprev=1.0_rk
    real(rk), save                                            :: menprev
    real(rk)                                                  :: dfs
    real(rk)                                                  :: ten
    integer                                                   :: aux_dat=20
    integer                                                   :: auxscl_dat=222
    integer(ik) :: l
    

    dfs=0.05


    if(.not.allocated(sclvar)) allocate(sclvar(1:nn(4)))
    if(.not.allocated(lambdascl)) allocate(lambdascl(1:nn(4)))
    if(.not.allocated(PEl)) allocate(PEl(1:nn(4)))
    if(.not.allocated(nrgs)) allocate(nrgs(1:nn(4)))
    if(.not.allocated(dissips)) allocate(dissips(1:nn(4)))

    if(ntimestep==1) call timing(t1)
    call timing(t2)

    maxdivv=incompressibility(nn,fu,nu1)    
    if(MHD) maxdivb=incompressibility(nn,fu,nb1)

    !Calculate kinetic energy
    call msvalue(nn,fu,nrgs)
    KE=0.5_rk*nrgs(nu1)
    !Root-mean-square velocity
    rmsu=sqrt(2.0_rk/3.0_rk*KE)

    !Calculate magnetic energy 
    if(MHD) then
       ME=0.5_rk*nrgs(nb1)
       rmsb=sqrt(2.0_rk/3.0_rk*ME)
    else
       ME=0.0_rk
    end if

    sclvar(:)=0.0_rk
    !calculate passive scalar variance
    if(PASSIVE_SCALAR) then
       sclvar(nsclf:nscll)=sqrt(nrgs(nsclf:nscll))
    end if

    !Total energy
    ten=KE+ME

    !Calculate dissipation rates
    call mean_dissipation(nn,fu,dissips)
    emean=dissips(nu1)
    mkin_hel=mean_kinetic_helicity(nn,fu)
    mkhdis=mean_kinetic_helicity_dissipation(nn,fu)
    if(MHD)  then
       emeanb=dissips(nb1)
       mcross_hel=mean_cross_helicity(nn,fu)
       mchdis=mean_cross_helicity_dissipation(nn,fu)
       mmh=mean_magnetic_helicity(nn,fu)
       mmhdis=mean_magnetic_helicity_dissipation(nn,fu)
    end if
    if(PASSIVE_SCALAR) then
       emeanscl(nsclf:nscll)=dissips(nsclf:nscll)
    else
       emeanscl=0.0_rk
    end if

    !Caclulate integral length scale
    call integral_length_scale(nn,fu,ils,lambda,&
         &nespec)

    !Eddy turnover time
    if(rmsu==0.0_rk) rmsu=1.0_rk
    ETT=ils/(3*rmsu)
    !Taylor microscale
    if(emean==0.0_rk) emean=1.0_rk
    if(.not.MHD) lambda=sqrt(15._rk*visc(nu1)*rmsu**2/emean)
    if(MHD) then
       lambda=2.0_rk*PI*sqrt((KE+ME)/(emean/visc(nu1)+emeanb/visc(nb1)))
    end if
    !Taylor microscale Reynolds number
    REl=rmsu*lambda/visc(nu1)

    !Tolerance for Kaneda et al. (2004) forcing
    tol=1.0e-2_rk

    !Kaneda et al. (2004) forcing: keep kinetic energy within given bounds
    !by negative viscosity

    if(FORCED.and.VISCOUS.and.VARIABLE_FORCING) then
       if(KE<KENTAR-tol.and.KE<=kenprev) then
          fscale(nu1:nu3)=fscale(nu1:nu3)+dfs
       else if(KE>KENTAR+tol.and.KE>=kenprev) then
          fscale(nu1:nu3)=fscale(nu1:nu3)-dfs
       end if
    end if


    reprev=REl
    kenprev=KE

    !Kaneda et al. (2004) forcing: keep magnetic energy within given bounds
    !by negative resistivity
    if(MHD.and.FORCED_MHD.and.VARIABLE_FORCING) then
       if(ME<MENTAR-tol.and.ME<=menprev) then
          fscale(nb1:nb3)=fscale(nb1:nb3)+dfs
       else if(ME>MENTAR+tol.and.ME>=menprev) then
          fscale(nb1:nb3)=fscale(nb1:nb3)-dfs
       end if
    end if

    menprev=ME


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FIX THIS: to multiple scalars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(PASSIVE_SCALAR) then
       do l=nsclf,nscll
          if(emeanscl(l)==0.0_rk) emeanscl(l)=1.0_rk       
          !Taylor microscale of the passive scalar field
          lambdascl(l)=sqrt(3._rk*visc(l)*sclvar(l)**2/emeanscl(l))
          !Taylor microscale Peclet number
          PEl(l)=rmsu*lambdascl(l)/visc(l)
          !Kaneda et al. (2004) forcing: keep passive scalar mean-square  within
          !given bounds by negative diffusivity
          if(FORCED.and.VARIABLE_FORCING) then
             if(sclvar(l)<1.0-tol.and.sclvar(l)<=sclvarprev(l)) then
                fscale(l)=fscale(l)+dfs
             else if(sclvar(l)>1.0+tol.and.sclvar(l)>=sclvarprev(l)) then
                fscale(l)=fscale(l)-dfs
             end if
             sclvarprev(:)=sclvar(:)
          end if
       end do
    end if



    !Broadcast
#ifdef _MPI_
    sbuf(1:nn(4))=fscale(1:nn(4))
    call mpi_bcast(sbuf,int(nn(4)),MPIRK,MPIROOT,MPI_COMM_WORLD,mpierr)
    fscale(1:nn(4))=sbuf(1:nn(4))
#endif

    percent=0.0_rk
    !Calculate estimated time left for the simulation
    if(TIMESTEPS==0) then
       if(t/=tstart) then
          percent=((t-tstart)/TMAX)*100
          etl=(((t2-t1))/(t-tstart))*(TMAX-(t-tstart))
       end if
    else
       if(t2/=t1.and.ntimestep/=1) then
          percent=real(ntimestep)/TIMESTEPS*100.
          etl=(((t2-t1))/(ntimestep-1))*(TIMESTEPS-ntimestep+1)
       end if
    end if

    if(t2-t1>small) then
       et=(t2-t1)
    else
       et=0.0
    end if
    
     

    edays=int(et/(60._rk*60._rk*24._rk),ik)
    et=et-edays*(60._rk*60._rk*24._rk)
    ehours=int(et/(60._rk*60._rk),ik)
    et=et-ehours*(60._rk*60._rk)
    emins=int(et/60._rk,ik)
    esecs=et-emins*60._rk

    if(mpirank==MPIROOT) then
       ldays=int(floor(etl/(60._rk*60._rk*24._rk)),ik)
       etl=etl-ldays*(60._rk*60._rk*24._rk)
       lhours=int(floor(etl/(60._rk*60._rk)),ik)
       etl=etl-lhours*(60._rk*60._rk)
       lmins=int(floor(etl/60._rk),ik)
       lsecs=etl-lmins*60._rk
    end if

    !Kolmogorov microscale
    eta=(visc(nu1)**3/emean)**(1.0_rk/4.0_rk)


    if(ntimestep==1) then 
       tntimestep=2
    else
       tntimestep=ntimestep
    end if

    !Print timestep information
    if(mpirank==MPIROOT) then
       write(*,'(a)') ' '
       write(*,'(a)') '************************************************&
            &********************************'
       write(*,'(a)') '------------------------------------------------&
            &--------------------------------'
       write(*,&
            &'(a,f7.3,a,i3,a,i3,a,i3,a,i3,a,f8.3,a,i3,a,i3,a,i3,a,i3,a)') &
            & '| ', percent, ' %| t:',edays,' :',ehours,' :',emins,' :',&
            &esecs,' |t/ts:',((t2-t1))/(tntimestep-1), ' |ETL:', ldays,&
            &' :',lhours,' :',lmins,' :',lsecs, '       |' 
       write(*,'(a)') '------------------------------------------------&
            &--------------------------------'
       write(*,&
            &'(a,i5,a,f8.3,3(a,e9.2),a,f7.3,a)')  '| ',ntimestep, ' |', t,&
            &'|mc:',maxdivv,' |maxw:', MAXVORT,' |maxu: ', MAXVEL, '  |TE:', &
            &(KE+ME)*100._rk,'   |'
       write(*,'(a)') '------------------------------------------------&
            &--------------------------------'
       write(*,&
            & '((a,f7.3),(a,f7.3),2(a,e10.3),(a,f7.3),a)')  '| kmax*eta:',&
            & kmax*eta,' |Rel:',REl, ' | FHD:', fscale(nu1),' |dt: ', dt,&
            & ' |KE:',KE, '   |'
       write(*,'(a)') '------------------------------------------------&
            &--------------------------------'
       if(PASSIVE_SCALAR) then
          do l=nsclf,nscll
             write(*,'(a,i3)') 'Passive scalar #', l-nsclf+1
             etaoc=(visc(l)**3/emean)**(1.0_rk/4.0_rk)
             write(*,&
                  &'((a,f7.3),(a,f7.3),(a,f7.3),(a,e9.2),a)') &
                  & '| kmax*etaoc:',kmax*etaoc , ' |Pel:', PEl(l),&
                  &' |kmax*ls:', kmax*lambdascl(l), ' |sclvar:', sclvar(l),&
                  &'           |'
             write(*,'(a)') '------------------------------------------------&
                  &--------------------------------'
          end do
       end if

       if(MHD) then
          write(*,&
               & '(2(a,e8.1),(a,f8.2),a,e8.1,a)') '| maxdivb: ', maxdivb,&
               &'| maxj: ', MAXJ,'| emeanb: ', emeanb, '| mmh: ', mmh, &
               &'           |'
          write(*,'(a)') '------------------------------------------------&
               &--------------------------------'
          write(*,&
               & '(4(a,e8.1),a)') '| mch: ',mcross_hel/sqrt(4.0_rk*ME*KE),&
               &'| menergy: ', ME,'| maxlorentz: ', MAXLF,'| FMHD:',fscale(nb1),&
               &'        |'
          write(*,'(a)') '------------------------------------------------&
               &--------------------------------'

       end if

       write(*,'(5(a,f8.3),a)') '| RE= ', rmsu*ils/visc(nu1), '| eta= ', eta,&
            &'| lambda= ', lambda, '| L= ', ils,'| ETT=', ETT, '     |'
       write(*,'(a)') '------------------------------------------------&
            &--------------------------------'
       write(*,'(a)') '************************************************&
            &********************************'

       !Write data files
       write(hydro_dat,'(11(e36.8,a1))') t,'|',KE,'|',mkin_hel,'|',emean,'|',ils,'|',lambda,'|',eta,'|',&
            &rmsu*ils/visc(nu1),'|',REl,'|',mkhdis,'|',fscale(nu1)

       !FIX THIS: to multiple scalars
       if(PASSIVE_SCALAR) then
!!$          write(passive_dat,'(11e37.8)') t,meanscl,sclvar,maxscl,minscl,&
!!$               &emeanscl,maxgrad,PEl,&
!!$               &lambdascl,etaoc,fscale(nsclf:nscll)          
       end if

       if(MHD) then
          write(magnetic_dat,'(7e37.8)') t,ME,emeanb,mmh,mmhdis,mcross_hel,&
               &mchdis
       end if

       !Print inertial particle timescale
       if(PARTICLES.and.INERTIAL) then
          TP=STK*sqrt(visc(nu1)/emean)
          print *, 'Tp=',TP
       end if


    end if


    !Write data files
    if(mpirank==MPIROOT) then
       open(newunit=aux_dat, file='aux.dat', form='formatted',action='write')
       write(aux_dat,'(a,e20.10)')  'Viscosity:                   ',visc(nu1)
       write(aux_dat,'(a,e20.10)')  'Kolmogorov microscale:       ', eta
       write(aux_dat,'(a,e20.10)')  'Mean energy dissipation rate:', emean
       close(aux_dat,status='keep')
       if(PASSIVE_SCALAR) then
          open(newunit=auxscl_dat, file='auxscl.dat', form='formatted',&
               &action='write')
          write(auxscl_dat,'(a,3e20.10)')  'Prandtl number                     &
               &          :',PR
          write(auxscl_dat,'(a,3e20.10)')  'Obukhov-Corsin micorscale          &
               &          :',etaoc
          write(auxscl_dat,'(a,3e20.10)')  'Mean passive scalar variance dissip&
               &ation rate:',emeanscl
          close(auxscl_dat,status='keep')
       end if
    end if

    return
  end subroutine print_progress

  subroutine timing(time)
    use types
    use mpivars
#ifdef _OPENMP_
    use omp_lib
#endif
    implicit none
    real(rk), intent(OUT) :: time
!!!!!!!
    !Timing
!!!!!!!
    !MPI timing
#ifdef _MPI_

    time=mpi_wtime()

#else
    !OpenMP timing
#ifdef _OPENMP_
    time=omp_get_wtime()
#else
    !Serial timing
    call cpu_time(time)
#endif
#endif

    return

  end subroutine timing



end module input_output
