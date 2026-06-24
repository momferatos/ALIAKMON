!!$     ___ __                                       
!!$ (  / _ \\ \        /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!$Copyright (c) 2009-2020 Georgios Momferatos
module validation
  use types
  use parameters
  use numerics, only: mean_kinetic_helicity, mean_kinetic_helicity_dissipation,&
       &mean_cross_helicity, mean_cross_helicity_dissipation,&
       &mean_magnetic_helicity, mean_magnetic_helicity_dissipation,&
       &ambipolar_diffusion_dissipation,curl,fourier,vector_potential,msvalue,&
       &mean_dissipation,cross_product
  implicit none
  !
  ! validation routines 
  !
  
contains

  subroutine dissipation_test
    use types
    use data, only: fu,nn
    use mpivars
    implicit none
    !
    !check dissipation rates, useful for validation
    !
    real(rk)       :: de,dedt,ken,men,ten,mkh,mch
    real(rk)       :: err,tdis
    real(rk)       :: dmkh,dmkhdt,mkherr,dmmh,dmmhdt,mmherr
    real(rk)       :: dmch,dmchdt,mcherr
    real(rk), dimension(1:nn(4)) :: nrgs, dissips
    real(rk), save :: e1,e2,e3,mkh1,mkh2,mkh3,mmh1,mmh2,mmh3,mch1,mch2
    real(rk), save :: mch3,vdis1,odis1,addis1,einj
    integer(ik), save :: icheck=1
    
    !calculate energies
    call msvalue(nn,fu,nrgs)
    ken=0.5*nrgs(nu1)
    if(MHD) then
       men=0.5*nrgs(nb1)
       !calculate mean magnetic helicity
       mmh=mean_magnetic_helicity(nn,fu)
       !calculate mean cross-helicity
       mch=mean_cross_helicity(nn,fu)
    else
       men=0.0_rk
    end if
    
    mkh=mean_kinetic_helicity(nn,fu)
    
    !calculate total energy
    ten=ken+men



    !if the subroutine is called for the first time
    if(icheck==1) then
       e1=ten
       mkh1=mkh
       mmh1=mmh
       mch1=mch
       icheck=2
       return
       !if the subroutine is called for the second time
    else if(icheck==2) then
       e2=ten
       mkh2=mkh
       mmh2=mmh
       mch2=mch
       icheck=3
       !if the simulation is viscous
       if(VISCOUS) then
          !calculate mean dissipation rates
          call mean_dissipation(nn,fu,dissips)
          vdis1=dissips(1)
          odis1=dissips(1)/MAGNETIC_PR
          !calculate mean kinetic helicity dissipation
          mkhdis=mean_kinetic_helicity_dissipation(nn,fu)
          !if the simulation is forced
          if(FORCED) then
             !calculate energy injection by the forcing term
             !einj=energy_injection(nn,fu)
          else
             einj=0.0_rk
          end if
       else
          vdis1=0.0_rk
       end if
       !Magnetohydrodynamics
       if(MHD) then
          !If the simulation is resistive
          if(RESISTIVE) then             
             !Calculate mean magnetic helicity dissipation 
             mmhdis=mean_magnetic_helicity_dissipation(nn,fu)
             !Calculate mean cross-helicity dissipation 
             mchdis=mean_cross_helicity_dissipation(nn,fu)
          end if

          !Ambipolar diffusion
          if(AMB_DIFF) then
             !Calculate ambipolar diffusion dissipation
             addis1=ambipolar_diffusion_dissipation(nn,fu)            
          else
             addis1=0.0_rk
          end if
       else
          odis1=0.0_rk
       end if
       return
    else

       !Total energy
       e3=ten
       de=(e3-e1)
       dedt=de/(2.*dt)
       e1=e2
       e2=e3

       !Kinetic helicity
       mkh3=mkh
       dmkh=(mkh3-mkh1)
       dmkhdt=dmkh/(2.*dt)
       mkh1=mkh2
       mkh2=mkh3

       !Magnetic helicity
       mmh3=mmh
       dmmh=(mmh3-mmh1)
       dmmhdt=dmmh/(2.*dt)
       mmh1=mmh2
       mmh2=mmh3

       !Cross-helicity
       mch3=mch
       dmch=(mch3-mch1)
       dmchdt=dmch/(2.*dt)
       mch1=mch2
       mch2=mch3

    end if

    !print test results
    if(mpirank==MPIROOT) then
       tdis=(vdis1+odis1+addis1)
       err=dedt+tdis-einj
       mkherr=dmkhdt+mkhdis
       mmherr=dmmhdt+mmhdis
       mcherr=dmchdt+mchdis

       !If there is dissipation
       if(VISCOUS.or.RESISTIVE) then
          print '(a,e10.3,f9.3,a)',  "E test: ", err, err/dedt*100,'%'
          if(.not.MHD) then
             print '(a,e10.3,f9.3,a,e10.3)',  "MKH test: ",mkherr, &
                  &mkherr/dmkhdt*100,'%'
          end if
          if(MHD) then
             print '(a,e10.3,f9.3,a,e10.3)',  "MMH test: ",&
                  &mmherr, mmherr/dmmhdt*100,&
                  &'%'  
             print '(a,e10.3,f9.3,a,e10.3)', "MCH test: ",&
                  &mcherr, mcherr/dmchdt*100,'%'   
          end if

          write(distest_dat,'(6e17.8)') err,err/dedt*100,mmherr,&
               &err/dmmhdt*100,mcherr,mcherr/dmchdt*100
          !If there is no dissipation
       else
          err=ten-teinitial
          print '(a,2e12.3,a)', 'TE',err,err/teinitial*100,'%' 
          if(.not.MHD) then
             err=mkh-mkhinitial
             print '(a,2e12.3,a)', 'MKH', err,err/mkhinitial*100,'%'
          else
             err=mch-mchinitial
             print '(a,2e12.3,a)', 'MCH', err,err/mchinitial*100,'%'
             err=mmh-mmhinitial
             print '(a,2e12.3,a)', 'MMH', err,err/mmhinitial*100,'%'
          end if
       end if

    end if


    if(VISCOUS) then
       call mean_dissipation(nn,fu,dissips)
       vdis1=dissips(1)
       odis1=dissips(1)/MAGNETIC_PR
       mkhdis=mean_kinetic_helicity_dissipation(nn,fu)
       !einj=energy_injection(nn,fu)
    else
       vdis1=0.0_rk
       odis1=0.0_rk
       mkhdis=0.0_rk
    end if
    if(MHD) then
       if(RESISTIVE) then
          mmhdis=mean_magnetic_helicity_dissipation(nn,fu)
          mchdis=mean_cross_helicity_dissipation(nn,fu)
       end if
       if(AMB_DIFF) then
          addis1=ambipolar_diffusion_dissipation(nn,fu)            
       else
          addis1=0.0_rk
       end if
    else
       odis1=0.0_rk
    end if

    return

  end subroutine dissipation_test

  subroutine energy_test(nn,u,fu,rhs)
    use types
    use mpivars
    use parameters,only:MHD
    use data, only:arr_en_1,nu1,nu3,nb1,nb3,arr_en_2, arr_en_3
    implicit none
    integer(ik), dimension(1:4),intent(in) :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT) :: u
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN) :: fu
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN) :: rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !check energy conservation in the ideal limit, useful for validation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk)                                               :: enerr,hderr,mhderr
    real(rk)                                               :: mcherr,mcherr1,mcherr2
    real(rk)                                               :: mmherr,mmherr1,mmherr2
    real(rk)                                               :: mkherr,mkherr1, mkherr2
    integer(ik)                                            :: i,j,k,l
    integer(ik)                                            :: n1,n2,n3

    n1=nn(1)
    n2=nn(2)
    n3=nn(3)

    !Copy arrays from Fourier space to physical space
    !$omp parallel do
    do l=1,nn(4) ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
       u(i,j,k,l)=fu(i,j,k,l)
       arr_en_1(i,j,k,l)=0.0_rk
       arr_en_2(i,j,k,l)=0.0_rk
       arr_en_3(i,j,k,l)=rhs(i,j,k,l)
    end do; end do ; end do ; end do 
    !$omp end parallel do
    !Vorticity
    call curl(nn,arr_en_2,fu,nu1)
    !curl of Navier-Stokes rhs
    call curl(nn,arr_en_1,arr_en_3,nu1)

    !Perform Fourier transforms
    call fourier(nn,-1_ik,u)
    call fourier(nn,-1_ik,arr_en_2)
    call fourier(nn,-1_ik,arr_en_3)
    call fourier(nn,-1_ik,arr_en_1)

    !Calculate conservation errors 
    hderr=0.0
    mkherr1=0.0
    mkherr2=0.0
    !$omp parallel do reduction(+:hderr,mkherr1,mkherr2)
    do k=1,n3
       do j=1,n2
          do i=1,n1
             ! velocity . rhs of momentum equation
             hderr=hderr+dot_product(u(i,j,k,nu1:nu3),arr_en_3(i,j,k,nu1:nu3))
             !Mean kinetic helicity
             !vorticity . rhs of momentum equation
             mkherr1=mkherr1+dot_product(arr_en_2(i,j,k,nu1:nu3),arr_en_3(i,j,k,nu1:nu3))
             !velocity * rhs of momentum equation
             mkherr2=mkherr2+dot_product(u(i,j,k,nu1:nu3),arr_en_1(i,j,k,nu1:nu3))
          end do
       end do
    end do
    !$omp end parallel do
    hderr=hderr/(n1*n2*gn3)
    mkherr1=mkherr1/(n1*n2*gn3)
    mkherr2=mkherr2/(n1*n2*gn3)

    mhderr=0.0_rk
    !Magnetohydrodynamics
    if(MHD) then
       !Calculate conservation errors
       mhderr=0.0
       mcherr1=0.0
       mcherr2=0.0
       !$omp parallel do reduction(+:mhderr,mcherr1,mcherr2)
       do k=1,n3
          do j=1,n2
             do i=1,n1
                !Magnetic helicity
                !magnetic field . rhs of induction equation
                mhderr=mhderr+dot_product(u(i,j,k,nb1:nb3),arr_en_3(i,j,k,nb1:nb3))
                !Mean cross-helicity
                !magnetic field . rhs of momentum equation
                mcherr1=mcherr1+dot_product(u(i,j,k,nb1:nb3),arr_en_3(i,j,k,nu1:nu3))
                !velocity . magnetc field
                mcherr2=mcherr2+dot_product(u(i,j,k,nu1:nu3),u(i,j,k,nb1:nb3))
             end do
          end do
       end do
       !$omp end parallel do
       mhderr=mhderr/(n1*n2*gn3)
       mcherr1=mcherr1/(n1*n2*gn3)
       mcherr2=mcherr2/(n1*n2*gn3)

       !Mean magnetic helicity
       ! B = nabla x A (arr_en_1 = A)
       call vector_potential(nn,arr_en_1,fu,nb1)
       call fourier(nn,-1_ik,arr_en_1,nb1)
       ! arr_en_2 = A x B
       call cross_product(nn,arr_en_2,arr_en_1,fu,nb1)
       mmherr1=0.0_rk
       mmherr2=0.0_rk
       !$omp parallel do reduction(+:mmherr1,mmherr2)
       do k=1,n3
          do j=1,n2
             do i=1,n1
                !vector potential * rhs of induction equationm
                mmherr1=mmherr1+dot_product(arr_en_1(i,j,k,nb1:nb3),arr_en_3(i,j,k,nb1:nb3))
                !magnetic field * (vector potential x magnetic field)
                mmherr2=mmherr2+dot_product(u(i,j,k,nb1:nb3),arr_en_2(i,j,k,nb1:nb3))
             end do
          end do
       end do
       !$omp end parallel do
       mmherr1=mmherr1/(n1*n2*gn3)
       mmherr2=mmherr2/(n1*n2*gn3)
    end if

    enerr=hderr+mhderr
    mkherr=mkherr1+mkherr2
    mcherr=mcherr1+mcherr2
    mmherr=mmherr1+mmherr2
#ifdef _MPI_
    call mpi_allreduce(MPI_IN_PLACE,enerr,1,MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,mkherr,1,MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,mcherr,1,MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,mmherr,1,MPIRK,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
    !Print errors
    if(mpirank==mpiroot) then
       print *, 'energy error: ', enerr
       if(.not.MHD) print *, 'mkh error: ', mkherr
       if(MHD) then
          print *, 'mch error   : ', mcherr
          print *, 'mmh error: ', mmherr
       end if
    end if

    arr_en_1=0.0
    arr_en_2=0.0
    arr_en_3=0.0
    return

  end subroutine energy_test

end module validation
