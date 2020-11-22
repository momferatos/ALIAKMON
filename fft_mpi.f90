!!$     ___ __                                       
!!$ (  / _ \\ \        /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!$Copyright (c) 2009-2020 Georgios Momferatos
module fft_mpi
#ifdef _MPI_
  use types
  use mpi
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3-mpi.f03'
#ifdef _DOUBLE_
  real(c_double), dimension(:,:,:), pointer :: array
  complex(c_double_complex), dimension(:,:,:), pointer   :: carray
#else
  real(c_float), dimension(:,:,:), pointer :: array
  complex(c_float_complex), dimension(:,:,:), pointer :: carray
#endif
  type(c_ptr) :: plan, iplan, p
contains
  subroutine fft_mpi_alloc(n1, n2, gn3, lksize, lkstart)
    use omp_lib
    use parameters, only: dim1
    implicit none
    integer(ik), intent(in) :: n1, n2, gn3
    integer(ik), intent(out) :: lksize, lkstart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(c_int) :: nt, error
    integer(c_intptr_t) :: nn1, nn2, gnn3
    integer(c_intptr_t) :: alloc_local, local_n1, local_n2, local_1_start, local_2_start
    integer(ik) :: n3

#ifdef _DOUBLE_
    call fftw_mpi_init()
#else
    call fftwf_mpi_init()
#endif

    
#ifdef _OPENMP_
    !$omp parallel    
    nt = int(omp_get_num_threads(), c_int)
    !$omp end parallel
#ifdef _DOUBLE_
    error = fftw_init_threads()
#else
    error = fftwf_init_threads()
#endif
#endif

    nn1 = int(n1, c_intptr_t)
    nn2 = int(n2, c_intptr_t)
    gnn3 = int(gn3, c_intptr_t)


#ifdef _DOUBLE_
    alloc_local = fftw_mpi_local_size_3d_transposed(nn1,nn2,gnn3,MPI_COMM_WORLD,local_n1,local_1_start,local_n2,local_2_start);
#else
    alloc_local = fftwf_mpi_local_size_3d_transposed(nn1,nn2,gnn3,MPI_COMM_WORLD,local_n1,local_1_start,local_n2,local_2_start);
#endif

    

    lksize=local_n1;
    n3=lksize
    lkstart=local_1_start

#ifdef _DOUBLE_    
    p = fftw_alloc_real(int(dim1(n1) * n2 * n3, C_SIZE_T))
#else
    p = fftwf_alloc_real(int(dim1(n1) * n2 * n3, C_SIZE_T))
#endif

    call c_f_pointer(p, array, shape=(/ dim1(n1), n2, lksize /))
    call c_f_pointer(p, carray, shape=(/ dim1(n1)/2, n2, lksize /))


#ifdef _OPENMP_
#ifdef _DOUBLE_
    call fftw_plan_with_nthreads(nt)
#else
    call fftwf_plan_with_nthreads(nt)
#endif
#endif
  
#ifdef _DOUBLE_
    plan = fftw_mpi_plan_dft_r2c_3d(nn1, nn2, gnn3, array, carray,MPI_COMM_WORLD,ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT))
    iplan = fftw_mpi_plan_dft_c2r_3d(nn1, nn2, gnn3, carray, array,MPI_COMM_WORLD,ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN))
#else
    plan = fftwf_mpi_plan_dft_r2c_3d(gnn3, nn2, nn1, array, carray,MPI_COMM_WORLD,ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT))
    iplan = fftwf_mpi_plan_dft_c2r_3d(gnn3, nn2, nn1, carray, array,MPI_COMM_WORLD,ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN))
#endif

    return

  end subroutine fft_mpi_alloc

  subroutine fft_mpi_dealloc
    implicit none


#ifdef _DOUBLE_
    call fftw_free(p)
    call fftw_destroy_plan(plan)
    call fftw_destroy_plan(iplan)
#ifdef _OPENMP_
    call fftw_cleanup_threads();
#endif
#else
    call fftwf_free(p)
    call fftwf_destroy_plan(plan)
    call fftwf_destroy_plan(iplan)
#ifdef _OPENMP_
    call fftwf_cleanup_threads();
#endif
#endif

    return

  end subroutine fft_mpi_dealloc

  subroutine fft_mpi_fourier(nn, gn3, dir, fu, nfs, nfe)
    use parameters, only : dim1
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    integer(ik), intent(in) :: gn3
    integer(ik), intent(in)       :: dir
    integer(ik), optional :: nfs, nfe
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(inout) :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk)                  :: scale
    integer(ik) :: nnfs,nnfe,nfi
    integer(ik) :: i,j,k
    
    if(dir == -1) then
       scale = 1.0_rk / real(nn(1)*nn(2)*gn3, rk)
    else
       scale = 1.0_rk
    end if

    nnfs=1
    nnfe=nn(4)
    if(present(nfs)) nnfs=nfs
    if(present(nfe)) nnfe=nfe

    
    do nfi=nnfs,nnfe
       if(dir == 1) then
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             array(i,j,k) = fu(i,j,k, nfi)
          end do ; end do ; end do 
          !$omp end parallel do
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=nn(1)+1,dim1(nn(1))
             array(i,j,k) = 0.0_rk
          end do ; end do ; end do
          !$omp end parallel do
       else 
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
             array(i,j,k) = fu(i,j,k, nfi)
          end do ; end do ; end do
          !$omp end parallel do
       end if
       
#ifdef _DOUBLE_
       if(dir == 1) then
          call fftw_execute_dft_r2c(plan, array, carray)
       else
          call fftw_execute_dft_c2r(iplan, carray, array)
       end if
#else
       if(dir == 1) then
          call fftwf_execute_dft_r2c(plan, array, carray)
       else
          call fftwf_execute_dft_c2r(iplan, carray, array)
       end if
#endif


       if(dir == -1) then
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             fu(i,j,k,nfi) = scale * array(i,j,k)
          end do ; end do ; end do
          !$omp end parallel do
       else
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
             fu(i,j,k,nfi) = array(i,j,k)
           end do ; end do ; end do
          !$omp end parallel do
       end if


    end do

    return

  end subroutine fft_mpi_fourier
  
#endif
end module fft_mpi

