!$
!!$ (  / _ \\ \       /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!!$   Copyright (c) 2009-2020 George Momferatos

module fft_heffte
  use mpi
  use mpivars
#ifdef _CUFFT_
  use heffte_cufft
#elif defined _MKL_
  use heffte_mkl
#else
  use heffte_fftw
#endif
  use types
  use parameters
  use iso_c_binding
  implicit none
  !
  ! Heffte wrapper
  !
  !#ifdef _DOUBLE_
  real(c_double), dimension(:), allocatable, target :: input
  complex(c_double_complex), dimension(:), allocatable, target :: output
  complex(c_double_complex), dimension(:), allocatable :: work
  !$acc declare deviceptr(input, output, work)
  real(c_double), dimension(:,:,:), pointer :: r_array
  complex(c_double_complex), dimension(:,:,:), pointer :: c_array
!!$#else
!!$  real(c_float), dimension(:), allocatable, target :: input
!!$  complex(c_float_complex), dimension(:), allocatable, target :: output
!!$  complex(c_float_complex), dimension(:), allocatable :: work
!!$  !$acc declare deviceptr(input, output, work) 
!!$  real(c_float), dimension(:,:,:), pointer :: r_array
!!$  complex(c_float_complex), dimension(:,:,:), pointer :: c_array
!!$#endif
  integer(c_int)       :: slice_direction, pencil_direction,&
       & r2c_direction
  integer(c_int) :: il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3
  integer(c_long) :: size_in, size_out, size_work
#ifdef _CUFFT_
  type(heffte_fft3d_r2c_cufft) :: fft
  !type(heffte_fft3d_r2c_cufft), dimension(:), allocatable :: fftloc
#elif defined _MKL_
  type(heffte_fft3d_r2c_mkl) :: fft
#else
  type(heffte_fft3d_r2c_fftw) :: fft
#endif
  logical :: first
  !Interfaces of external subroutines in C

  interface
     subroutine heffte_init_slabs(n1,n2,n3,slice_direction&
          &,r2c_direction, il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1&
          &,oh2,oh3) bind(c)
       use iso_c_binding
       implicit none
       integer(c_long),value :: n1,n2,n3
       integer(c_int),value  :: slice_direction, r2c_direction
       integer(c_int)        :: il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3&
            &,oh1,oh2,oh3
     end subroutine heffte_init_slabs
  end interface

  interface
     subroutine heffte_init_pencils(n1,n2,n3,pencil_direction&
          &,r2c_direction, il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1&
          &,oh2,oh3) bind(c)
       use iso_c_binding
       implicit none
       integer(c_long),value :: n1,n2,n3
       integer(c_int),value  :: pencil_direction, r2c_direction
       integer(c_int)        :: il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3&
            &,oh1,oh2,oh3
     end subroutine heffte_init_pencils
  end interface


  interface
     subroutine heffte_set_num_device(numd) bind(c)
       use iso_c_binding
       implicit none
       integer(c_int), value :: numd
     end subroutine heffte_set_num_device
  end interface

contains

  subroutine fft_heffte_alloc(n1,n2,n3, ljsize, ljstart, lksize&
       &,lkstart)
    use iso_c_binding
    use mpi
    implicit none
    integer(ik), intent(in)  :: n1,n2,n3
    integer(ik), intent(out) :: ljsize, ljstart
    integer(ik), intent(out) :: lksize, lkstart
!!!!!!!!!!!!!!!!!!!!!!!!
    !Allocate FFT structures
!!!!!!!!!!!!!!!!!!!!!!!!
    first = .true.
    r2c_direction=0
    slice_direction=2
    pencil_direction = 0

    select case(FFT_DECOMPOSITION)
    case(SLABS)
       call heffte_init_slabs(int(n1,c_long),int(n2,c_long),int(n3&
            &,c_long), slice_direction,r2c_direction, il1,il2,il3,ih1&
            &,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3)
    case(PENCILS)
       call heffte_init_pencils(int(n1,c_long),int(n2,c_long),int(n3&
            &,c_long), pencil_direction,r2c_direction, il1,il2,il3&
            &,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3)
    end select

    ljstart = il2
    ljsize = ih2 - il2 + 1
    lkstart = il3
    lksize = ih3 - il3 + 1

    return

!!$#ifdef _CUFFT_
!!$    fft = heffte_fft3d_r2c_cufft(il1, il2, il3, ih1, ih2, $
    !!&ih3, 0, 1, 2, ol1, ol2, ol3, oh1, oh2, oh3, 0, 1, 2, $
    !!&r2c_direction, MPI_COMM_WORLD, $         &.false., .true.,
    !!.false.)
!!$    fft = heffte_fft3d_r2c_cufft(il1,il2,il3,ih1,ih2,ih3, $
    !!&ol1,ol2,ol3,oh1,oh2,oh3 $
    !!&,r2c_direction,MPI_COMM_WORLD,use_alltoall=.false.)
!!$#elif defined _MKL_
!!$    fft = heffte_fft3d_r2c_mkl(il1, il2, il3, ih1, ih2, $
    !!&ih3, 0, 1, 2, ol1, ol2, ol3, oh1, oh2, oh3, 0, 1, 2, $
    !!&r2c_direction, MPI_COMM_WORLD, $         &.false., .true.,
    !!.false.)
!!$    fft = heffte_fft3d_r2c_mkl(il1,il2,il3,ih1,ih2,ih3, $
    !!&ol1,ol2,ol3,oh1,oh2,oh3 $
    !!&,r2c_direction,MPI_COMM_WORLD)
!!$#else
!!$    fft = heffte_fft3d_r2c_fftw(il1,il2,il3,ih1,ih2,ih3, $
    !!&ol1,ol2,ol3,oh1,oh2,oh3 $
    !!&,r2c_direction,MPI_COMM_WORLD)
!!$#endif



!!$acc enter data create(input(1:size_in), output(1:size_out),
    !!work(1:size_work))

    !call c_f_pointer(c_loc(input),r_array,shape=(/n1, ljsize,
    !lksize/))
    !call c_f_pointer(c_loc(output),c_array,shape=(/dim1(n1)/2,
    !ljsize, lksize/))

    return

  end subroutine fft_heffte_alloc

  subroutine fft_heffte_dealloc
    use iso_c_binding
    !
    !Deallocate FFT structures
    !
    !deallocate(input)
    !deallocate(output)
    !deallocate(work)
!!$acc exit data delete(input, output, work)

    !    call fft%release()

    return

  end subroutine fft_heffte_dealloc


  subroutine fft_heffte_fourier(nn, dir, fu, nfs, nfe)
    use parameters, only : dim1
    use iso_c_binding
    use omp_lib
    use cudafor
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    integer(ik), intent(in)       :: dir
    integer(ik), optional :: nfs, nfe
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)),&
         & intent(inout) :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: nnfs,nnfe,nfi
    integer(ik) :: i,j,k,iii
    complex(ck) :: ctmp
    integer :: nthreads, itmp, nthread
    real(rk) :: scale


    r2c_direction=0
    slice_direction=2
    pencil_direction = 0




    nnfs=1
    nnfe=nn(4)
    if(present(nfs)) nnfs=nfs
    if(present(nfe)) nnfe=nfe
    nthreads=nnfe-nnfs+1



    if(first) then

#ifdef _CUFFT_
       fft = heffte_fft3d_r2c_cufft(il1, il2, il3, ih1, ih2, ih3, 0,&
            & 1, 2, ol1, ol2, ol3, oh1, oh2, oh3, 0, 1, 2,&
            & r2c_direction, MPI_COMM_WORLD, .false., .true., .true.)
#elif defined _MKL_
       fft = heffte_fft3d_r2c_mkl(il1, il2, il3, ih1, ih2, ih3, 0, 1,&
            & 2, ol1, ol2, ol3, oh1, oh2, oh3, 0, 1, 2, r2c_direction&
            &, MPI_COMM_WORLD, .false., .true., .true.)
#endif
       size_in=fft%size_inbox()
       size_out=fft%size_outbox()
       size_work=fft%size_workspace()


       if(.not.allocated(input)) allocate(input(1:size_in))
       if(.not.allocated(output)) allocate(output(1:size_out))
       if(.not.allocated(work)) allocate(work(1:size_work))

       !$acc enter data create(input(1:size_in), output(1:size_out), work(1:size_work))

       !$acc data present(input(1:size_in), output(1:size_out), work(1:size_work))
       !$acc update device(input(1:size_in)) 
       !$acc host_data use_device(input(1:size_in), output(1:size_out), work(1:size_work))
       call fft%forward(input,output,work, scale_cufft_none)
       !$acc end host_data
       !$acc update self(output(1:size_out))

       !$acc update device(output(1:size_out)) 
       !$acc host_data use_device(input(1:size_in), output(1:size_out), work(1:size_work))
       call fft%backward(output, input, work, scale_cufft_full)
       !$acc end host_data
       !$acc update self(input(1:size_in))
       !$acc end data

       !$acc exit data delete(input, output, work)

       if(allocated(input)) deallocate(input)
       if(allocated(output)) deallocate(output)
       if(allocated(work)) deallocate(work)



       first = .false.



    end if

    !$omp parallel num_threads(nthreads) &
    !$omp& default(none) &
    !$omp& shared(fu, numdevices, n1, nn, dir, nnfs, size_in, &
    !$omp& size_out, size_work, ljsize, nnfe, lksize, &
    !$omp& scale_cufft_none, scale_cufft_full, mpirank) &
    !$omp& private(nthread, numdevice, input, output, work, c_array, iii, &
    !$omp& ctmp, cudaerror, r_array, nfi, i, j, k) &
    !$omp& shared(fft)

    nthread = omp_get_thread_num()

    nfi = nthread + 1

#ifdef _CUFFT_
    cudaerror = cudagetdevicecount(numdevices)
    if(numdevices > 1) then
       numdevice = int(mod(nthread, numdevices), i4b)
    else
       numdevice = 0
    end if
    cudaerror = cudasetdevice(numdevice)
    cudaerror = cudagetdevice(numdevice)
    call heffte_set_num_device(numdevice)
!!$    if(mpirank == mpiroot) print '(3(a,i3),a)', 'Thread # ', nthread, &
!!$         & 'uses GPU # ', numdevice, ' out of ', numdevices, ' visible.'
#endif

    allocate(input(1:size_in))
    allocate(output(1:size_out))
    allocate(work(1:size_work))

    call c_f_pointer(c_loc(input),r_array,shape=(/n1, ljsize, lksize/))
    call c_f_pointer(c_loc(output),c_array,shape=(/dim1(n1)/2, ljsize, lksize/))
    !$acc enter data create(input(1:size_in), output(1:size_out), work(1:size_work))

    !$acc data present(input(1:size_in), output(1:size_out), work(1:size_work))
    !$omp do ordered schedule(static, 1)
    do nfi=nnfs,nnfe

       if(dir == 1) then



          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             r_array(i,j,k)=fu(i,j,k,nfi)
          end do; end do ; end do 




          !$omp ordered
          !$acc update device(input(1:size_in)) 
          !$acc host_data use_device(input(1:size_in), output(1:size_out), work(1:size_work))
          call fft%forward(input,output,work, scale_cufft_none)
          !$acc end host_data
          !$acc update self(output(1:size_out))
          !$omp end ordered 



          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
             iii = (i-1) / 2 + 1
             ctmp = c_array(iii,j,k)
             fu(i,j,k,nfi) = real(ctmp)
             fu(i+1,j,k,nfi) = aimag(ctmp)
          end do; end do ; end do


       else 

          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
             iii = (i - 1) / 2  + 1
             c_array(iii,j,k) = cmplx(fu(i,j,k,nfi), fu(i+1,j,k,nfi))
          end do; end do ; end do



          !$omp ordered
          !$acc update device(output(1:size_out)) 
          !$acc host_data use_device(input(1:size_in), output(1:size_out), work(1:size_work))
          call fft%backward(output, input, work, scale_cufft_full)
          !$acc end host_data
          !$acc update self(input(1:size_in))
          !$omp end ordered



          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             fu(i,j,k,nfi) = r_array(i,j,k)
          end do; end do ; end do
       end if
    end do
    !$omp end do
    !$acc end data

    !$acc exit data delete(input, output, work)
    deallocate(input)
    deallocate(output)
    deallocate(work)


    !$omp end parallel


    return


  end subroutine fft_heffte_fourier




end module fft_heffte
