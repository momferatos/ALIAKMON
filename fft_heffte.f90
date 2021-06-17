!!$     ___ __                                       
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
  real(c_float), dimension(:), allocatable, target :: input
  complex(c_float_complex), dimension(:), allocatable, target :: output
  complex(c_float_complex), dimension(:), allocatable :: work
  !$acc declare present(input, output, work) 
  real(c_float), dimension(:,:,:), pointer :: r_array
  complex(c_float_complex), dimension(:,:,:), pointer :: c_array
  integer(c_long) :: size_in, size_out, size_work
  
#ifdef _CUFFT_
  type(heffte_fft3d_r2c_cufft) :: fft
#else
  type(heffte_fft3d_r2c_fftw) :: fft
#endif
  !Interfaces of external subroutines in C

  interface
     subroutine heffte_init_slabs(n1,n2,n3,slice_direction,r2c_direction,&
          &il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3) bind(c)
       use iso_c_binding
       implicit none
       integer(c_long),value :: n1,n2,n3
       integer(c_int),value  :: slice_direction, r2c_direction
       integer(c_int)        :: il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3
     end subroutine heffte_init_slabs
  end interface

  interface
     subroutine heffte_init_pencils(n1,n2,n3,pencil_direction,r2c_direction,&
          &il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3) bind(c)
       use iso_c_binding
       implicit none
       integer(c_long),value :: n1,n2,n3
       integer(c_int),value  :: pencil_direction, r2c_direction
       integer(c_int)        :: il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3
     end subroutine heffte_init_pencils
  end interface

contains

  subroutine fft_heffte_alloc(n1,n2,n3, ljsize, ljstart, lksize,lkstart)
    use iso_c_binding
    use mpi
    implicit none
    integer(ik), intent(in)  :: n1,n2,n3
    integer(ik), intent(out) :: ljsize, ljstart
    integer(ik), intent(out) :: lksize, lkstart
!!!!!!!!!!!!!!!!!!!!!!!!
    !Allocate FFT structures
!!!!!!!!!!!!!!!!!!!!!!!!
    integer(c_long)          :: nn1,nn2,nn3
    integer(c_long), target  :: nlksize,nlkstart
    integer(c_long) :: ndir
    integer(c_long) :: alloc
    integer(c_int)       :: slice_direction, pencil_direction, r2c_direction
    integer(c_int) :: il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3

    r2c_direction=0
    slice_direction=2
    pencil_direction = 0

    select case(FFT_SUBDIVISION)
    case(SLABS)
       call heffte_init_slabs(&
            &int(n1,c_long),int(n2,c_long),int(n3,c_long),&
            &slice_direction,r2c_direction,&
            &il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3)
    case(PENCILS)
       call heffte_init_pencils(&
            &int(n1,c_long),int(n2,c_long),int(n3,c_long),&
            &pencil_direction,r2c_direction,&
            &il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3)
    end select

    ljstart = il2
    ljsize = ih2 - il2 + 1
    lkstart = il3
    lksize = ih3 - il3 + 1

#ifdef _CUFFT_
    fft = heffte_fft3d_r2c_cufft(il1,il2,il3,ih1,ih2,ih3,&
         &ol1,ol2,ol3,oh1,oh2,oh3&
         &,r2c_direction,MPI_COMM_WORLD)
#else
    fft = heffte_fft3d_r2c_fftw(il1,il2,il3,ih1,ih2,ih3,&
         &ol1,ol2,ol3,oh1,oh2,oh3&
         &,r2c_direction,MPI_COMM_WORLD)
#endif
    size_in=fft%size_inbox()
    size_out=fft%size_outbox()
    size_work=fft%size_workspace()

    allocate(input(size_in))
    allocate(output(size_out))
    allocate(work(size_work))

    call c_f_pointer(c_loc(input),r_array,shape=(/n1, ljsize, lksize/))
    call c_f_pointer(c_loc(output),c_array,shape=(/dim1(n1)/2, ljsize, lksize/))
    !$acc enter data create(input(1:size_in), output(1:size_out), work(1:size_work))
    return

  end subroutine fft_heffte_alloc

  subroutine fft_heffte_dealloc
    use iso_c_binding
    !
    !Deallocate FFT structures
    !
    deallocate(input)
    deallocate(output)
    deallocate(work)
    !$acc exit data delete(input, output, work)

    call fft%release()

    return

  end subroutine fft_heffte_dealloc


  subroutine fft_heffte_fourier(nn, dir, fu, nfs, nfe)
    use parameters, only : dim1
    use iso_c_binding
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    integer(ik), intent(in)       :: dir
    integer(ik), optional :: nfs, nfe
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(inout) :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: nnfs,nnfe,nfi
    integer(ik) :: i,j,k,iii
    integer(c_long)                     :: nn1,nn2,ngn3,nlksize,nlkstart,ndir,alloc
    type(c_ptr) :: gpu_input_c, gpu_output_c, gpu_workspace_c
    type(c_ptr) :: gpu_input, gpu_output, gpu_workspace
    type(c_funptr) :: func_ptr
    complex(ck) :: ctmp


    nnfs=1
    nnfe=nn(4)
    if(present(nfs)) nnfs=nfs
    if(present(nfe)) nnfe=nfe

    !$acc data present(input(1:size_in), output(1:size_out), work(1:size_work))
    if(dir == 1) then
       do nfi=nnfs,nnfe

          !$omp parallel do 
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             r_array(i,j,k)=fu(i,j,k,nfi)
          end do; end do ; end do 
          !$omp end parallel do

#ifdef _CUFFT_
          !$acc update device(input(1:size_in)) 
          !$acc host_data use_device(input(1:size_in), output(1:size_out), work(1:size_work))
          call fft%forward(input,output,work, scale_cufft_none)
          !$acc end host_data
          !$acc update host(output(1:size_out))
#else
          call fft%forward(input,output,work, scale_fftw_none)
#endif

          !omp parallel do private(iii, ctmp)
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
             iii = (i-1) / 2 + 1
             ctmp = c_array(iii,j,k)
             fu(i,j,k,nfi) = real(ctmp)
             fu(i+1,j,k,nfi) = aimag(ctmp)
          end do; end do ; end do
          !!omp end parallel do
       end do


    else 

       do nfi=nnfs,nnfe


          !$omp parallel do private(iii) 
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
             iii = (i - 1) / 2  + 1
             c_array(iii,j,k) = cmplx(fu(i,j,k,nfi), fu(i+1,j,k,nfi))
          end do; end do ; end do
          !omp end parallel do
#ifdef _CUFFT_
          !$acc update device(output(1:size_out)) 
          !$acc host_data use_device(input(1:size_in), output(1:size_out), work(1:size_work))
          call fft%backward(output, input, work, scale_cufft_full)
          !$acc end host_data
          !$acc update host(input(1:size_in))
#else
          call fft%backward(output,input,work, scale_fftw_full)
#endif



          !$omp parallel do 
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             fu(i,j,k,nfi) = r_array(i,j,k)
          end do; end do ; end do
          !$omp end parallel do


       end do

    end if
    !$acc end data


    return

  end subroutine fft_heffte_fourier




end module fft_heffte
