!!$     ___ __                                       
!!$ (  / _ \\ \       /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!!$   Copyright (c) 2009-2020 George Momferatos

module fft_cuda
  use mpi
  use mpivars
  use heffte_cufft
  use types
  use parameters
  use iso_c_binding
  implicit none
  !
  ! Heffte CUDA Fortran wrapper
  !

  
  real(c_float), dimension(:), allocatable, target :: input
  complex(c_float_complex), dimension(:), allocatable, target :: output
  complex(c_float_complex), dimension(:), allocatable  :: work
  real(c_float), dimension(:,:,:), pointer  :: r_array
  complex(c_float_complex), dimension(:,:,:), pointer :: c_array
  type(heffte_fft3d_r2c_cufft) :: fft
  !$acc declare present(input, output, work)

  !Interfaces of external subroutines in C
  interface
     subroutine fft_cuda_cpp(n1,n2,gn3,lksize,lkstart,r_data,c_data,dir,alloc) bind(c)
       use iso_c_binding
       implicit none
       integer(c_long),value :: n1,n2,gn3
       integer(c_long)       :: lksize, lkstart
       type(c_ptr), value    :: r_data, c_data
       integer(c_long),value :: dir
       integer(c_long),value :: alloc
     end subroutine fft_cuda_cpp

  end interface

  interface
     subroutine fft_alloc_cpp(n1,n2,gn3,slice_direction,r2c_direction,il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3) bind(c)
       use iso_c_binding
       implicit none
       integer(c_long),value :: n1,n2,gn3
       integer(c_int),value :: slice_direction, r2c_direction
       integer(c_int)       :: il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3
     end subroutine fft_alloc_cpp
  end interface

  
contains

  subroutine fft_cuda_alloc(n1,n2,gn3,lksize,lkstart)
    use iso_c_binding
    use mpi
    implicit none
    integer(ik), intent(in)  :: n1,n2,gn3
    integer(ik), intent(out) :: lksize, lkstart
!!!!!!!!!!!!!!!!!!!!!!!!
    !Allocate FFT structures
!!!!!!!!!!!!!!!!!!!!!!!!
    integer(c_long)          :: nn1,nn2,ngn3
    integer(c_long), target  :: nlksize,nlkstart
    integer(c_long) :: ndir
    integer(c_long) :: alloc, size_in, size_out, size_work
    integer(c_int)       :: slice_direction, r2c_direction
    integer(c_int) :: il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3
    
    r2c_direction=0
    slice_direction=2
    
    call fft_alloc_cpp(int(n1,c_long),int(n2,c_long),int(gn3,c_long),&
         &slice_direction,r2c_direction,il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3)

    lkstart = il3
    lksize = ih3-il3+1
    
    
    fft = heffte_fft3d_r2c_cufft(il1,il2,il3,ih1,ih2,ih3,&
         &ol1,ol2,ol3,oh1,oh2,oh3&
         &,r2c_direction,MPI_COMM_WORLD)

    size_in = fft%size_inbox()
    size_out = fft%size_outbox()
    size_work = fft%size_workspace()
    allocate(input(1:size_in))
    allocate(output(1:size_out))
    allocate(work(1:size_work))
    call c_f_pointer(c_loc(input),r_array,shape=(/n1, n2, lksize/))
    call c_f_pointer(c_loc(output),c_array,shape=(/dim1(n1)/2, n2, lksize/))
    !$acc enter data create(input(1:size_in), output(1:size_out), work(1:size_work))
    
    
    
    return

  end subroutine fft_cuda_alloc

  subroutine fft_cuda_dealloc
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

  end subroutine fft_cuda_dealloc

!!$  subroutine func(i,o,w)
!!$    use iso_c_binding
!!$    implicit none
!!$    real(c_float), dimension(:) :: i
!!$    complex(c_float_complex), dimension(:) :: o, w
!!$    call fft%forward(i,o,w)
!!$    return
!!$  end subroutine func
  
  subroutine fft_cuda_fourier(nn, gn3, dir, fu, nfs, nfe)
    use iso_c_binding
    use parameters, only : dim1
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    integer(ik), intent(in) :: gn3
    integer(ik), intent(in)       :: dir
    integer(ik), optional :: nfs, nfe
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(inout) :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: nnfs,nnfe,nfi
    integer(ik) :: i,j,k,iii
    integer(c_long)                     :: nn1,nn2,ngn3,nlksize,nlkstart,ndir,alloc
    real(rk) :: scale
    type(c_ptr) :: gpu_input_c, gpu_output_c, gpu_workspace_c
    integer(c_long) :: size_in, size_out, size_work
    type(c_ptr) :: gpu_input, gpu_output, gpu_workspace
    type(c_funptr) :: func_ptr
    complex(ck) :: ctmp
!!$    real(c_float), dimension(:,:,:), pointer, contiguous :: r_array
!!$    complex(c_float_complex), dimension(:,:,:), pointer, contiguous :: c_array
!!$      
    scale = 1.0_rk / real(nn(1) * nn(2) * gn3)
    nnfs=1
    nnfe=nn(4)
    if(present(nfs)) nnfs=nfs
    if(present(nfe)) nnfe=nfe

    size_in=fft%size_inbox()
    size_out=fft%size_outbox()
    size_work=fft%size_workspace()
!!$
!!$    call c_f_pointer(c_loc(input),r_array,shape=[n1,n2,lksize])
!!$    call c_f_pointer(c_loc(output),c_array,shape=[dim1(n1)/2,n2,lksize])

    !$acc data present(input(1:size_in), output(1:size_out), work(1:size_work))
    if(dir == 1) then
       do nfi=nnfs,nnfe

          
          !$omp parallel do 
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             r_array(i,j,k)=fu(i,j,k,nfi)
          end do; end do ; end do 
          !$omp end parallel do
          
          
         
          !$acc update device(input(1:size_in)) 
          !$acc host_data use_device(input(1:size_in), output(1:size_out), work(1:size_work))
          call fft%forward(input,output,work,scale_cufft_none)
          !$acc end host_data
          !$acc update host(output(1:size_out))
         
          
          
          
          !$omp parallel do private(iii, ctmp)
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
             iii = (i-1) / 2 + 1
             ctmp = c_array(iii,j,k)
             fu(i,j,k,nfi) = real(ctmp)
             fu(i+1,j,k,nfi) = aimag(ctmp)
          end do; end do ; end do
          !$omp end parallel do
          
          
       end do
       

    else 
       
       do nfi=nnfs,nnfe

          
          !$omp parallel do private(iii)
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
             iii = (i - 1) / 2  + 1
             c_array(iii,j,k) = cmplx(fu(i,j,k,nfi), fu(i+1,j,k,nfi))
          end do; end do ; end do
          !$omp end parallel do


          !$acc update device(output(1:size_out)) 
          !$acc host_data use_device(input(1:size_in), output(1:size_out), work(1:size_work))
          call fft%backward(output,input,work, scale_cufft_full)
          !$acc end host_data
          !$acc update host(input(1:size_in))
         
          
         
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             fu(i,j,k,nfi) = r_array(i,j,k)
          end do; end do ; end do
          !$omp end parallel do

          
          
       end do
       
    end if
    !$acc end data


    return

  end subroutine fft_cuda_fourier




end module fft_cuda
