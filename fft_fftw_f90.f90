!!$     ___ __                                       
!!$ (  / _ \\ \       /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!!$   Copyright (c) 2009-2020 George Momferatos

module fft_fftw
  use mpi
  use mpivars
  use heffte_fftw
  use types
  use parameters
  use iso_c_binding
  implicit none
  !
  ! Heffte FFTW Fortran wrapper
  !

  real(c_float), dimension(:,:,:), pointer :: r_array
  complex(c_float_complex), dimension(:,:,:), pointer :: c_array
  real(c_float), dimension(:), allocatable, target, private :: input
  complex(c_float_complex), dimension(:), allocatable, target, private  :: output
  complex(c_float_complex), dimension(:), allocatable, private :: work
  !$acc declare create(input, output, work)
  type(heffte_fft3d_r2c_fftw) :: fft

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
     subroutine fft_cuda_alloc_cpp(n1,n2,gn3,slice_direction,r2c_direction,il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3) bind(c)
       use iso_c_binding
       implicit none
       integer(c_long),value :: n1,n2,gn3
       integer(c_int),value :: slice_direction, r2c_direction
       integer(c_int)       :: il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3
     end subroutine fft_cuda_alloc_cpp
  end interface

! void fft_cuda_load_cpp(float *input_f, float *gpu_input_f, float *gpu_output_f, float *gpu_workspace_f, 

  interface
     subroutine fft_cuda_load_cpp(input, gpu_input, gpu_output, gpu_workspace, size_inbox, size_outbox, size_workspace,fft) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value    :: input, gpu_input, gpu_output, gpu_workspace
       integer(c_long), value :: size_inbox, size_outbox, size_workspace
       type(c_funptr) :: fft
     end subroutine fft_cuda_load_cpp
  end interface

  interface
     subroutine fft_cuda_unload_cpp(output, gpu_output, size_outbox) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value    :: output, gpu_output
       integer(c_long), value :: size_outbox
     end subroutine fft_cuda_unload_cpp
  end interface
  
contains

  subroutine fft_fftw_alloc(n1,n2,gn3,lksize,lkstart)
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
    integer(c_long) :: alloc
    integer(c_int)       :: slice_direction, r2c_direction
    integer(c_int) :: il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3
    
    r2c_direction=0
    slice_direction=2
    
    call fft_cuda_alloc_cpp(int(n1,c_long),int(n2,c_long),int(gn3,c_long),&
         &slice_direction,r2c_direction,il1,il2,il3,ih1,ih2,ih3,ol1,ol2,ol3,oh1,oh2,oh3)

    lkstart = il3
    lksize = ih3-il3+1
    
    
    fft = heffte_fft3d_r2c_fftw(il1,il2,il3,ih1,ih2,ih3,&
         &ol1,ol2,ol3,oh1,oh2,oh3&
         &,r2c_direction,MPI_COMM_WORLD)

    allocate(input(fft%size_inbox()))
    allocate(output(fft%size_outbox()))
    allocate(work(fft%size_workspace()))

    !allocate(pinput(fft%size_inbox()))
    !allocate(poutput(fft%size_outbox()))
    !input=0.0

    !stop
    !allocate(r_array(1:dim1(n1),1:n2,1:lksize))
    !allocate(c_array(1:n1,1:n2,1:lksize))

    
    call c_f_pointer(c_loc(input),r_array,shape=(/n1, n2, lksize/))
    call c_f_pointer(c_loc(output),c_array,shape=(/dim1(n1)/2, n2, lksize/))
    
    return

  end subroutine fft_fftw_alloc

  subroutine fft_fftw_dealloc
    use iso_c_binding
    !
    !Deallocate FFT structures
    !
    deallocate(input)
    deallocate(output)
    deallocate(work)
    
    return

  end subroutine fft_fftw_dealloc

!!$  subroutine func(i,o,w)
!!$    use iso_c_binding
!!$    implicit none
!!$    real(c_float), dimension(:) :: i
!!$    complex(c_float_complex), dimension(:) :: o, w
!!$    call fft%forward(i,o,w)
!!$    return
!!$  end subroutine func
  
  subroutine fft_fftw_fourier(nn, gn3, dir, fu, nfs, nfe)
    use iso_c_binding
    use parameters, only : dim1
    use iso_c_binding
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
    
    scale = 1.0_rk / real(nn(1) * nn(2) * gn3)
    !    func_ptr = c_funloc(func)
    nnfs=1
    nnfe=nn(4)
    if(present(nfs)) nnfs=nfs
    if(present(nfe)) nnfe=nfe

    size_in=fft%size_inbox()
    size_out=fft%size_outbox()
    size_work=fft%size_workspace()


    
    if(dir == 1) then
       do nfi=nnfs,nnfe

         !!$omp parallel do 
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             r_array(i,j,k)=fu(i,j,k,nfi)
          end do; end do ; end do 
          !$omp end parallel do

          
          !$acc data  copyin(input) copyout(output) create(work)
          !$acc host_data use_device(input, output, work)
          call fft%forward(input,output,work, scale_fftw_none)
          !$acc end host_data
          !$acc end data

          
          !!omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
             iii = (i-1) / 2 + 1
             ctmp = c_array(iii,j,k)
             fu(i,j,k,nfi) = real(ctmp)
             fu(i+1,j,k,nfi) = aimag(ctmp)
          end do; end do ; end do
          !!omp end parallel do
       end do
       !$acc exit data delete(work)

    else 
       !$acc enter data create(work)
       do nfi=nnfs,nnfe


          !!omp parallel do 
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))-1,2
             iii = (i - 1) / 2  + 1
             c_array(iii,j,k) = cmplx(fu(i,j,k,nfi), fu(i+1,j,k,nfi))
          end do; end do ; end do
          !!omp end parallel do
                    
          !$acc data  copyin(output) copyout(output) create(work)
          !$acc host_data use_device(input, output, work)
          call fft%backward(output,input,work, scale_fftw_full)
          !$acc end host_data
          !$acc end data 
          
          
          !!omp parallel do 
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             fu(i,j,k,nfi) = r_array(i,j,k)
          end do; end do ; end do
          !!omp end parallel do

          
       end do
       
    end if
    


    return

  end subroutine fft_fftw_fourier




end module fft_fftw
