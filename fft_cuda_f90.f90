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
  use types
  use parameters
  use iso_c_binding
  implicit none
  !
  ! Heffte CUDA Fortran wrapper
  !
  real(rks), dimension(:,:,:), allocatable :: r_array
  real(rks), dimension(:,:,:), allocatable :: c_array


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
  
contains

  subroutine fft_cuda_alloc(n1,n2,gn3,lksize,lkstart)
    use iso_c_binding
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

    !Convert to C types
    nn1=int(n1,c_long)
    nn2=int(n2,c_long)
    ngn3=int(gn3,c_long)
    nlksize=int(lksize,c_long)
    nlkstart=int(lkstart,c_long)
    ndir=int(1,c_long)
    alloc=int(1,c_long)
    
    !Call C function
    call fft_cuda_cpp(nn1,nn2,ngn3,nlksize,nlkstart,c_loc(r_array),c_loc(c_array),ndir,alloc) 
    
    !Get back local slice size and starting index
    lksize=int(nlksize,ik)
    lkstart=int(nlkstart,ik)

    allocate(r_array(1:n1,1:n2,1:lksize))
    allocate(c_array(1:dim1(n1),1:n2,1:lksize))

    !allocate(r_array(1:lksize,1:nn2,1:nn1))
    !allocate(c_array(1:lksize,1:nn2,1:dim1(nn1)))

    return

  end subroutine fft_cuda_alloc

  subroutine fft_cuda_dealloc
    use iso_c_binding
    !
    !Deallocate FFT structures
    !
    deallocate(r_array)
    deallocate(c_array)
    return
    
  end subroutine fft_cuda_dealloc

  subroutine fft_cuda_fourier(nn, gn3, dir, fu, nfs, nfe)
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
    integer(ik) :: i,j,k
    integer(c_long)                     :: nn1,nn2,ngn3,nlksize,nlkstart,ndir,alloc
    real(rk) :: scale

    scale = 1.0_rk !/ real(n1*n2*gn3,rk)

    !Convert to C types
    nn1=int(n1,c_long)
    nn2=int(n2,c_long)
    ngn3=int(gn3,c_long)
    ndir=int(dir,c_long)
    alloc=int(0,c_long)

    
    nnfs=1
    nnfe=nn(4)
    if(present(nfs)) nnfs=nfs
    if(present(nfe)) nnfe=nfe

    do nfi=nnfs,nnfe
       if(dir == 1) then
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             r_array(i,j,k) = fu(i,j,k, nfi)
          end do; end do ; end do 
          !$omp end parallel do
       else 
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
             c_array(i,j,k) = fu(i,j,k, nfi)
          end do; end do ; end do
          !$omp end parallel do
       end if

       call fft_cuda_cpp(nn1,nn2,ngn3,nlksize,nlkstart,c_loc(r_array),c_loc(c_array),int(dir,c_long),alloc)

       if(dir == 1) then
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
             fu(i,j,k,nfi) = c_array(i,j,k)
          end do; end do ; end do
          !$omp end parallel do
       else
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             fu(i,j,k,nfi) = r_array(i,j,k)
          end do; end do ; end do
          !$omp end parallel do
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=nn(1)+1,dim1(nn(1))
             !fu(i,j,k,nfi) = 0.0_rk
          end do; end do ; end do
          !$omp end parallel do
       end if


    end do

    return

  end subroutine fft_cuda_fourier


end module fft_cuda
