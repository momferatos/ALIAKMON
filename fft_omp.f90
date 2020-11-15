!!$     ___ __                                       
!!$ (  / _ \\ \        /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!$Copyright (c) 2009-2017 Georgios Momferratos
!$
!$This program is free software; you can redistribute it and/or modify
!$it under the terms of the GNU General Public License as published by
!$the Free Software Foundation; either version 2 of the License, or
!$(at your option) any later version.
!$ 
!$This program is distributed in the hope that it will be useful,
!$but WITHOUT ANY WARRANTY; without even the implied warranty of
!$MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!$GNU General Public License for more details.
!$
!$You should have received a copy of the GNU General Public License
!$along with this program; if not, write to the Free Software
!$Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!$
!$
!$

module fft_omp
#ifdef _OPENMP
  use types
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
#ifdef _DOUBLE_
  real(c_double), dimension(:,:,:), pointer :: array
  complex(c_double_complex), dimension(:,:,:), pointer   :: carray
#else
  real(c_float), dimension(:,:,:), pointer :: array
  complex(c_float_complex), dimension(:,:,:), pointer :: carray
#endif
  type(c_ptr) :: plan, iplan, p
contains


  subroutine fft_omp_alloc(n1, n2, n3)
    use omp_lib
    use parameters, only: dim1
    implicit none
    integer(ik) :: n1, n2, n3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(c_int) :: nt, error
    integer(c_int) :: nn1, nn2, nn3
    !$omp parallel    
    nt = omp_get_num_threads()
    !$omp end parallel
    nn1 = int(n1, c_int)
    nn2 = int(n2, c_int)
    nn3 = int(n3, c_int)

#ifdef _DOUBLE_    
    p = fftw_alloc_real(int(dim1(n1) * n2 * n3, C_SIZE_T))
#else
    p = fftwf_alloc_complex(int(dim1(n1) * n2 * n3, C_SIZE_T))
#endif

    call c_f_pointer(p, array, shape=[ dim1(n1), n2, n3 ])
    call c_f_pointer(p, carray, shape=[ dim1(n1) / 2, n2, n3 ])

#ifdef _DOUBLE_
    error = fftw_init_threads()
    call fftw_plan_with_nthreads(nt)
    plan = fftw_plan_dft_r2c_3d(nn3, nn2, nn1, array,  carray,FFTW_MEASURE)
    iplan = fftw_plan_dft_c2r_3d(nn3, nn2, nn1, carray, array,FFTW_MEASURE)
#else
    error = fftwf_init_threads()
    call fftwf_plan_with_nthreads(nt)
    plan = fftwf_plan_dft_r2c_3d(nn3, nn2, nn1, array,  carray,FFTW_MEASURE);
    iplan = fftwf_plan_dft_c2r_3d(nn3, nn2, nn1, carray, array,FFTW_MEASURE);
#endif

    return

  end subroutine fft_omp_alloc

  subroutine fft_omp_dealloc
    implicit none

    call fftw_free(p)

#ifdef _DOUBLE_
    call fftw_destroy_plan(plan)
    call fftw_destroy_plan(iplan)
#else
    call fftwf_destroy_plan(plan)
    call fftwf_destroy_plan(iplan)
#endif

    return

  end subroutine fft_omp_dealloc

  subroutine fft_omp_fourier(nn, dir, fu, nfs, nfe)
    use parameters, only : dim1
    implicit none
    integer(ik), dimension(1:4), intent(in) :: nn
    integer(ik), intent(in)       :: dir
    integer(ik), optional :: nfs,nfe
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(inout) :: fu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(rk)                  :: scale
    integer(ik) :: nnfs,nnfe,nfi
    integer(ik) :: i,j,k

    scale = 1.0_rk / real(nn(1)*nn(2)*nn(3), rk)

    nnfs=1
    nnfe=nn(4)
    if(present(nfs)) nnfs=nfs
    if(present(nfe)) nnfe=nfe


    do nfi=nnfs,nnfe
       if(dir == 1) then
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
             array(i,j,k) = fu(i,j,k, nfi)
          end do; end do ; end do 
          !$omp end parallel do
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=nn(1)+1,dim1(nn(1))
             array(i,j,k) = 0.0_rk
          end do; end do ; end do
          !$omp end parallel do
       else 
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
             array(i,j,k) = fu(i,j,k, nfi)
          end do; end do ; end do
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
          end do; end do ; end do
          !$omp end parallel do
       else
          !$omp parallel do
          do k=1,nn(3) ; do j=1,nn(2) ; do i=1,dim1(nn(1))
             fu(i,j,k,nfi) = array(i,j,k)
          end do; end do ; end do
          !$omp end parallel do
       end if

    end do

    return

  end subroutine fft_omp_fourier

#endif
end module fft_omp

