module vtk
  use types
  use parameters
  implicit none

contains

  subroutine output_scalar_vtk_2d_file_ascii(nn,u,filename,dataname,&
       &comment, keys)
    implicit none
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(ik), intent(IN), dimension(1:4)                         :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)   :: u
    character(*)                                    :: filename, dataname,&
         &comment
    character(len=64), dimension(:) :: keys
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(ik) :: i, j, k, n1, n2, n3
    real(rk) :: umin=0.0_rk, umax=1.0_rk
    character(128) :: fmt
    integer :: vtk_file=17
    real(rk) :: LBOX = 2.0 * PI

    
    if(rks==dp) then
       fmt='(t1,e22.16)'
    else
       fmt='(t1,e16.8)'
    end if

    n1=nn(1)
    n2=nn(2)
    n3=nn(3)

    open(vtk_file,file=filename,form='formatted',action='write',err=100)
    write(vtk_file, '(t1,a)') '# vtk DataFile Version 2.0'
    write(vtk_file, '(t1,a)') comment
    write(vtk_file, '(t1,a)') 'ASCII'

    write(vtk_file, '(t1,a)') 'DATASET STRUCTURED_GRID'
    write(vtk_file, '(t1,a,3i8)') 'DIMENSIONS ', n1, n2, 1
    write(vtk_file, '(t1,a,i15,a)') 'POINTS ', n1*n2, ' float'



    do j=1,n2
       do i=1,n1
          write(vtk_file, '(t1,2f15.6,f3.0)') real((i-1)*LBOX/(n1-1)), &
               &real((j-1)*LBOX/(n2-1)), 0.
       end do
    end do


    write(vtk_file, '(t1,a,i15)') 'POINT_DATA ', n1*n2
    if(rks==dp) then
       write(vtk_file, '(t1,3a)') 'SCALARS ', dataname, ' double 1'
    else
       write(vtk_file, '(t1,3a)') 'SCALARS ', dataname, ' float 1'
    end if

    write(vtk_file, '(t1,a)') 'LOOKUP_TABLE default'

    do j=1,n2
       do i=1,n1
          write(vtk_file, fmt) u(i,j,n3/2,1)
       end do
    end do


    close(vtk_file,status='keep')

    return

100 stop 'Error. Cannot open file.'


  end subroutine output_scalar_vtk_2d_file_ascii

  subroutine output_scalar_vtk_2d_file(nn,u,filename,tt,keys)
    implicit none
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(ik), dimension(1:4), intent(IN)                    :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN)   :: u
    character(*)                                    :: filename
    real(rk) :: tt
    character(len=*), dimension(1:) :: keys
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integer(ik) :: i, j, k, n1, n2, n3, nf
    real(rk) :: umin=0.0_rk, umax=1.0_rk
    character(128) :: fmt
    integer :: vtk_file=17
    real(rk) :: LBOX = 2.0 * PI
    character(len=1024) :: comment
    character(len=128) :: dataname
    n1=nn(1)
    n2=nn(2)
    n3=nn(3)

    
    write(comment, '(a,e19.7)') 'time: ', tt

    if(rks==dp) then
       fmt='(t1,e22.16)'
    else
       fmt='(t1,e16.8)'
    end if



    open(vtk_file,file=filename,form='formatted',action='write',err=100)
    write(vtk_file, '(t1,a)') '# vtk DataFile Version 2.0'
    write(vtk_file, '(t1,a)') trim(comment)
    write(vtk_file, '(t1,a)') 'BINARY'

    write(vtk_file, '(t1,a)') 'DATASET STRUCTURED_GRID'
    write(vtk_file, '(t1,a,3i8)') 'DIMENSIONS ', n1, n2, 1
    write(vtk_file, '(t1,a,i15,a)') 'POINTS ', n1*n2, ' float '
    close(vtk_file, status='keep')

    open(vtk_file, file=filename,form='unformatted',action='write',&
         &position='append', access='stream', convert='big_endian', err=100) 
    do j=1,n2
       do i=1,n1
          write(vtk_file) real((i-1)*LBOX/(n1-1), sp), &
               &real((j-1)*LBOX/(n2-1), sp), real(0., sp)
       end do
    end do
    close(vtk_file, status='keep')
    open(vtk_file,file=filename,form='formatted',action='write',&
         &position='append', err=100)
    write(vtk_file, *)
    write(vtk_file, '(t1,a,i15)') 'POINT_DATA ', n1*n2
    close(vtk_file,status='keep')
    do nf=nu1,nu3+nscl
       open(vtk_file,file=filename,form='formatted',action='write',&
            &position='append', err=100)
       write(vtk_file, '(t1,3a)') 'SCALARS ', trim(keys(nf)), ' float 1'


       write(vtk_file, '(t1,a)') 'LOOKUP_TABLE default'
       close(vtk_file, status='keep')

       open(vtk_file, file=filename,form='unformatted',action='write',&
            &position='append', access='stream', convert='big_endian', err=100) 
       do j=1,n2
          do i=1,n1
             write(vtk_file) real(u(i,j,n3/2,nf), sp)
          end do
       end do
       close(vtk_file,status='keep')
       open(vtk_file,file=filename,form='formatted',action='write',&
            &position='append', err=100)
       write(vtk_file,*)
       close(vtk_file,status='keep')
    end do

    return

100 stop 'Error. Cannot open file.'


  end subroutine output_scalar_vtk_2d_file
  
end module vtk
