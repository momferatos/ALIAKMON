!!$     ___ __                                       
!!$ (  / _ \\ \        /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!$Copyright (c) 2009-2020 Georgios Momferatos
module hdf5_aliakmon
  use types
#ifdef _MPI_
  use mpi
#endif
  use mpivars,only:mpirank,mpisize,mpiroot,lkstart,ljstart
  use parameters
  use data
  use hdf5
  implicit none

  interface write_hdf5_file
     module procedure write_hdf5_compute_file
     module procedure write_hdf5_pp_file
     module procedure write_hdf5_slice_file
  end interface write_hdf5_file


  integer(hid_t) :: file_id       ! File identifier 
  integer(hid_t) :: dset_id       ! Dataset identifier 
  integer(hid_t) :: filespace     ! Dataspace identifier in file 
  integer(hid_t) :: memspace      ! Dataspace identifier in memory
  integer(hid_t) :: plist_id      ! Property list identifier
  integer(hid_t) :: dcpl
  

  integer(hsize_t), dimension(1:4) :: dimsf ! Dataset dimensions.
  !     INTEGER, DIMENSION(7) :: dimsfi = (/5,8,0,0,0,0,0/)
  integer(hsize_t), dimension(1:4) :: dimsfi 

  integer(hsize_t), dimension(1:4) :: icount  
  integer(hsize_t), dimension(1:4) :: offset
  
  integer :: rank ! Dataset rank 

  integer :: error, error_n  ! Error flags

  logical :: gzip_avail
  integer :: filter_info, filter_info_both

  integer(hsize_t), dimension(4) :: vecchunk = [3, 32, 32, 32]
  integer(hsize_t), dimension(3) :: sclchunk = [32, 32, 32]
  integer(hsize_t), dimension(2) :: slicechunk = [8, 8]
  
  real(rks), dimension(:,:), allocatable :: h5_slice_data
  real(rks), dimension(:,:,:), allocatable :: h5_scalar_data
  real(rks), dimension(:,:,:,:), allocatable :: h5_vector_data

contains

  subroutine add_timestamp(filename, time)
    implicit none
    character(len=*), intent(IN) :: filename
    real(rk), intent(IN) :: time
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t) :: hdf5_file
    character(len = 4), parameter :: dsetname = 'time'
    integer :: trank
    integer(hid_t) :: dset_id, dcpl
    integer        :: error
    integer(hsize_t), dimension(1) :: data_dims
    real(8), dimension(1) :: time_arr

    time_arr(1) = time

    ! set dataset rank                                                          
    trank = 1
    ! Set dataset dimensions for data                                           
    data_dims(1) = 1

    call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, hdf5_file, error)

    ! Create the data space for the  dataset.                                   
    call h5screate_simple_f(trank, data_dims, filespace, error)
    ! Create the dataset                                                        
    call h5dcreate_f(hdf5_file, trim(dsetname), H5T_NATIVE_DOUBLE, filespace, &
         &dset_id, error)

    ! write the dataset                                                         

    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time_arr, data_dims, error)

    ! cleanup
    
    call h5sclose_f(filespace, error)

    call h5dclose_f(dset_id, error)
    
    call h5fclose_f(hdf5_file, error)
    
    return
    
  end subroutine add_timestamp
  
  subroutine write_hdf5_vector_dataset(dataset_name, dataset, nn)
    implicit none
    integer(ik), dimension(1:4), intent(in)   :: nn          
    character(len=*), intent(IN)            :: dataset_name
    real(rks), dimension(1:3,1:nn(1),1:nn(2),1:nn(3)), intent(IN) :: dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rank = 4

    dimsf(1) = 3
    dimsf(2) = nn(1)
    dimsf(3) = gn2
    dimsf(4) = gn3

    icount(1) = 3
    icount(2) = nn(1)
    icount(3) = nn(2)
    icount(4) = nn(3)

    offset(1) = 0
    offset(2) = 0
    offset(3) = ljstart-1
    offset(4) = lkstart-1

    

    !
    ! Create the data space for the  dataset. 
    !
    call h5screate_simple_f(rank, dimsf, filespace, error)

    if(gzip_avail .and. COMPRESSION_LEVEL /= 0) then
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, error)
       call h5pset_deflate_f(dcpl, COMPRESSION_LEVEL, error)
       call h5pset_chunk_f(dcpl, rank, vecchunk, error)
    end if
    !
    ! Create the dataset with default properties.
    !

    if(rks==sp) then
       call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_REAL, filespace, &
            dset_id, error, dcpl_id=dcpl)
    else
       call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_DOUBLE, filespace, &
            dset_id, error, dcpl_id=dcpl)
    end if
    
#ifdef _MPI_
    call h5sclose_f(filespace, error)
    !
    ! Each process defines dataset in memory and writes it to the hyperslab
    ! in the file. 
    !
    call h5screate_simple_f(rank, icount, memspace, error) 
    ! 
    ! Select hyperslab in the file.
    !
    call h5dget_space_f(dset_id, filespace, error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
         &icount, error)

    !
    ! Create property list for collective dataset write
    !
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    !
    ! Write the dataset collectively. 
    !
    if(rks==sp) then
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, dataset, dimsfi, error, &
            file_space_id = filespace, mem_space_id = memspace, &
            &xfer_prp = plist_id)
    else
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dimsfi, error, &
            file_space_id = filespace, mem_space_id = memspace, &
            &xfer_prp = plist_id)
    end if

    !
    ! Close dataspaces.
    !
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)

#else
    if(rks==sp) then
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, dataset, dimsfi, error)
    else
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dimsfi, error)
    end if
#endif
    !
    ! Write the dataset independently. 
    !
    !    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dimsfi, error, &
    !                     file_space_id = filespace, mem_space_id = memspace)
    !


    !
    ! Close the dataset
    !
    call h5dclose_f(dset_id, error)

    return

  end subroutine write_hdf5_vector_dataset

  subroutine write_hdf5_scalar_dataset(dataset_name, dataset, nn)
    implicit none
    integer(ik), dimension(1:4), intent(in)   :: nn          
    character(len=*), intent(IN)            :: dataset_name
    real(rks), dimension(1:nn(1),1:nn(2),1:nn(3)), intent(IN) :: dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rank = 3
    dimsf(1)=nn(1)
    dimsf(2)=gn2
    dimsf(3)=gn3

    icount(1) = dimsf(1)
    icount(2) = nn(2)
    icount(3) = nn(3)

    offset(1) = 0
    offset(2) = ljstart-1
    offset(3) = lkstart-1


    !
    ! Create the data space for the  dataset. 
    !
    call h5screate_simple_f(rank, dimsf, filespace, error)

    if(gzip_avail .and. COMPRESSION_LEVEL /= 0) then
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, error)
       call h5pset_deflate_f(dcpl, COMPRESSION_LEVEL, error)
       call h5pset_chunk_f(dcpl, rank, sclchunk, error)
    end if
       
    !
    ! Create the dataset with default properties.
    !

    if(rks==sp) then
       call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_REAL, filespace, &
            dset_id, error, dcpl_id=dcpl)
    else
       call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_DOUBLE, filespace, &
            dset_id, error, dcpl_id=dcpl)
    end if
#ifdef _MPI_
    call h5sclose_f(filespace, error)
    !
    ! Each process defines dataset in memory and writes it to the hyperslab
    ! in the file. 
    !

    call h5screate_simple_f(rank, icount, memspace, error) 
    ! 
    ! Select hyperslab in the file.
    !
    call h5dget_space_f(dset_id, filespace, error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
         &icount, error)

    !
    ! Create property list for collective dataset write
    !
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    !
    ! Write the dataset collectively. 
    !
    if(rks==sp) then
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, dataset, dimsfi, error, &
            file_space_id = filespace, mem_space_id = memspace, &
            &xfer_prp = plist_id)
    else
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dimsfi, error, &
            file_space_id = filespace, mem_space_id = memspace, &
            &xfer_prp = plist_id)
    end if

    !
    ! Close dataspaces.
    !
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)

#else
    if(rks==sp) then
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, dataset, dimsfi, error)
    else
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dimsfi, error)
    end if
#endif
    !
    ! Write the dataset independently. 
    !
    !    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dimsfi, error, &
    !                     file_space_id = filespace, mem_space_id = memspace)
    !


    !
    ! Close the dataset
    !
    call h5dclose_f(dset_id, error)

    return

  end subroutine write_hdf5_scalar_dataset

  subroutine write_hdf5_slice_dataset(dataset_name, dataset, nn)
    implicit none
    integer(ik), dimension(1:4), intent(in)   :: nn          
    character(len=*), intent(IN)            :: dataset_name
    real(rks), dimension(1:nn(2),1:nn(3)), intent(IN) :: dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rank = 2
    dimsf(1)=gn2
    dimsf(2)=gn3

    icount(1) = nn(2)
    icount(2) = nn(3)

    offset(1) = ljstart-1
    offset(2) = lkstart-1


    !
    ! Create the data space for the  dataset. 
    !
    call h5screate_simple_f(rank, dimsf, filespace, error)

    if(gzip_avail .and. COMPRESSION_LEVEL /= 0) then
       call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, error)
       call h5pset_deflate_f(dcpl, COMPRESSION_LEVEL, error)
       call h5pset_chunk_f(dcpl, rank, slicechunk, error)
    end if


    !
    ! Create the dataset with default properties.
    !

    if(rks==sp) then
       call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_REAL, filespace, &
            dset_id, error, dcpl_id=dcpl)
    else
       call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_DOUBLE, filespace, &
            dset_id, error, dcpl_id=dcpl)
    end if
#ifdef _MPI_
    call h5sclose_f(filespace, error)
    !
    ! Each process defines dataset in memory and writes it to the hyperslab
    ! in the file. 
    !

    call h5screate_simple_f(rank, icount, memspace, error) 
    ! 
    ! Select hyperslab in the file.
    !
    call h5dget_space_f(dset_id, filespace, error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
         &icount, error)

    !
    ! Create property list for collective dataset write
    !
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

    !
    ! Write the dataset collectively. 
    !
    if(rks==sp) then
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, dataset, dimsfi, error, &
            file_space_id = filespace, mem_space_id = memspace, &
            &xfer_prp = plist_id)
    else
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dimsfi, error, &
            file_space_id = filespace, mem_space_id = memspace, &
            &xfer_prp = plist_id)
    end if

    !
    ! Close dataspaces.
    !
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)

#else
    if(rks==sp) then
       call h5dwrite_f(dset_id, H5T_NATIVE_REAL, dataset, dimsfi, error)
    else
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dimsfi, error)
    end if
#endif
    !
    ! Write the dataset independently. 
    !
    !    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dimsfi, error, &
    !                     file_space_id = filespace, mem_space_id = memspace)
    !


    !
    ! Close the dataset
    !
    call h5dclose_f(dset_id, error)

    return

  end subroutine write_hdf5_slice_dataset

  subroutine write_hdf5_compute_file(nn,u,time,nfile)
    use data, only: ga
    use fvdom, only: calcqr
    implicit none
    integer(ik), intent(IN), dimension(1:4)                  :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3), 1:nn(4)), intent(IN) :: u
    real(rk), intent(IN)                                   :: time
    integer(ik), intent(IN)                                :: nfile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Writes HDF5 file in parallel !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(LEN=256)              :: filename ! File name
    character(len=64), dimension(1:nn(4)) :: datanames
    logical, dimension(1:nn(4) + 2) :: data_is_vector
    integer(ik)                     :: ndatanames
    integer(ik) :: i,j,k,l
    write(filename,'(a,i6.6,a)') 'output.',nfile,'.h5'

!!$    dimsf(1)=nn(1)
!!$    dimsf(2)=gn2
!!$    dimsf(3)=gn3
!!$    dimsf(4)=3
!!$
!!$    dimsfi(:)=dimsf(:)

    !
    ! Initialize FORTRAN predefined datatypes
    !
    call h5open_f(error) 
    call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, gzip_avail, error)

    if(gzip_avail .and. COMPRESSION_LEVEL /= 0) then
       call h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, filter_info, error)
       filter_info_both=ior(H5Z_FILTER_ENCODE_ENABLED_F,&
            &H5Z_FILTER_DECODE_ENABLED_F)
       if(filter_info .ne. filter_info_both) gzip_avail = .false.
    end if

#ifdef _MPI_
    ! 
    ! Setup file access property list with parallel I/O access.
    !
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
#endif

    !
    ! Create the file collectively.
    ! 
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error,&
         & access_prp = plist_id)
    call h5pclose_f(plist_id, error)

    allocate(h5_scalar_data(1:nn(1),1:nn(2),1:nn(3)))
    allocate(h5_vector_data(1:3,1:nn(1),1:nn(2),1:nn(3)))

    data_is_vector(:) = .false.
    ndatanames = 0
    !$omp parallel do
    do l=1,3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       h5_vector_data(l,i,j,k)=u(i,j,k,l)
    end do; end do ; end do ; end do
    !$omp end parallel do
    call write_hdf5_vector_dataset('/u',h5_vector_data,nn)
    ndatanames = ndatanames + 1
    datanames(ndatanames) = 'u'
    data_is_vector(ndatanames) = .true.
    !write passive scalar
    if(PASSIVE_SCALAR) then
       h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))=u(1:nn(1),1:nn(2),1:nn(3),nsclf)
       call write_hdf5_scalar_dataset('/scl',h5_scalar_data,nn)
       ndatanames = ndatanames + 1
       datanames(ndatanames) = 'scl'
    end if

    if(RADIATION) then
       call calcqr
       h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))=ga(1:nn(1),1:nn(2),1:nn(3))
       call write_hdf5_scalar_dataset('/G', h5_scalar_data, nn)
       ndatanames = ndatanames + 1
       datanames(ndatanames) = 'G'
       !$omp parallel do
       do l=1,3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          h5_vector_data(l,i,j,k)=qr(i,j,k,l)
       end do; end do ; end do ; end do
       !$omp end parallel do
       call write_hdf5_vector_dataset('/q', h5_vector_data, nn)
       ndatanames = ndatanames + 1
       datanames(ndatanames) = 'q'
       data_is_vector(ndatanames) = .true.
    end if


!!$    !write mhd
!!$    if(MHD) then
!!$       h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))=u(1:nn(1),1:nn(2),1:nn(3),nb1)
!!$       call write_hdf5_scalar_dataset('/b1',h5_scalar_data,nn)
!!$
!!$       h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))=u(1:nn(1),1:nn(2),1:nn(3),nb1)
!!$       call write_hdf5_scalar_dataset('/b2',h5_scalar_data,nn)
!!$
!!$       h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))=u(1:nn(1),1:nn(2),1:nn(3),nb1)
!!$       call write_hdf5_scalar_dataset('/b3',h5_scalar_data,nn)
!!$    end if

    ! Deallocate data buffer.
    !
    deallocate(h5_scalar_data)
    deallocate(h5_vector_data)


    call h5pclose_f(plist_id, error)

    !
    ! Close the file.
    !
    call h5fclose_f(file_id, error)

    if(mpirank==mpiroot) then
       call write_xdmf_file(nn(1),nn(2),ndatanames,&
            &datanames, data_is_vector, filename,nfile)
       call add_timestamp(filename, time)
    end if

    !
    ! Close FORTRAN predefined datatypes.
    !
    call h5close_f(error)


    
!!$    datanames(nu2)='u2'
!!$    datanames(nu3)='u3'
!!$    if(PASSIVE_SCALAR) then
!!$       datanames(nsclf:nscll)='scl'
!!$    end if
!!$    if(MHD) then
!!$       datanames(nb1)='b1'
!!$       datanames(nb2)='b2'
!!$       datanames(nb3)='b3'
!!$    end if



    return

  end subroutine write_hdf5_compute_file

  subroutine write_hdf5_slice_file(nn,nfields,slice,datanames,time,nfile)
    use data, only: ga
    use fvdom, only: calcqr
    implicit none
    integer(ik), intent(IN), dimension(1:4)                :: nn
    integer(ik), intent(IN)                                :: nfields
    real(rks), dimension(1:nfields,1:nn(2), 1:nn(3)), intent(IN) :: slice
    character(len=64), dimension(1:nfields) :: datanames
    real(rk), intent(IN)                                   :: time
    integer(ik), intent(IN)                                :: nfile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Writes HDF5 file in parallel !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(LEN=256)              :: filename ! File name
    integer(ik) :: i,j,k,l
    integer(ik) :: nf

    write(filename,'(a,i6.6,a)') 'slice.',nfile,'.h5'


    !
    ! Initialize FORTRAN predefined datatypes
    !
    call h5open_f(error)
    
    call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, gzip_avail, error)

    if(gzip_avail .and. COMPRESSION_LEVEL /= 0) then
       call h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, filter_info, error)
       filter_info_both=ior(H5Z_FILTER_ENCODE_ENABLED_F,&
            &H5Z_FILTER_DECODE_ENABLED_F)
       if(filter_info .ne. filter_info_both) gzip_avail = .false.
    end if
    
#ifdef _MPI_
    ! 
    ! Setup file access property list with parallel I/O access.
    !
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
#endif

    !
    ! Create the file collectively.
    ! 
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error,&
         & access_prp = plist_id)
    call h5pclose_f(plist_id, error)

    allocate(h5_slice_data(1:nn(2),1:nn(3)))

    do nf=1,nfields

       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2)
          h5_slice_data(j,k) = slice(nf,j,k)
       end do; end do
       !$omp end parallel do

       call write_hdf5_slice_dataset(trim(datanames(nf)),h5_slice_data,nn)


    end do

    if(RADIATION) then
       call calcqr
       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2)
          h5_slice_data(j,k) = ga(1,j,k)
       end do; end do
       !$omp end parallel do
       
       call write_hdf5_slice_dataset('/G',h5_slice_data, nn)
       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2)
          h5_slice_data(j,k) = sqrt(dot_product(qr(1, j, k, 1:3),&
               qr(1, j, k, 1:3)))
       end do; end do
       !$omp end parallel do
       
       call write_hdf5_slice_dataset('/q',h5_slice_data, nn)
    end if
    
    ! Deallocate data buffer.
    !
    deallocate(h5_slice_data)

    call h5pclose_f(plist_id, error)

    !
    ! Close the file.
    !
    call h5fclose_f(file_id, error)

    if(mpirank==mpiroot) then
       call add_timestamp(filename, time)
       call write_xdmf_slice_file(nn(1),nn(2),nfields,&
            &datanames,filename,nfile)
    end if
    
    !
    ! Close FORTRAN predefined datatypes.
    !
    call h5close_f(error)
        
    return

  end subroutine write_hdf5_slice_file

  subroutine write_xdmf_file(n1,n2,ndatanames,datanames,&
       &data_is_vector,h5filename,nfile)
    use types
    implicit none
    integer(ik), intent(IN)          :: n1, n2
    integer(ik), intent(IN)          :: ndatanames
    character(len=*), dimension(:)   :: datanames
    logical, dimension(:)            :: data_is_vector
    character(len=*)                 :: h5filename
    integer(ik), intent(IN),optional :: nfile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !writes xdmf file corresponding to hdf5 file!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(LEN=256)     :: xdmf_filename ! File name
    integer(ik)            :: xdmf_file
    integer(ik)            :: i

    if(present(nfile)) then
       write(xdmf_filename,'(a,i6.6,a)') 'output.',nfile,'.xmf'
    else
       write(xdmf_filename,'(2a)') trim(h5filename),'.xmf'
    end if

    open(newunit=xdmf_file,file=trim(xdmf_filename),form='formatted')

    write(xdmf_file,'(a)') '<?xml version="1.0" encoding="utf-8"?>'
    write(xdmf_file,'(a)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" &
         &Version="3.0">'
    write(xdmf_file,'(a)') '  <Domain>'
    write(xdmf_file,'(a)') '    <Grid Name="Grid">'
    write(xdmf_file,'(a)') '      <Geometry Origin="" Type="ORIGIN_DXDYDZ">'
    write(xdmf_file,'(a)') '        <DataItem DataType="Float" &
         &Dimensions="3"            Format="XML"'
    write(xdmf_file,'(a)') '		  Precision="8">0 0 0</DataItem>'
    write(xdmf_file,'(a)') '        <DataItem DataType="Float" &
         &Dimensions="3"            Format="XML"'
    write(xdmf_file,'(a)') '		  Precision="8">1 1 1</DataItem>'
    write(xdmf_file,'(a)') '      </Geometry>'
    write(xdmf_file,'(a,3i6,a)') '      <Topology NumberOfElements="',n1,gn2,gn3,'" &
         &Type="3DCoRectMesh"/>'
    do i=1,ndatanames
       if(data_is_vector(i)) then
          write(xdmf_file,'(5a)') '      <Attribute Center="Node" Name="',&
               &trim(datanames(i)),'"             Type="Vector">'
          write(xdmf_file,'(a,i6,a,4i6,a)') '        <DataItem DataType="Float" &
               &Precision="',rks,'" Dimensions="'&
               &,n1,gn2,gn3,3,'"'
          write(xdmf_file,'(5a)') &
               &'            		  Format="HDF">',&
               &trim(h5filename),':/',trim(datanames(i)),'</DataItem>'
          write(xdmf_file,'(a)') '      </Attribute>'
       else
          write(xdmf_file,'(3a)') '      <Attribute Center="Node" Name="',&
               &trim(datanames(i)),'"             Type="Scalar">'
          write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
               &Precision="',rks,'" Dimensions="'&
               &,n1,gn2,gn3, '"'
          write(xdmf_file,'(5a)') &
               &'            		  Format="HDF">',&
               &trim(h5filename),':/',trim(datanames(i)),'</DataItem>'
          write(xdmf_file,'(a)') '      </Attribute>'
       end if
    end do

    write(xdmf_file,'(a)') '    </Grid>'
    write(xdmf_file,'(a)') '  </Domain>'
    write(xdmf_file,'(a)') '</Xdmf>'

    close(xdmf_file,status='keep')

  end subroutine write_xdmf_file

  subroutine write_xdmf_slice_file(n1,n2,ndatanames,datanames,h5filename,nfile)
    use types
    implicit none
    integer(ik), intent(IN)          :: n1, n2
    integer(ik), intent(IN)          :: ndatanames
    character(len=*), dimension(:)   :: datanames
    character(len=*)                 :: h5filename
    integer(ik), intent(IN),optional :: nfile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !writes xdmf file corresponding to hdf5 file!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(LEN=256)     :: xdmf_filename ! File name
    integer(ik)            :: xdmf_file
    integer(ik)            :: i


    write(xdmf_filename,'(a,i6.6,a)') 'slice.',nfile,'.xmf'


    open(newunit=xdmf_file,file=trim(xdmf_filename),form='formatted')

    write(xdmf_file,'(a)') '<?xml version="1.0" encoding="utf-8"?>'
    write(xdmf_file,'(a)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" &
         &Version="3.0">'
    write(xdmf_file,'(a)') '  <Domain>'
    write(xdmf_file,'(a)') '    <Grid Name="Grid">'
    write(xdmf_file,'(a)') '      <Geometry Origin="" Type="ORIGIN_DXDYDZ">'
    write(xdmf_file,'(a)') '        <DataItem DataType="Float" &
         &Dimensions="3"            Format="XML"'
    write(xdmf_file,'(a)') '		  Precision="8">0 0 0</DataItem>'
    write(xdmf_file,'(a)') '        <DataItem DataType="Float" &
         &Dimensions="3"            Format="XML"'
    write(xdmf_file,'(a)') '		  Precision="8">1 1 1</DataItem>'
    write(xdmf_file,'(a)') '      </Geometry>'
    write(xdmf_file,'(a,3i6,a)') '      <Topology NumberOfElements="',1,gn2,gn3,'" &
         &Type="3DRectMesh"/>'
    do i=1,ndatanames

       write(xdmf_file,'(3a)') '      <Attribute Center="Node" Name="',&
            &trim(datanames(i)),'"             Type="Scalar">'
       write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
            &Precision="',rks,'" Dimensions="'&
            &,1,gn2,gn3, '"'
       write(xdmf_file,'(5a)') &
            &'            		  Format="HDF">',&
            &trim(h5filename),':/',trim(datanames(i)),'</DataItem>'
       write(xdmf_file,'(a)') '      </Attribute>'

    end do

    if(RADIATION) then
       write(xdmf_file,'(3a)') '      <Attribute Center="Node" Name="',&
            &'G','"             Type="Scalar">'
       write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
            &Precision="',rks,'" Dimensions="'&
            &,1,gn2,gn3, '"'
       write(xdmf_file,'(5a)') &
            &'            		  Format="HDF">',&
            &trim(h5filename),':/','G','</DataItem>'
       write(xdmf_file,'(a)') '      </Attribute>'

       write(xdmf_file,'(3a)') '      <Attribute Center="Node" Name="',&
            &'q','"             Type="Scalar">'
       write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
            &Precision="',rks,'" Dimensions="'&
            &,1,gn2,gn3, '"'
       write(xdmf_file,'(5a)') &
            &'            		  Format="HDF">',&
            &trim(h5filename),':/','q','</DataItem>'
       write(xdmf_file,'(a)') '      </Attribute>'
    end if

    write(xdmf_file,'(a)') '    </Grid>'
    write(xdmf_file,'(a)') '  </Domain>'
    write(xdmf_file,'(a)') '</Xdmf>'

    close(xdmf_file,status='keep')

  end subroutine write_xdmf_slice_file

  subroutine write_hdf5_pp_file(nn,u,u2,filename,time)
    implicit none
    integer(ik), dimension(1:4), intent(IN)                        :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN) :: u
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(IN) :: u2
    character(len=*), intent(IN)                           :: filename
    real(rk), intent(IN) :: time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Writes HDF5 file in parallel !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=64), dimension(16) :: datanames
    integer(ik)                     :: ndatanames
    real(rks), dimension(:,:,:), allocatable :: h5_scalar_data
    integer(ik) :: i,j,k

    dimsf(1)=nn(1)
    dimsf(2)=gn2
    dimsf(3)=gn3
    dimsf(4)=3

    dimsfi(:)=dimsf(:)

    !
    ! Initialize FORTRAN predefined datatypes
    !
    call h5open_f(error) 
#ifdef _MPI_
    ! 
    ! Setup file access property list with parallel I/O access.
    !
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
#endif
    !
    ! Create the file collectively.
    ! 
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error,&
         & access_prp = plist_id)

    call h5pclose_f(plist_id, error)

    allocate(h5_scalar_data(1:nn(1),1:nn(2),1:nn(3)))

    ndatanames = 0
    if(OUTPUT_W) then
       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          h5_scalar_data(i,j,k)=u(i,j,k,nu1)
       end do; end do ; end do
       !$omp end parallel do
       call write_hdf5_scalar_dataset('/w1',h5_scalar_data,nn)

       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          h5_scalar_data(i,j,k)=u(i,j,k,nu2)
       end do; end do ; end do
       !$omp end parallel do
       call write_hdf5_scalar_dataset('/w2',h5_scalar_data,nn)

       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          h5_scalar_data(i,j,k)=u(i,j,k,nu3)
       end do; end do ; end do
       !$omp end parallel do
       call write_hdf5_scalar_dataset('/w3',h5_scalar_data,nn)

       datanames(1)='w1'
       datanames(2)='w2'
       datanames(3)='w3'
       ndatanames = ndatanames + 3

    end if

    !write energy dissipation rate
    if(OUTPUT_DISS) then
       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          h5_scalar_data(i,j,k)=u2(i,j,k,1)
       end do; end do ; end do
       !$omp end parallel do
       call write_hdf5_scalar_dataset('/diss',h5_scalar_data,nn)
       ndatanames=ndatanames+1
       datanames(ndatanames)='diss'
    end if

    !write passive scalar dissipation
    if(PASSIVE_SCALAR.and.OUTPUT_SCL_DISS) then
       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          h5_scalar_data(i,j,k)=u(i,j,k,nsclf)
       end do; end do ; end do
       !$omp end parallel do
       call write_hdf5_scalar_dataset('/scldiss',h5_scalar_data,nn)
       ndatanames=ndatanames+1
       datanames(ndatanames)='scldiss'
    end if

    !write current
    if(MHD.and.OUTPUT_J) then
       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          h5_scalar_data(i,j,k)=u(i,j,k,nb1)
       end do; end do ; end do
       !$omp end parallel do
       call write_hdf5_scalar_dataset('/j1',h5_scalar_data,nn)

       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          h5_scalar_data(i,j,k)=u(i,j,k,nb2)
       end do; end do ; end do
       !$omp end parallel do
       call write_hdf5_scalar_dataset('/j2',h5_scalar_data,nn)

       !$omp parallel do
       do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
          h5_scalar_data(i,j,k)=u(i,j,k,nb3)
       end do; end do ; end do
       !$omp end parallel do
       call write_hdf5_scalar_dataset('/j3',h5_scalar_data,nn)
       datanames(ndatanames+1)='j1'
       datanames(ndatanames+2)='j2'
       datanames(ndatanames+3)='j3'
       ndatanames=ndatanames+3
    end if

    ! Deallocate data buffer.
    !
    deallocate(h5_scalar_data)


    call h5pclose_f(plist_id, error)

    !
    ! Close the file.
    !
    call h5fclose_f(file_id, error)

    if(mpirank==mpiroot) then
       call add_timestamp(filename, time)
       call write_xdmf_pp_file
    end if

    !
    ! Close FORTRAN predefined datatypes.
    !
    call h5close_f(error)


    return

  contains

    subroutine write_xdmf_pp_file
      use types
      implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !writes xdmf file corresponding to hdf5 file!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      character(LEN=1024) :: xdmf_filename ! File name
      integer(ik)         :: xdmf_file

      write(xdmf_filename,'(2a)') trim(filename),'.xmf'
      open(newunit=xdmf_file,file=trim(xdmf_filename),form='formatted')

      write(xdmf_file,'(a)') '<?xml version="1.0" encoding="utf-8"?>'
      write(xdmf_file,'(a)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" &
           &Version="3.0">'
      write(xdmf_file,'(a)') '  <Domain>'
      write(xdmf_file,'(a)') '    <Grid Name="Grid">'
      write(xdmf_file,'(a)') '      <Geometry Origin="" Type="ORIGIN_DXDYDZ">'
      write(xdmf_file,'(a)') '        <DataItem DataType="Float" &
           &Dimensions="3"            Format="XML"'
      write(xdmf_file,'(a)') '		  Precision="8">0 0 0</DataItem>'
      write(xdmf_file,'(a)') '        <DataItem DataType="Float" &
           &Dimensions="3"            Format="XML"'
      write(xdmf_file,'(a)') '		  Precision="8">1 1 1</DataItem>'
      write(xdmf_file,'(a)') '      </Geometry>'
      write(xdmf_file,'(a,3i6,a)') '      <Topology NumberOfElements="',n1,gn2,gn3,'" &
           &Type="3DCoRectMesh"/>'
      if(OUTPUT_W) then
         write(xdmf_file,'(a)') '      <Attribute Center="Node" Name="w1" &
              &Type="Scalar">'
         write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
              &Precision="',rks,'" Dimensions="',n1,gn2,gn3,'"'
         write(xdmf_file,'(3a)') '		  Format="HDF">',trim(filename),&
              &':/w1</DataItem>'
         write(xdmf_file,'(a)') '      </Attribute>'
         write(xdmf_file,'(a)') '      <Attribute Center="Node" Name="w2" &
              &Type="Scalar">'
         write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
              &Precision="',rks,'" Dimensions="',n1,gn2,gn3,'"'
         write(xdmf_file,'(3a)') '		  Format="HDF">',trim(filename),&
              &':/w2</DataItem>'
         write(xdmf_file,'(a)') '      </Attribute>'
         write(xdmf_file,'(a)') '      <Attribute Center="Node" Name="w3" &
              &Type="Scalar">'
         write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
              &Precision="',rks,'" Dimensions="',n1,gn2,gn3,'"'
         write(xdmf_file,'(3a)') '		  Format="HDF">',trim(filename),&
              &':/w3</DataItem>'
         write(xdmf_file,'(a)') '      </Attribute>'
      end if
      if(OUTPUT_DISS) then
         write(xdmf_file,'(a)') '      <Attribute Center="Node" Name="diss" &
              &Type="Scalar">'
         write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
              &Precision="',rks,'" Dimensions="',n1,gn2,gn3,'"'
         write(xdmf_file,'(3a)') '		  Format="HDF">',&
              &trim(filename), ':/diss</DataItem>'
         write(xdmf_file,'(a)') '      </Attribute>'
      end if
      if(PASSIVE_SCALAR .and. OUTPUT_SCL_DISS) then
         write(xdmf_file,'(a)') '      <Attribute Center="Node" Name="scldiss" &
              &Type="Scalar">'
         write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
              &Precision="',rks,'" Dimensions="',n1,gn2,gn3,'"'
         write(xdmf_file,'(3a)') '		  Format="HDF">',&
              &trim(filename), ':/scldiss</DataItem>'
         write(xdmf_file,'(a)') '      </Attribute>'
      end if
      if(MHD .and. OUTPUT_J) then
         write(xdmf_file,'(a)') '      <Attribute Center="Node" Name="j1" &
              &Type="Scalar">'
         write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
              &Precision="',rks,'" Dimensions="',n1,gn2,gn3,'"'
         write(xdmf_file,'(3a)') '		  Format="HDF">',&
              &trim(filename), ':/j1</DataItem>'
         write(xdmf_file,'(a)') '      </Attribute>'
         write(xdmf_file,'(a)') '      <Attribute Center="Node" Name="j2" &
              &Type="Scalar">'
         write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
              &Precision="',rks,'" Dimensions="',n1,gn2,gn3,'"'
         write(xdmf_file,'(3a)') '		  Format="HDF">',&
              &trim(filename), ':/j2</DataItem>'
         write(xdmf_file,'(a)') '      </Attribute>'
         write(xdmf_file,'(a)') '      <Attribute Center="Node" Name="j3" &
              &Type="Scalar">'
         write(xdmf_file,'(a,i6,a,3i6,a)') '        <DataItem DataType="Float" &
              &Precision="',rks,'" Dimensions="',n1,gn2,gn3,'"'
         write(xdmf_file,'(3a)') '		  Format="HDF">',&
              &trim(filename), ':/j3</DataItem>'
         write(xdmf_file,'(a)') '      </Attribute>'
      end if
      write(xdmf_file,'(a)') '    </Grid>'
      write(xdmf_file,'(a)') '  </Domain>'
      write(xdmf_file,'(a)') '</Xdmf>'

      close(xdmf_file,status='keep')

      return

    end subroutine write_xdmf_pp_file

  end subroutine write_hdf5_pp_file

  subroutine read_hdf5_file(nn,u,filename)
    implicit none
    integer(ik), dimension(1:4), intent(IN)                         :: nn
    real(rks), dimension(1:dim1(nn(1)),1:nn(2),1:nn(3),1:nn(4)), intent(OUT) :: u
    character(*)                                           :: filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Writes HDF5 file in parallel !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(ik) :: i,j,k,l
    integer(hid_t) :: file_id       ! File identifier 
    integer(hid_t) :: dset_id       ! Dataset identifier 
    integer(hid_t) :: filespace     ! Dataspace identifier in file 
    integer(hid_t) :: memspace      ! Dataspace identifier in memory
    integer(hid_t) :: plist_id      ! Property list identifier 


    integer(hsize_t), dimension(1:4) :: dimsf ! Dataset dimensions.
    !     INTEGER, DIMENSION(7) :: dimsfi = (/5,8,0,0,0,0,0/)
    integer(hsize_t), dimension(1:4) :: dimsfi 

    integer(hsize_t), dimension(1:4) :: icount  
    integer(hsize_t), dimension(1:4) :: offset 
    integer :: rank ! Dataset rank 

    integer :: error  ! Error flags

    real(rks), dimension(:,:,:), allocatable :: h5_scalar_data

!!$    rank=4
!!$    dimsf(1)=3
!!$    dimsf(2)=nn(1)
!!$    dimsf(3)=gn2
!!$    dimsf(4)=gn3
!!$    dimsfi(:)=dimsf(:)

    !
    ! Initialize FORTRAN predefined datatypes
    !
    call h5open_f(error)

    call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, gzip_avail, error)

    if(gzip_avail .and. COMPRESSION_LEVEL /= 0) then
       call h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, filter_info, error)
       filter_info_both=ior(H5Z_FILTER_ENCODE_ENABLED_F,&
            &H5Z_FILTER_DECODE_ENABLED_F)
       if(filter_info .ne. filter_info_both) gzip_avail = .false.
    end if
    
#ifdef _MPI_
    ! 
    ! Setup file access property list with parallel I/O access.
    !
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
#endif
    !
    ! Create the file collectively.
    ! 
    call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error,&
         & access_prp = plist_id)

    call h5pclose_f(plist_id, error)

    allocate(h5_scalar_data(1:nn(1),1:nn(2),1:nn(3)))
    allocate(h5_vector_data(1:3,1:nn(1),1:nn(2),1:nn(3)))

    call read_hdf5_vector_dataset('/u',h5_vector_data,nn)
    !$omp parallel do
    do l=1,3 ; do k=1,nn(3) ; do j=1,nn(2) ; do i=1,nn(1)
       u(i,j,k,l)=h5_vector_data(l,i,j,k)
    end do; end do; end do; end do
    !$omp end parallel do

!!$    call read_hdf5_scalar_dataset('/u1',h5_scalar_data,nn)
!!$    u(1:nn(1),1:nn(2),1:nn(3),nu1)=h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))
!!$
!!$    call read_hdf5_scalar_dataset('/u2',h5_scalar_data,nn)
!!$    u(1:nn(1),1:nn(2),1:nn(3),nu2)=h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))
!!$
!!$    call read_hdf5_scalar_dataset('/u3',h5_scalar_data,nn)
!!$    u(1:nn(1),1:nn(2),1:nn(3),nu3)=h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))

!!$    !write passive scalar
!!$    if(PASSIVE_SCALAR) then
!!$       call read_hdf5_scalar_dataset('/scl',h5_scalar_data,nn)
!!$       u(1:nn(1),1:nn(2),1:nn(3),nsclf)=h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))
!!$    end if

    !write mhd
    if(MHD) then
       call read_hdf5_scalar_dataset('/b1',h5_scalar_data,nn)
       u(1:nn(1),1:nn(2),1:nn(3),nb1)=h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))

       call read_hdf5_scalar_dataset('/b2',h5_scalar_data,nn)
       u(1:nn(1),1:nn(2),1:nn(3),nb2)=h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))

       call read_hdf5_scalar_dataset('/b3',h5_scalar_data,nn)
       u(1:nn(1),1:nn(2),1:nn(3),nb3)=h5_scalar_data(1:nn(1),1:nn(2),1:nn(3))
    end if

    ! Deallocate data buffer.
    !
    deallocate(h5_scalar_data)
    deallocate(h5_vector_data)

    call h5pclose_f(plist_id, error)

    !
    ! Close the file.
    !
    call h5fclose_f(file_id, error)

    !
    ! Close FORTRAN predefined datatypes.
    !
    call h5close_f(error)

    return

  contains

    subroutine read_hdf5_scalar_dataset(dataset_name, dataset,nn)
      implicit none
      integer(ik), dimension(1:4), intent(in) :: nn
      character(*), intent(IN) :: dataset_name
      real(rks), dimension(1:nn(1),1:nn(2),1:nn(3)), intent(OUT) :: dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !
      ! Create the data space for the  dataset. 
      !
      !call h5screate_simple_f(rank, dimsf, filespace, error)

      !
      ! Create the dataset with default properties.
      !

      call h5dopen_f(file_id, dataset_name, &
           dset_id, error)

#ifdef _MPI_
      !
      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file. 
      !
      rank = 3
      icount(1) = dimsf(1)
      icount(2) = nn(2)
      icount(3) = nn(3)
      offset(1) = 0
      offset(2) = ljstart-1
      offset(3) = lkstart-1

      call h5screate_simple_f(rank, icount, memspace, error) 
      ! 
      ! Select hyperslab in the file.
      !
      call h5dget_space_f(dset_id, filespace, error)
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
           &icount, error)

      !
      ! Create property list for collective dataset write
      !
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

      !
      ! Write the dataset collectively. 
      !
      if(rks==sp) then
         call h5dread_f(dset_id, H5T_NATIVE_REAL, dataset, dimsfi, error, &
              file_space_id = filespace, mem_space_id = memspace, &
              &xfer_prp = plist_id)
      else
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dimsfi, error, &
              file_space_id = filespace, mem_space_id = memspace, &
              &xfer_prp = plist_id)
      end if
#else
      if(rks==sp) then
         call h5dread_f(dset_id, H5T_NATIVE_REAL, dataset, dimsfi, error)
      else
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dimsfi, error)
      end if
#endif
      !
      ! Write the dataset independently. 
      !
      !    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dimsfi, error, &
      !                     file_space_id = filespace, mem_space_id = memspace)
      !

      !
      ! Close dataspaces.
      !
      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)

      !
      ! Close the dataset
      !
      call h5dclose_f(dset_id, error)

      return

    end subroutine read_hdf5_scalar_dataset

    subroutine read_hdf5_vector_dataset(dataset_name, dataset,nn)
      implicit none
      integer(ik), dimension(1:4), intent(in) :: nn
      character(*), intent(IN) :: dataset_name
      real(rks), dimension(1:3,1:nn(1),1:nn(2),1:nn(3)), intent(OUT) :: dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !
      ! Create the data space for the  dataset. 
      !
      !call h5screate_simple_f(rank, dimsf, filespace, error)

      !
      ! Create the dataset with default properties.
      !

      rank = 4

      icount(1) = 3
      icount(2) = nn(1)
      icount(3) = nn(2)
      icount(4) = nn(3)

      offset(1) = 0
      offset(2) = 0
      offset(3) = ljstart-1
      offset(4) = lkstart-1
      
      
      call h5dopen_f(file_id, dataset_name, &
           dset_id, error)



#ifdef _MPI_
      !
      ! Each process defines dataset in memory and writes it to the hyperslab
      ! in the file. 
      !



      call h5screate_simple_f(rank, icount, memspace, error)


      ! 
      ! Select hyperslab in the file.
      !
      call h5dget_space_f(dset_id, filespace, error)


      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, &
           &icount, error)

      !
      ! Create property list for collective dataset write
      !
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, error)
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

      !
      ! Write the dataset collectively. 
      !
      if(rks==sp) then
         call h5dread_f(dset_id, H5T_NATIVE_REAL, dataset, dimsfi, error, &
              file_space_id = filespace, mem_space_id = memspace, &
              &xfer_prp = plist_id)
      else
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dimsfi, error, &
              file_space_id = filespace, mem_space_id = memspace, &
              &xfer_prp = plist_id)
      end if
#else
      if(rks==sp) then
         call h5dread_f(dset_id, H5T_NATIVE_REAL, dataset, dimsfi, error)
      else
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dataset, dimsfi, error)
      end if
#endif
      !
      ! Write the dataset independently. 
      !
      !    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dimsfi, error, &
      !                     file_space_id = filespace, mem_space_id = memspace)
      !

      !
      ! Close dataspaces.
      !
      call h5sclose_f(filespace, error)
      call h5sclose_f(memspace, error)

      !
      ! Close the dataset
      !
      call h5dclose_f(dset_id, error)

      return

    end subroutine read_hdf5_vector_dataset

  end subroutine read_hdf5_file

end module hdf5_aliakmon
