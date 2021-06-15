#include "heffte.h"
#include <iostream>
#include "complex.h"
#include "mpi.h"
/*!
 * \brief HeFFTe example 5, using the cuFFT backend.
 *
 * This example is near identical to the first (fftw) example,
 * the main difference is the use of the cufft backend.
 * The interface and types for the cufft backend work the same,
 * with the exception that the array must sit on the GPU device
 * and if the optional vector interface is used then the vector
 * containers are of type heffte::cuda::vector.
 */

//subroutine fft_mpi_alloc(n1, n2, gn3, lksize, lkstart)

int me;
int num_ranks;
MPI_Comm comm = MPI_COMM_WORLD;

extern "C" void fft_alloc_cpp(long n1, long n2, long gn3, int slice_direction, int r2c_direction, int *il1, int *il2, 
				    int *il3, int *ih1, int *ih2, int *ih3,int *ol1, int *ol2, int *ol3, int *oh1, int *oh2, int *oh3);


void fft_alloc_cpp(long n1, long n2, long gn3, int slice_direction, int r2c_direction, int *il1, int *il2, 
			int *il3, int *ih1, int *ih2, int *ih3,int *ol1, int *ol2, int *ol3, int *oh1, int *oh2, int *oh3) {
  //me = heffte::mpi::comm_rank(comm);
  //num_ranks = heffte::mpi::comm_size(comm);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    // using problem with size 10x20x30 problem
    //heffte::box3d<> real_indexes({0, 0, 0}, {gn3 - 1, n2 - 1, n1 - 1});
    //heffte::box3d<> complex_indexes({0, 0, 0}, {gn3 - 1, n2 - 1, (int)floor(n1 / 2) + 1 - 1});

    std::array<int,3> a = {0,0,0};
     std::array<int,3> b = {n1-1, n2-1, gn3-1};
     std::array<int,3> c = {0,0,0};
     std::array<int,3> d = {floor(n1 / 2) + 1 - 1, n2 - 1, gn3 - 1};
      heffte::box3d<> real_indexes(a,b);
    heffte::box3d<> complex_indexes(c,d);
    
    // the dimension where the data will shrink
    //int r2c_direction = 0;
    //int slice_direction = 2;
    // check if the complex indexes have correct dimension
    assert(real_indexes.r2c(r2c_direction) == complex_indexes);

    /*  
    // create a processor grid with minimum surface (measured in number of indexes)
    std::array<int,3> proc_grid = heffte::proc_setup_min_surface(real_indexes, num_ranks);

    // split all indexes across the processor grid, defines a set of boxes
    std::vector<heffte::box3d<>> real_boxes    = heffte::split_world(real_indexes,    proc_grid);
    std::vector<heffte::box3d<>> complex_boxes = heffte::split_world(complex_indexes, proc_grid);
    */
    
    std::array< int, 3 > const order = {0, 1, 2};
    std::array<int, 3> proc_grid = heffte::proc_setup_min_surface(real_indexes, num_ranks);

    std::vector<heffte::box3d<>> real_boxes1    = heffte::split_world(real_indexes,    proc_grid);
    std::vector<heffte::box3d<>> complex_boxes1 = heffte::split_world(complex_indexes, proc_grid);

    std::vector<heffte::box3d<>> real_boxes  = heffte::make_slabs(real_indexes,num_ranks,0,1,real_boxes1,order);
    //real_boxes  = heffte::make_slabs(real_indexes,num_ranks,0,1,real_boxes,order);
    std::vector<heffte::box3d<>> complex_boxes  = heffte::make_slabs(complex_indexes,num_ranks,0,1,complex_boxes1,order);
    //complex_boxes  = heffte::make_slabs(complex_indexes,num_ranks,0,1,complex_boxes,order);
    

    // pick the box corresponding to this rank
    heffte::box3d<> const inbox  = real_boxes[me];
    //*lkstart = (long)inbox.low[slice_direction];
    //*lksize = (long)(inbox.high[slice_direction] - inbox.low[slice_direction] + 1);
    heffte::box3d<> const outbox = complex_boxes[me];

    *il1 = (int)inbox.low[0] + 1;
    *il2 = (int)inbox.low[1] + 1;
    *il3 = (int)inbox.low[2] + 1;
    *ih1 = (int)inbox.high[0] + 1;
    *ih2 = (int)inbox.high[1] + 1;
    *ih3 = (int)inbox.high[2] + 1;

    *ol1 = (int)outbox.low[0] + 1;
    *ol2 = (int)outbox.low[1] + 1;
    *ol3 = (int)outbox.low[2] + 1;
    *oh1 = (int)outbox.high[0] + 1;
    *oh2 = (int)outbox.high[1] + 1;
    *oh3 = (int)outbox.high[2] + 1;

#ifdef _CUDA_
        if (heffte::gpu::device_count() > 1){
        // on a multi-gpu system, distribute the devices across the mpi ranks
        heffte::gpu::device_set(heffte::mpi::comm_rank(comm) % heffte::gpu::device_count());
    }
#endif
    return;
    
}

