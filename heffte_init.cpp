#include "heffte.h"
#include <iostream>
#include "complex.h"
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "mpi.h"
#ifdef _CUFFT_
#include "cuda.h"
#endif

int me;
int num_ranks;
MPI_Comm comm = MPI_COMM_WORLD;

extern "C" int numdevice, mpisize, mpirank;

extern "C" void heffte_init_slabs(long n1, long n2, long n3,
			      int slice_direction, int r2c_direction,
			      int *il1, int *il2, int *il3,
			      int *ih1, int *ih2, int *ih3,
			      int *ol1, int *ol2, int *ol3,
			      int *oh1, int *oh2, int *oh3);

extern "C" void heffte_init_pencils(long n1, long n2, long n3,
			      int pencil_direction, int r2c_direction,
			      int *il1, int *il2, int *il3,
			      int *ih1, int *ih2, int *ih3,
			      int *ol1, int *ol2, int *ol3,
			      int *oh1, int *oh2, int *oh3);

extern "C" void heffte_set_num_device(int numd);

void num_of_pencils(int mpi_size, int *pencils_y, int *pencils_z);
 

void heffte_set_num_device(int numd) {
#ifdef _CUDA_
  heffte::gpu::device_set(numd);
#endif 
  return;
}

void heffte_init_slabs(long n1, long n2, long n3,
			 int slice_direction, int r2c_direction,
			 int *il1, int *il2, int *il3,
			 int *ih1, int *ih2, int *ih3,
			 int *ol1, int *ol2, int *ol3,
			 int *oh1, int *oh2, int *oh3)  {
  me = mpirank;
  num_ranks = mpisize;
  //MPI_Comm_rank(MPI_COMM_WORLD, &me);
  //MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
  //me = heffte::mpi::comm_rank(comm);
  //num_ranks = heffte::mpi::comm_size(comm);
  
  std::array<int,3> a = {0,0,0};
  std::array<int,3> b = {(int)(n1-1), (int)(n2-1), (int)(n3-1)};
  std::array<int,3> c = {0,0,0};
  std::array<int,3> d = {(int)(floor(n1 / 2)), (int)(n2 - 1), (int)(n3 - 1)};
  heffte::box3d<> real_indexes(a,b);
  heffte::box3d<> complex_indexes(c,d);
    
    
  // check if the complex indexes have correct dimension
  assert(real_indexes.r2c(r2c_direction) == complex_indexes);

  
    
  std::array< int, 3 > const order = {0, 1, 2};
  std::array<int, 3> proc_grid = heffte::proc_setup_min_surface(real_indexes,
								num_ranks);
  
  std::vector<heffte::box3d<>> real_boxes1    =
    heffte::split_world(real_indexes,    proc_grid);
  std::vector<heffte::box3d<>> complex_boxes1 =
    heffte::split_world(complex_indexes, proc_grid);

  heffte::rank_remap rank_remap;
  
  std::vector<heffte::box3d<>> real_boxes  =
    heffte::make_slabs(real_indexes,num_ranks,0,1,real_boxes1,order,rank_remap);
  
  std::vector<heffte::box3d<>> complex_boxes  =
    heffte::make_slabs(complex_indexes,num_ranks,0,1,complex_boxes1,order,
		       rank_remap);
    
  // pick the box corresponding to this rank
  heffte::box3d<> const inbox  = real_boxes[me];
  
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
  
#ifdef _CUFFT_
  //if (heffte::gpu::device_count() > 1){
	  // on a multi-gpu system, distribute the devices across the mpi ranks
  //  	  heffte::gpu::device_set(numdevice);
	  //	}
#endif
    return;
    
}

void heffte_init_pencils(long n1, long n2, long n3,
			 int pencil_direction, int r2c_direction,
			 int *il1, int *il2, int *il3,
			 int *ih1, int *ih2, int *ih3,
			 int *ol1, int *ol2, int *ol3,
			 int *oh1, int *oh2, int *oh3) {
  me = heffte::mpi::comm_rank(comm);
  num_ranks = heffte::mpi::comm_size(comm);
  
  std::array<int,3> a = {0,0,0};
  std::array<int,3> b = {(int)(n1-1), (int)(n2-1), (int)(n3-1)};
  std::array<int,3> c = {0,0,0};
  std::array<int,3> d = {(int)floor((int)(n1 / 2)), (int)(n2 - 1), (int)(n3 - 1)};
  heffte::box3d<> real_indexes(a,b);
  heffte::box3d<> complex_indexes(c,d);
    
    
  // check if the complex indexes have correct dimension
  assert(real_indexes.r2c(r2c_direction) == complex_indexes);

  
    
  std::array< int, 3 > const order = {0, 1, 2};
  std::array<int, 3> proc_grid_tmp;
  std::array<int, 2> proc_grid;
  
  int pencils_y;
  int pencils_z;
  num_of_pencils(num_ranks, &pencils_y, &pencils_z);
  
  proc_grid[0] = pencils_y;
  proc_grid[1] = pencils_z;
  proc_grid_tmp[0] = 1;
  proc_grid_tmp[1] = pencils_y;
  proc_grid_tmp[2] = pencils_z;


  std::vector<heffte::box3d<>> real_boxes1    =
    heffte::split_world(real_indexes,    proc_grid_tmp);
  std::vector<heffte::box3d<>> complex_boxes1 =
    heffte::split_world(complex_indexes, proc_grid_tmp);

  heffte::rank_remap rank_remap;
    
  std::vector<heffte::box3d<>> real_boxes  =
    heffte::make_pencils(real_indexes,proc_grid,pencil_direction,
			 real_boxes1,order,rank_remap);
  
  std::vector<heffte::box3d<>> complex_boxes  =
    heffte::make_pencils(complex_indexes,proc_grid,pencil_direction,
			 complex_boxes1,order,rank_remap);
  
  // pick the box corresponding to this rank
  heffte::box3d<> const inbox  = real_boxes[me];
  
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
  
#ifdef _CUFFT_
  //        if (heffte::gpu::device_count() > 1){
	  // on a multi-gpu system, distribute the devices across the mpi ranks
  //  	  heffte::gpu::device_set(numdevice);
	  //	}
#endif
    return;
    
}

void num_of_pencils(int mpi_size, int *pencils_y, int *pencils_z) {
  int n, m, nfacs, sqrt_mpi_size;
  int *facs;
  
  nfacs=0;
  n = 2;
  m = mpi_size;
  while(m != 1)  {
    if(m % n == 0) {
      nfacs++;
      m /= n;
    } else
      n++;
  }
  
  facs = (int *)calloc(nfacs, sizeof(int));

  nfacs=0;
  n = 2;
  m = mpi_size;
  while(m != 1)  {
    if(m % n == 0) {
      facs[nfacs++] = n;
      m /= n;
    } else
      n++;
  }

  sqrt_mpi_size = int(floor(sqrt(float(num_ranks))));
  
  *pencils_y = 1;
  for(n = 0; n < nfacs; n++) {
    *pencils_y *= facs[n];
    if(*pencils_y >= sqrt_mpi_size) break;
  }

  *pencils_z = mpi_size / *pencils_y;

  assert(*pencils_y * *pencils_z == mpi_size);

  if(*pencils_y != 1 && *pencils_z != 1) 
    if(me == 0) printf("Using pencil FFT decomposition: %d x %d.\n", *pencils_y,
	   *pencils_z);  
  
  return;
}

