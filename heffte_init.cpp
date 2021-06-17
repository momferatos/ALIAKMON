#include "heffte.h"
#include <iostream>
#include "complex.h"
#include <math.h>

int me;
int num_ranks;
MPI_Comm comm = MPI_COMM_WORLD;

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


void heffte_init_slabs(long n1, long n2, long n3,
			 int slice_direction, int r2c_direction,
			 int *il1, int *il2, int *il3,
			 int *ih1, int *ih2, int *ih3,
			 int *ol1, int *ol2, int *ol3,
			 int *oh1, int *oh2, int *oh3)  {
  me = heffte::mpi::comm_rank(comm);
  num_ranks = heffte::mpi::comm_size(comm);
  
  std::array<int,3> a = {0,0,0};
  std::array<int,3> b = {n1-1, n2-1, n3-1};
  std::array<int,3> c = {0,0,0};
  std::array<int,3> d = {floor(n1 / 2) + 1 - 1, n2 - 1, n3 - 1};
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
  
  std::vector<heffte::box3d<>> real_boxes  =
    heffte::make_slabs(real_indexes,num_ranks,0,1,real_boxes1,order);
  
  std::vector<heffte::box3d<>> complex_boxes  =
    heffte::make_slabs(complex_indexes,num_ranks,0,1,complex_boxes1,order);
  
  
  
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
        if (heffte::gpu::device_count() > 1){
	  // on a multi-gpu system, distribute the devices across the mpi ranks
	  heffte::gpu::device_set(heffte::mpi::comm_rank(comm)
				  % heffte::gpu::device_count());
	}
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
  std::array<int,3> b = {n1-1, n2-1, n3-1};
  std::array<int,3> c = {0,0,0};
  std::array<int,3> d = {floor(n1 / 2) + 1 - 1, n2 - 1, n3 - 1};
  heffte::box3d<> real_indexes(a,b);
  heffte::box3d<> complex_indexes(c,d);
    
    
  // check if the complex indexes have correct dimension
  assert(real_indexes.r2c(r2c_direction) == complex_indexes);

  
    
  std::array< int, 3 > const order = {0, 1, 2};
  std::array<int, 3> proc_grid_tmp;
  std::array<int, 2> proc_grid;
  int sqr_mpi_size = int(floor(sqrt(float(num_ranks))));
  int pencils_y = sqr_mpi_size;
  int pencils_z = num_ranks / sqr_mpi_size;
  if(pencils_y * pencils_z != num_ranks) pencils_y += 1;
  assert(pencils_y * pencils_z == num_ranks);
  proc_grid[0] = pencils_y;
  proc_grid[1] = pencils_z;
  proc_grid_tmp[0] = 1;
  proc_grid_tmp[1] = pencils_y;
  proc_grid_tmp[2] = pencils_z;

  printf("%d: %d x %d\n", me, pencils_y, pencils_z);
  std::vector<heffte::box3d<>> real_boxes1    =
    heffte::split_world(real_indexes,    proc_grid_tmp);
  std::vector<heffte::box3d<>> complex_boxes1 =
    heffte::split_world(complex_indexes, proc_grid_tmp);
  
  std::vector<heffte::box3d<>> real_boxes  =
    heffte::make_pencils(real_indexes,proc_grid,pencil_direction,real_boxes1,order);
  
  std::vector<heffte::box3d<>> complex_boxes  =
    heffte::make_pencils(complex_indexes,proc_grid,pencil_direction,complex_boxes1,order);
  
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
        if (heffte::gpu::device_count() > 1){
	  // on a multi-gpu system, distribute the devices across the mpi ranks
	  heffte::gpu::device_set(heffte::mpi::comm_rank(comm)
				  % heffte::gpu::device_count());
	}
#endif
    return;
    
}

