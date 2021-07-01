#include "heffte.h"
#include <iostream>
#include "complex.h"
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "cuda.h"
#include "mpi.h"
#include "hwloc.h"

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

extern "C" void print_gpus();
  
void num_of_pencils(int mpi_size, int *pencils_y, int *pencils_z);
 
#define ABORT_ON_ERROR(func)                          \
  { CUresult res;                                     \
    res = func;                                       \
    if (CUDA_SUCCESS != res) {                        \
        printf("%s returned error=%d\n", #func, res); \
        abort();                                      \
    }                                                 \
  }
static hwloc_topology_t topology = NULL;
static int gpuIndex = 0;
static hwloc_obj_t gpus[16] = {0};
 
/**
 * This function searches for all the GPUs that are hanging off a NUMA
 * node.  It walks through each of the PCI devices and looks for ones
 * with the NVIDIA vendor ID.  It then stores them into an array.
 * Note that there can be more than one GPU on the NUMA node.
 */
 
static void find_gpus(hwloc_topology_t topology, hwloc_obj_t parent, hwloc_obj_t child) {
    hwloc_obj_t pcidev;
    pcidev = hwloc_get_next_child(topology, parent, child);
    if (NULL == pcidev) {
        return;
    } else if (0 != pcidev->arity) {
        /* This device has children so need to look recursively at them */
        find_gpus(topology, pcidev, NULL);
        find_gpus(topology, parent, pcidev);
    } else {
        if (pcidev->attr->pcidev.vendor_id == 0x10de) {
            gpus[gpuIndex++] = pcidev;
        }
        find_gpus(topology, parent, pcidev);
    }
}

void print_gpus()
{
    int rank, retval, length;
    char procname[MPI_MAX_PROCESSOR_NAME+1];
    const unsigned long flags = HWLOC_TOPOLOGY_FLAG_IO_DEVICES | HWLOC_TOPOLOGY_FLAG_IO_BRIDGES;
    hwloc_cpuset_t newset;
    hwloc_obj_t node, bridge;
    char pciBusId[16];
    CUdevice dev;
    char devName[256];
 
    
    rank = me;
    if (MPI_SUCCESS != MPI_Get_processor_name(procname, &length)) {
        strcpy(procname, "unknown");
    }
 
    /* Now decide which GPU to pick.  This requires hwloc to work properly.
     * We first see which CPU we are bound to, then try and find a GPU nearby.
     */
    retval = hwloc_topology_init(&topology);
    assert(retval == 0);
    retval = hwloc_topology_set_flags(topology, flags);
    assert(retval == 0);
    retval = hwloc_topology_load(topology);
    assert(retval == 0);
    newset = hwloc_bitmap_alloc();
    retval = hwloc_get_last_cpu_location(topology, newset, 0);
    assert(retval == 0);
 
    /* Get the object that contains the cpuset */
    node = hwloc_get_first_largest_obj_inside_cpuset(topology, newset);
 
    /* Climb up from that object until we find the HWLOC_OBJ_NODE */
    while (node->type != HWLOC_OBJ_NODE) {
        node = node->parent;
    }
 
    /* Now look for the HWLOC_OBJ_BRIDGE.  All PCI busses hanging off the
     * node will have one of these */
    bridge = hwloc_get_next_child(topology, node, NULL);
    while (bridge->type != HWLOC_OBJ_BRIDGE) {
        bridge = hwloc_get_next_child(topology, node, bridge);
    }
 
    /* Now find all the GPUs on this NUMA node and put them into an array */
    find_gpus(topology, bridge, NULL);
 
    ABORT_ON_ERROR(cuInit(0));
    /* Now select the first GPU that we find */
    if (gpus[0] == 0) {
      printf("%s:, No GPU found\n", procname);
        exit(1);
    } else {
        sprintf(pciBusId, "%.2x:%.2x:%.2x.%x", gpus[0]->attr->pcidev.domain, gpus[0]->attr->pcidev.bus,
        gpus[0]->attr->pcidev.dev, gpus[0]->attr->pcidev.func);
        ABORT_ON_ERROR(cuDeviceGetByPCIBusId(&dev, pciBusId));
        ABORT_ON_ERROR(cuDeviceGetName(devName, 256, dev));
        printf("rank=%d (%s): Selected GPU=%s, name=%s\n", rank, procname, pciBusId, devName);
    }
 

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
  
  std::vector<heffte::box3d<>> real_boxes  =
    heffte::make_pencils(real_indexes,proc_grid,pencil_direction,
			 real_boxes1,order);
  
  std::vector<heffte::box3d<>> complex_boxes  =
    heffte::make_pencils(complex_indexes,proc_grid,pencil_direction,
			 complex_boxes1,order);
  
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

