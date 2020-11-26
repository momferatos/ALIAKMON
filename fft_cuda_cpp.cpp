#include "heffte.h"
#include <iostream>
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

extern "C"  void fft_cuda_cpp(long n1, long n2, long gn3, long *lksize, long *lkstart, float *rdata, float *cdata, long dir, long alloc);

void fft_cuda_cpp(long n1, long n2, long gn3, long *lksize, long *lkstart, float *rdata, float *cdata, long dir, long alloc) {
  // wrapper around MPI_Comm_rank() and MPI_Comm_size(), using this is optional
     me = heffte::mpi::comm_rank(comm);
     num_ranks = heffte::mpi::comm_size(comm);
    
    // using problem with size 10x20x30 problem
    //heffte::box3d<> real_indexes({0, 0, 0}, {gn3 - 1, n2 - 1, n1 - 1});
    //heffte::box3d<> complex_indexes({0, 0, 0}, {gn3 - 1, n2 - 1, (int)floor(n1 / 2) + 1 - 1});

    heffte::box3d<> real_indexes({0, 0, 0}, {n1-1, n2-1, gn3-1});
    heffte::box3d<> complex_indexes({0, 0, 0}, {(int)floor(n1 / 2) + 1 - 1, n2 - 1, gn3 - 1});
    
    // the dimension where the data will shrink
    int r2c_direction = 0;
    int slice_direction = 2;
    // check if the complex indexes have correct dimension
    assert(real_indexes.r2c(r2c_direction) == complex_indexes);

    /*  
    // create a processor grid with minimum surface (measured in number of indexes)
    std::array<int,3> proc_grid = heffte::proc_setup_min_surface(real_indexes, num_ranks);

    // split all indexes across the processor grid, defines a set of boxes
    std::vector<heffte::box3d<>> real_boxes    = heffte::split_world(real_indexes,    proc_grid);
    std::vector<heffte::box3d<>> complex_boxes = heffte::split_world(complex_indexes, proc_grid);
    */
    
    std::array< int, 3 > const order = {0,1,2};
    std::array<int, 3> proc_grid = heffte::proc_setup_min_surface(real_indexes, num_ranks);

    std::vector<heffte::box3d<>> real_boxes1    = heffte::split_world(real_indexes,    proc_grid);
    std::vector<heffte::box3d<>> complex_boxes1 = heffte::split_world(complex_indexes, proc_grid);

    std::vector<heffte::box3d<>> real_boxes  = heffte::make_slabs(real_indexes,num_ranks,0,1,real_boxes1,order);
    //real_boxes  = heffte::make_slabs(real_indexes,num_ranks,0,1,real_boxes,order);
    std::vector<heffte::box3d<>> complex_boxes  = heffte::make_slabs(complex_indexes,num_ranks,0,1,complex_boxes1,order);
    //complex_boxes  = heffte::make_slabs(complex_indexes,num_ranks,0,1,complex_boxes,order);
    

    // pick the box corresponding to this rank
    heffte::box3d<> const inbox  = real_boxes[me];
    *lkstart = (long)inbox.low[slice_direction];
    *lksize = (long)(inbox.high[slice_direction] - inbox.low[slice_direction] + 1);
    heffte::box3d<> const outbox = complex_boxes[me];
    
    if(alloc == 1) return;

    //if (heffte::gpu::device_count() > 1){
    //    // on a multi-gpu system, distribute the devices across the mpi ranks
    //    heffte::gpu::device_set(heffte::mpi::comm_rank(comm) % heffte::gpu::device_count());
    //}

    // define the heffte class and the input and output geometry
    // heffte::plan_options can be specified just as in the backend::fftw
    heffte::fft3d_r2c<heffte::backend::cufft> fft(inbox, outbox, r2c_direction, comm);
    //float fac = (n2) * (*lksize);
    //printf("%d,%d,%d,%d\n",me, dir, inbox.high[slice_direction], outbox.high[slice_direction]);
    //std::iota(input.begin(), input.end(), 0); // put some data in the input

    // load the input into the GPU memory
    // this is equivalent to cudaMalloc() followed by cudaMemcpy()
    // the destructor of heffte::gpu::vector will call cudaFree()
        
    // perform forward/backward fft using arrays and the user-created workspace
    if(dir==1) {
      // create some input on the CPU
      std::vector<float> input(fft.size_inbox());
            
      //copy input data	
      /*
      int l_c = n1;
      int l_f = 2 * (n1 / 2 + 1);
      int b_c = n2;
      int b_f = n2;
      int h_c = *lksize;
      int h_f = *lksize;
      for(int a=0; a < l_c; a++)
	for(int b=0; b < b_c; b++)
	  for(int c=0; c < h_c; c++) {
	    int i_c =  a*(b_c*h_c) + b*h_c + c;
	    int i_f =  c*(b_c*l_f) + b*l_f + c;
	    input.data()[i_c] = data[i_c];
	      }
      */
#pragma omp parallel for
      for(int i=0; i < fft.size_inbox(); i++)
	input.data()[i] = rdata[i];
      // allocate memory on the device for the output
      heffte::gpu::vector<float> gpu_input = heffte::gpu::transfer::load(input);
      heffte::gpu::vector<std::complex<float>> gpu_output(fft.size_outbox());
      
      // allocate scratch space, this is using the public type alias buffer_container
      // and for the cufft backend this is heffte::gpu::vector
      // for the CPU backends (fftw and mkl) the buffer_container is std::vector
      heffte::fft3d_r2c<heffte::backend::cufft>::buffer_container<std::complex<float>> workspace(fft.size_workspace());
      static_assert(std::is_same<decltype(gpu_output), decltype(workspace)>::value,
                  "the containers for the output and workspace have different types");
      fft.forward(gpu_input.data(), gpu_output.data(), workspace.data(),heffte::scale::none);
     
      // move the result back to the CPU for comparison purposes
      std::vector<std::complex<float>> output = heffte::gpu::transfer::unload(gpu_output);
#pragma omp parallel for
      for(int i=0, j=0; i < fft.size_outbox(); i++,j+=2) {
      	cdata[j] = output.data()[i].real();
      	cdata[j+1] = output.data()[i].imag();
      }
      	    // optional step, free the workspace since the inverse will use the vector API
      //input = std::vector<float>();
      //gpu_input = heffte::gpu::vector<float>();
      //gpu_output = heffte::gpu::vector<std::complex<float>>();
      //workspace = heffte::fft3d_r2c<heffte::backend::cufft>::buffer_container<std::complex<float>>();
      //output = std::vector<std::complex<float>>();
    } else {
	// create some input on the CPU
	std::vector<std::complex<float>> input(fft.size_outbox());
    
	//copy input data
#pragma omp parallel for 
	for(int i=0, j=0; i < fft.size_outbox(); i++,j+=2) {
	  input.data()[i].real(cdata[j]);
	  input.data()[i].imag(cdata[j+1]);
	}
     
        // allocate memory on the device for the output
	heffte::gpu::vector<std::complex<float>> gpu_input = heffte::gpu::transfer::load(input);
	heffte::gpu::vector<float> gpu_output(fft.size_inbox());
	// allocate scratch space, this is using the public type alias buffer_container
	// and for the cufft backend this is heffte::gpu::vector
	// for the CPU backends (fftw and mkl) the buffer_container is std::vector
	heffte::fft3d_r2c<heffte::backend::cufft>::buffer_container<std::complex<float>> workspace(fft.size_workspace());
	static_assert(std::is_same<decltype(gpu_input), decltype(workspace)>::value,
                  "the containers for the output and workspace have different types");
	fft.backward(gpu_input.data(), gpu_output.data(), workspace.data(), heffte::scale::full);
	// move the result back to the CPU for comparison purposes
	std::vector<float> output = heffte::gpu::transfer::unload(gpu_output);
    
	//copy input data
	/*int l_c = n1;
      int l_f = 2 * (n1 / 2 + 1);
      int b_c = n2;
      int b_f = n2;
      int h_c = *lksize;
      int h_f = *lksize;
      for(int a=0; a < l_c; a++)
	for(int b=0; b < b_c; b++)
	  for(int c=0; c < h_c; c++) {
	    int i_c =  a*(b_c*h_c) + b*h_c + c;
	    int i_f =  c*(b_c*l_f) + b*l_f + c;
	    data[i_f] = output.data()[i_c];
	      }
	*/
	
	
#pragma omp parallel for 
	for(int i=0; i < fft.size_inbox(); i++) 
	  rdata[i] = output.data()[i];
	
    }


    // compute the inverse FFT transform using the container API
    //heffte::gpu::vector<float> gpu_inverse = fft.backward(gpu_output);

    
	
    return;

    

      // print the error for each MPI rank
    //std::cout << std::scientific;
    //for(int i=0; i<num_ranks; i++){
    //    MPI_Barrier(comm);
    //    if (me == i) std::cout << "rank " << i << " computed error: " << err << std::endl;
} 

/*
int main(int argc, char** argv){
  long n1, n2, gn3, lksize, lkstart, i;
  float *data, *data2, err;
  n1=n2=gn3=128;
    MPI_Init(&argc, &argv);
    fft_cuda_c(n1, n2, gn3, &lksize, &lkstart, data, 1, true);
    
    data = (float *)calloc((size_t)n1*n2*lksize,sizeof(float));
    data2 = (float *)calloc((size_t)n1*n2*lksize,sizeof(float));
    for(int i=0; i<n1*n2*lksize; i++){
      float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
      data2[i] = data[i] = r;
	}
    fft_cuda_c(n1, n2, gn3, &lksize, &lkstart, data, 1, false);
    fft_cuda_c(n1, n2, gn3, &lksize, &lkstart, data, -1, false);

    // compute the error between the input and the inverse
    err = 0.0;
    for(size_t i=0; i<n1*n2*lksize; i++)
      err = std::max(err, std::abs(data[i] - data2[i]));

    std::cout << std::scientific;
    for(int i=0; i<num_ranks; i++){
       MPI_Barrier(comm);
        if (me == i) std::cout << "rank " << i << " computed error: " << err << std::endl;
    }
    MPI_Finalize();

    return 0;
}
*/
