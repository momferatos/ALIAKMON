MPIFC = mpifort
MPICC = mpicc
MPICXX = mpic++

ARCH=haswell

GFLAGS = -cpp -g -O0 #-Ktrap=denorm,divz,inexact,inv,ovf,unf
OPTFLAGS=-fast -fastsse -tp=$(ARCH)
CXXFLAGS = $(OPTFLAGS)

FFLAGS= -cpp
PREC = # -D _DOUBLE_
OMPFLAGS= -D _OPENMP_ -mp=multicore
MPIFLAGS= -D _MPI_
CUDAFLAGS = -D _CUDA_ -acc=gpu -gpu=cc75,cuda11.3 -cuda
ROCMFLAGS =
ACCELFLAGS  = $(CUDAFLAGS)
PARFLAGS= $(OPTFLAGS) $(OMPFLAGS) $(MPIFLAGS) $(ACCELFLAGS) 

LIBS = /home/giorgos/libs
FFTWROOT=$(LIBS)/fftw
HEFFTEROOT=$(LIBS)/heffte
CUDAROOT=/usr/local/cuda-11.3
HDF5ROOT=$(LIBS)/hdf5

INCLUDE= -I $(FFTWROOT)/include -I $(HDF5ROOT)/include -I $(HEFFTEROOT)/include -I $(CUDAROOT)/include 

LIB= -L $(FFTWROOT)/lib -L $(HDF5ROOT)/lib -L $(HEFFTEROOT)/lib -L $(CUDAROOT)/lib64 

LDFLAGS =-lpthread -lm -ldl -lhdf5_hl -lhdf5hl_fortran -lhdf5_fortran -lhdf5 -lheffte -lhefftefftwfortran -lhefftecufftfortran -lcufft -lcudart -lstdc++ -lmpi_cxx -lfftw3f -lfftw3f_mpi -lfftw3f_threads

FFTWOBJS = parameters.o data.o hdf5.o fft_alloc_cpp.o fft_fftw_f90.o fft_omp.o vtk.o numerics.o validation.o initial_conditions.o input_output.o aliakmon.o

CUDAOBJS = parameters.o data.o hdf5.o fft_alloc_cpp.o fft_cuda_f90.o fft_omp.o vtk.o numerics.o validation.o initial_conditions.o input_output.o aliakmon.o

all: aliakmon.fftw

aliakmon.cuda: $(CUDAOBJS)
	$(MPIFC) -o $@.exe $^ $(FFLAGS) $(PARFLAGS) $(LIB) $(LDFLAGS)
	
aliakmon.fftw: $(FFTWOBJS)
	$(MPIFC) -o $@.exe $^ $(FFLAGS) $(PARFLAGS) $(LIB) $(LDFLAGS) 

%.o: %.f90
	$(MPIFC) -c -o $@ $< $(FFLAGS) $(PREC) $(PARFLAGS) $(INCLUDE) 

%.o: %.cpp
	$(MPICXX) -c -o $@ $< $(CXXFLAGS) $(PREC) $(PARFLAGS) $(INCLUDE) 

clean:
	rm -f *.o *.mod aliakmon.exe 

