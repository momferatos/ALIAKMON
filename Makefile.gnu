MPIFC = mpifort
MPICC = mpicc
MPICXX = mpic++

DEBUGFLAGS = -cpp -g -O0 # -Ktrap=denorm,divz,inexact,inv,ovf,unf
OPTFLAGS= # -fast -fastsse -tp=$(ARCH)
BUILDFLAGS = $(OPTFLAGS)

OMPFLAGS= -D _OPENMP_ -fopenmp
MPIFLAGS= -D _MPI_
PARFLAGS= $(OMPFLAGS) $(MPIFLAGS) $(ACCELFLAGS) 


CXXFLAGS = $(MKLFLAGS) $(BUILDFLAGS) $(PARFLAGS)

FCFLAGS = -cpp -Wunused -Wunused-parameter $(MKLFLAGS) $(BUILDFLAGS) $(PARFLAGS)

FFTWROOT=$(LIBSROOT)/fftw-gnu
HEFFTEROOT=$(LIBSROOT)/heffte-gnu
CUDAROOT=/usr/local/cuda-11.3
HDF5ROOT=$(LIBSROOT)/hdf5-gnu

INCLUDE= -I $(FFTWROOT)/include -I $(HDF5ROOT)/include -I $(HEFFTEROOT)/include -I $(CUDAROOT)/include 

LIB= -L $(FFTWROOT)/lib -L $(HDF5ROOT)/lib -L $(HEFFTEROOT)/lib -L $(CUDAROOT)/lib64 

LDFLAGS =-lpthread -lm -ldl -lhdf5_hl -lhdf5hl_fortran -lhdf5_fortran -lhdf5 -lheffte -lhefftefftwfortran $(CUDALDFLAGS) $(MKLLDFLAGS) -lstdc++ -lmpi_cxx -lfftw3f -lfftw3f_mpi -lfftw3f_threads

OBJS=parameters.o data.o hdf5.o heffte_init.o fft_heffte.o vtk.o numerics.o validation.o initial_conditions.o input_output.o aliakmon.o	

all: aliakmon

aliakmon: $(OBJS)
	$(MPIFC) -o $@.$(BACKEND).exe $^ $(FCFLAGS) $(PRECISION) $(LIB) $(LDFLAGS)

%.o: %.f90
	$(MPIFC) -c -o $@ $< $(FCFLAGS) $(PRECISION) $(INCLUDE) 

%.o: %.cpp
	$(MPICXX) -c -o $@ $< $(CXXFLAGS) $(PRECISION) $(INCLUDE) 
clean:
	rm -f *.o *.mod aliakmon.$(BACKEND).exe 

