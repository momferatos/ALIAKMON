MPIFC = mpifort
MPICC = mpicc
MPICXX = mpic++

DEBUGFLAGS = -cpp -g -O0 -Mbounds # -Ktrap=denorm,divz,inexact,inv,ovf,unf
OPTFLAGS=-fast -fastsse # -tp=zen2
BUILDFLAGS = $(OPTFLAGS)

OMPFLAGS= -D _OPENMP_ -mp=multicore
MPIFLAGS= -D _MPI_
PARFLAGS= $(MPIFLAGS) $(ACCELFLAGS) $(OMPFLAGS)


CXXFLAGS = -Mmkl $(MKLFLAGS) $(FFTWFLAGS) $(BUILDFLAGS) $(PARFLAGS)

FCFLAGS = -cpp $(MKLFLAGS) $(FFTWFLAGS) $(BUILDFLAGS) $(PARFLAGS)

FFTWROOT=$(LIBSROOT)/fftw
HEFFTEROOT=$(LIBSROOT)/heffte
CUDAROOT=$(CUDA_ROOT)
HDF5ROOT=$(LIBSROOT)/hdf5
INCLUDE= -I $(FFTWROOT)/include -I $(HDF5ROOT)/include -I $(HEFFTEROOT)/include  -I $(CUDAROOT)/include 

LIB= -L $(FFTWROOT)/lib -L $(HDF5ROOT)/lib -L $(HEFFTEROOT)/lib  -L $(CUDAROOT)/lib64 -L $(SZIPROOT)/lib

LDFLAGS =-lpthread -lm -ldl -lhdf5_hl -lhdf5hl_fortran -lhdf5_fortran -lhdf5 -lheffte -lhefftefftwfortran -lhefftestockfortran $(CUDALDFLAGS) $(MKLLDFLAGS) -lstdc++ -lmpi_cxx -lfftw3f -lfftw3f_mpi -lfftw3f_threads -lz 

OBJS=parameters.o data.o heffte_init.o fft_heffte.o numerics.o fvdom.o hdf5.o validation.o initial_conditions.o input_output.o  aliakmon.o	

all: aliakmon

aliakmon: $(OBJS)
	$(MPIFC) -o $@.$(BACKEND).exe $^ $(FCFLAGS) $(PRECISION) $(LIB) $(LDFLAGS)

%.o: %.f90
	$(MPIFC) -c -o $@ $< $(FCFLAGS) $(PRECISION) $(INCLUDE) 

%.o: %.cpp
	$(MPICXX) -c -o $@ $< $(CXXFLAGS) $(PRECISION) $(INCLUDE) 
clean:
	rm -f *.o *.mod aliakmon.$(BACKEND).exe 

