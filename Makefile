MPIFC = mpifort
MPICC = mpicc
MPICXX = mpic++

ARCH=haswell

DEBUGFLAGS = -cpp -g -O0 -Ktrap=denorm,divz,inexact,inv,ovf,unf
OPTFLAGS= -fast -fastsse -tp=$(ARCH)
BUILDFAGS = $(DEBUGFLAGS)
CXXFLAGS = $(OPTFLAGS)

FFLAGS= -cpp
PREC = # -D _DOUBLE_
OMPFLAGS= -D _OPENMP_ -mp=multicore
MPIFLAGS= -D _MPI_
PARFLAGS= $(BUILDFLAGS) $(OMPFLAGS) $(MPIFLAGS) $(ACCELFLAGS) 


FFTWROOT=${LIBSROOT}/fftw
HEFFTEROOT=${LIBSROOT}/heffte
CUDAROOT=/usr/local/cuda-11.3
HDF5ROOT=${LIBSROOT}/hdf5

INCLUDE= -I $(FFTWROOT)/include -I $(HDF5ROOT)/include -I $(HEFFTEROOT)/include -I $(CUDAROOT)/include 

LIB= -L $(FFTWROOT)/lib -L $(HDF5ROOT)/lib -L $(HEFFTEROOT)/lib -L $(CUDAROOT)/lib64 

LDFLAGS =-lpthread -lm -ldl -lhdf5_hl -lhdf5hl_fortran -lhdf5_fortran -lhdf5 -lheffte -lhefftefftwfortran $(CUDALDFLAGS) -lstdc++ -lmpi_cxx -lfftw3f -lfftw3f_mpi -lfftw3f_threads

all: aliakmon

aliakmon: $(OBJS)
	$(MPIFC) -o $@.$(BACKEND).exe $^ $(FFLAGS) $(PARFLAGS) $(LIB) $(LDFLAGS)

%.o: %.f90
	$(MPIFC) -c -o $@ $< $(FFLAGS) $(PREC) $(PARFLAGS) $(INCLUDE) 

%.o: %.cpp
	$(MPICXX) -c -o $@ $< $(CXXFLAGS) $(PREC) $(PARFLAGS) $(INCLUDE) 

clean:
	rm -f *.o *.mod aliakmon.exe 

