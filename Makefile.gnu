MPIFC = mpifort
MPICC = mpicc
MPICXX = mpic++

DEBUGFLAGS = -fpp -g 
OPTFLAGS=
BUILDFLAGS = -std=f2018 -g -Wunused

OMPFLAGS= -D _OPENMP_ -fopenmp
MPIFLAGS= -D _MPI_
PARFLAGS= $(MPIFLAGS) $(OMPFLAGS)


CXXFLAGS = -std=c++11 $(MKLFLAGS) $(FFTWFLAGS) $(BUILDFLAGS) $(PARFLAGS)

FCFLAGS = -cpp $(MKLFLAGS) $(FFTWFLAGS) $(BUILDFLAGS) $(PARFLAGS)


HEFFTEROOT=$(LIBSROOT)/heffte-gnu
INCLUDE= -I $(HDF5ROOT)/include -I $(HEFFTEROOT)/include
HDF5ROOT=$(HDF5_DIR)

LIB= -L $(MKLROOT)/lib -L $(HDF5ROOT)/lib -L $(HEFFTEROOT)/lib -L $(SZIPROOT)/lib

LDFLAGS =-lpthread -lm -ldl -lhdf5_hl -lhdf5hl_fortran -lhdf5_fortran -lhdf5 -lheffte -lhefftestockfortran -lstdc++ -lz -lsz

OBJS=parameters.o data.o heffte_init.o fft_heffte.o numerics.o hdf5.o validation.o initial_conditions.o input_output.o  aliakmon.o	

all: aliakmon

aliakmon: $(OBJS)
	$(MPIFC) -o $@.$(BACKEND).exe $^ $(FCFLAGS) $(PRECISION) $(LIB) $(LDFLAGS)

%.o: %.f90
	$(MPIFC) -c -o $@ $< $(FCFLAGS) $(PRECISION) $(INCLUDE) 

%.o: %.cpp
	$(MPICXX) -c -o $@ $< $(CXXFLAGS) $(PRECISION) $(INCLUDE) 
clean:
	rm -f *.o *.mod aliakmon.$(BACKEND).exe 

