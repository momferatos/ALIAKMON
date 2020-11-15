MPIFC=mpifort

FFTWROOT=/home/giorgos
HDF5ROOT=/home/giorgos

GOPT = -cpp -g -traceback -O0 -Ktrap=denorm,divz,inexact,inv,ovf,unf #-O0 -g -fp-model=strict -fp-model=except -check bounds -stand=f08 -convert big_endian
FOPT= -cpp -fast -fastsse -tp=host -m64 #-ftz -O3 -xCORE-AVX-I -ip -ipo -diag-file=aliakmon.diag -stand=f08 -convert big_endian
PREC = #-D_DOUBLE_
OMPOPT= -D_OPENMP_ -mp
MPIOPT= -D_MPI_ 
PAROPT= $(OMPOPT) $(MPIOPT)
GPAROPT = $(OMPOPT) $(MPIOPT)
INCLUDE=-I $(FFTWROOT)/include -I $(HDF5ROOT)/include/ #-I $(MKLROOT)/include/fftw/ 
LIB=-L $(HDF5ROOT)/lib/ -L $(HOME)/lib # -L $(MKLROOT)/lib/intel64_lin/ 
LIBS = -lfftw3_mpi -lfftw3_omp -lfftw3f_mpi -lfftw3f_omp -lfftw3f -lfftw3 -lhdf5_hl -lhdf5hl_fortran -lhdf5_fortran -lhdf5 -lsz -lz # -Wl,-Bstatic -Wl,--start-group $(HOME)/lib/libfftw3x_cdft_lp64.a $(HOME)/lib/libfftw3xf_intel.a ${MKLROOT}/lib/intel64/libmkl_cdft_core.a ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -Wl,-Bdynamic -liomp5 -lpthread -lm -ldl -lhdf5_hl -lhdf5hl_fortran -lhdf5_fortran -lhdf5 -lsz -lz

OBJS = parameters.o data.o hdf5.o fft_mpi.o fft_omp.o vtk.o numerics.o validation.o initial_conditions.o input_output.o aliakmon.o

all: aliakmon

%.o: %.f90
	$(MPIFC) -c -o $@ $< $(FOPT) $(PREC) $(PAROPT) $(INCLUDE) 

aliakmon: $(OBJS)
	$(MPIFC) -o $@.exe $^ $(FOPT) $(PAROPT) $(LIB) $(LIBS) 

debug: $(OBJS)
	$(MPIFC) -o $@ $^ $(GOPT) $(GPAROPT) $(LIB) $(LIBS) 

clean:
	rm -f *.o *.mod aliakmon.exe 

