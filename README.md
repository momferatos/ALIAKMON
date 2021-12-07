!!$     ___ __                                       
!!$ (  / _ \\ \       /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!!$   Copyright (c) 2009-2022 George Momferatos

In the following, $ALIAKMON_ROOT is the directory containing ALIAKMON's
source code

For building ALIAKMON, the NVHPC SDK is recommended:

https://developer.nvidia.com/hpc-sdk

#---------------------------------------------------------------
# Set environment variables - edit these lines accordingly
#---------------------------------------------------------------

NVVERSION=21.9

NVCOMPILERS=$HOME/opt/nvidia/hpc_sdk; export NVCOMPILERS
NVARCH=`uname -s`_`uname -m`; export NVARCH

# FFTW3., HDF5, and heFFTE are supposed to go here:
export LIBSROOT=$HOME/libs


MPISUBPATH=/openmpi4/openmpi-4.0.5

MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/$NVVERSION/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/$NVVERSION/compilers/bin:$PATH; export PATH

export PATH=$NVCOMPILERS/$NVARCH/$NVVERSION/comm_libs/$MPISUBPATH/bin/:$PATH

export LD_LIBRARY_PATH=$NVCOMPILERS/$NVARCH/$NVVERSION/compilers/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NVCOMPILERS/$NVARCH/$NVVERSION/comm_libs/$MPISUBPATH/lib/:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=$LIBSROOT/hdf5-NVHPC/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIBSROOT/heffte-NVHPC/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIBSROOT/fftw-NVHPC/lib/:$LD_LIBRARY_PATH

export UCX_MEMTYPE_CACHE=n

#---------------------------------------------------------------
# Build FFTW
#---------------------------------------------------------------

chdir into FFTw source directory
sh $ALIAKMON_ROOT/configure_fftw.sh
make
make check
make install

#---------------------------------------------------------------
# Build HDF5
#---------------------------------------------------------------

chdir into HDF5 source directory
sh $ALIAKMON_ROOT/configure_hdf5.sh
make
make check
make install

#---------------------------------------------------------------
# Build heFFTe
#---------------------------------------------------------------

git clone https://bitbucket.org/icl/heffte
cd heffte
mdkir build
cd build
(First edit the file $ALIAKMON_ROOT/heffte_cmake.sh)
sh $ALIAKMON_ROOT/heffte_cmake.sh
make
make test (Some tests fail simply because of the unavailability of enough MPI slots)
make install

#---------------------------------------------------------------
# build ALIAKMON-GPU with the HeFFTe cuFFT backend (CPU + GPU)
#---------------------------------------------------------------

cd $ALIAKMON_ROOT
source config.cufft
make clean 
make all

#---------------------------------------------------------------
# build ALIAKMON with the heFFTe FFTW backend (CPU)
#---------------------------------------------------------------

cd aliakmon
source config.fftw
make clean 
make all

#---------------------------------------------------------------
# build ALIAKMON with the heFFTe MKL backend (CPU)
#---------------------------------------------------------------

cd aliakmon
source config.mkl
make clean 
make all

#---------------------------------------------------------------
# build ALIAKMON with the heFFTe STOCK backend
#---------------------------------------------------------------

cd aliakmon
source config.stock
make clean 
make all

#---------------------------------------------------------------
# run ALIAKMON
#---------------------------------------------------------------

cd $ALIAKMON_ROOT
mkdir test_case
cp aliakmon.nml test_case
cd test_case
edit aliakon.nml (see comments inside the file)
mpirun -x OMP_NUM_THREADS=nthreads -x UCX_MEMTYPE_CACHE=n -np nprocs ../aliakmon.{cufft,fftw,mkl,stock}.exe

#---------------------------------------------------------------
# check results
#---------------------------------------------------------------

python $ALIAKMON_ROOT/python/plot.py $ALIAKMON_ROOT/test_case/slice.??????.py {'e','w', etc.} (see .py file for details)

or

open the .vtk or .xmf files with paraview. 

