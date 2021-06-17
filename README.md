!!$     ___ __                                       
!!$ (  / _ \\ \       /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!!$   Copyright (c) 2009-2020 George Momferatos

1) Set environment variables:

# change these lines accordingly
NVCOMPILERS=/home/giorgos/opt/nvidia/hpc_sdk; export NVCOMPILERS
NVARCH=`uname -s`_`uname -m`; export NVARCH
# FFTW3., HDF5, and heFFTE 2.1 are supposed to go here:
export LIBSROOT=/home/giorgos/libs

MPISUBPATH=/openmpi4/openmpi-4.0.5

MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/21.5/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/21.5/compilers/bin:$PATH; export PATH

export PATH=$NVCOMPILERS/$NVARCH/21.5/comm_libs/$MPISUBPATH/bin/:$PATH

export LD_LIBRARY_PATH=$NVCOMPILERS/$NVARCH/21.5/compilers/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NVCOMPILERS/$NVARCH/21.5/comm_libs/$MPISUBPATH/lib/:$LD_LIBRARY_PATH


export LD_LIBRARY_PATH=$LIBSROOT/hdf5/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIBSROOT/heffte/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LIBSROOT/fftw/lib/:$LD_LIBRARY_PATH

export UCX_MEMTYPE_CACHE=n

# Architecture
export ARCH=haswell

2) Build FFTW 3.3.9:

CFLAGS="-fast -fastsse -tp=$ARCH o-Mipa=fast" CXX=pgcpp CC=pgcc F77=pgf77 ./configure --prefix=$LIBS --enable-single --enable-parallel --enable-openmp --enable-fortran --enable-threads --prefix=$LIBSROOT/hdf5
make
make check
make install

3) Build HDF5 1.12.0

CPP=cpp CFLAGS="-fPIC -m64 -tp=$ARCH" CXXFLAGS="-fPIC -m64 -tp=$ARCH" FCFLAGS="-fPIC -m64 -tp=$ARCH" CC=mpicc CXX=mpic++ FC=mpif90 ./configure --enable-threadsafe --enable-fortran --enable-parallel --prefix=$LIBS --enable-unsupported --prefix=$LIBSROOT/hdf5
make
make test
make install

4) Build HeFFTe 2.1 (https://bitbucket.org/icl/heffte/)

mdkir build
cd build
sh ../build.sh

contents of build.sh:

------------------------------
export LDFLAGS='-lrt -lpthread'
cmake \
    -D CMAKE_C_COMPILER=`which pgcc` \
    -D CMAKE_CXX_COMPILER=`which pgc++`\
    -D MPI_CXX_COMPILER=`which mpic++` \
    -D CMAKE_CXX_FLAGS='-gpu=cc75 -gpu=cuda11.3 -acc'\
    -D MPI_C_COMPILER=`which mpicc` \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_SHARED_LIBS=ON     \
    -D CMAKE_INSTALL_PREFIX=$LIBSROOT/heffte/ \
    -D Heffte_ENABLE_FFTW=ON \
    -D FFTW_ROOT=/home/giorgos/libs/fftw/ \
    -D Heffte_ENABLE_CUDA=ON \
    -D CUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-11.3/ \
    -D Heffte_ENABLE_FORTRAN=ON \
    -D CUDA_NVCC_FLAGS='-ccbin nvc++ --gpu-architecture compute_75 --gpu-code compute_75' \
    -D Heffte_DISABLE_GPU_AWARE_MPI=OFF \
    ../heffte
------------------------------

make
make test 


the following tests fail for me, without impact on my code:

3 - heffte_reshape3d_np7 (Failed)
4 - heffte_reshape3d_np12 (Failed)
9 - heffte_fft3d_np8 (Failed)
10 - heffte_fft3d_np12 (Failed)
16 - heffte_fft3d_r2c_np8 (Failed)
17 - heffte_fft3d_r2c_np12 (Failed)

make install

5) build ALIAKMON-GPU with HeFFTe cuFFT backend (CPU + GPU)

cd aliakmon
source config.heffte.cufft
make

6) build ALIAKMON with heFFTe FFTW backend (100% CPU)

cd aliakmon
source config.heffte.fftw
make

7) buld ALIAKMON with Intel MKL (no heFFTe, so decomposition in pencils is
not available)

cd aliakmon
source config.mkl
make -f Makefile.mkl

8) run ALIAKMON-GPU or ALIAKMON

cd aliakmon/
mkdir test_case
cp aliakmon.nml test_case
cd test_case
edit aliakon.nml (see comments inside the file)
mpirun -x OMP_NUM_THREADS=n -x UCX_MEMTYPE_CACHE=n -np n ../aliakmon.{cufft,fftw,mkl}.exe

9) As a check, open the vtk files with paraview

