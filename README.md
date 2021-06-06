!!$     ___ __                                       
!!$ (  / _ \\ \       /                               
!!$   | |_| |\ \  _  __  ___  ___   _  __   __  _  __
!!$   |  _  | > \| |/  \/ / |/ / | | |/ / _ \ \| |/ /
!!$   | | | |/ ^ \ ( ()  <|   <| |_| | |_/ \_| | / / 
!!$   |_| |_/_/ \_\_)__/\_\_|\_\ ._,_|\___^___/|__/  
!!$                            |_|
!!$  
!!$   Copyright (c) 2009-2020 George Momferatos
1) NVHPC environment variables:

NVARCH=`uname -s`_`uname -m`; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/21.5/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/21.5/compilers/bin:$PATH; export PATH
#export PATH=$NVCOMPILERS/$NVARCH/21.5/comm_libs/openmpi4/openmpi-4.0.5/bin/:$PATH
export PATH=$NVCOMPILERS/$NVARCH/21.5/comm_libs/mpi/bin/:$PATH

export LD_LIBRARY_PATH=$NVCOMPILERS/$NVARCH/21.5/compilers/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NVCOMPILERS/$NVARCH/21.5/comm_libs/mpi/lib/:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=/home/giorgos/libs/hdf5/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/giorgos/libs/heffte/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/giorgos/libs/fftw/lib/:$LD_LIBRARY_PATH

export UCX_MEMTYPE_CACHE=n


3) Libraries go here:

export LIBS=/home/giorgos/libs/

2) Build FFTW 3.3.9:

CFLAGS="-fast -fastsse -tp sandybridge -Mipa=fast" CXX=pgcpp CC=pgcc F77=pgf77 ./configure --prefix=$LIBS --enable-single --enable-parallel --enable-openmp --enable-fortran --enable-threads
make
make check
make install

3) Build HDF5 1.12.0

CPP=cpp CFLAGS="-fPIC -m64 -tp=sandybridge" CXXFLAGS="-fPIC -m64 -tp=sandybridge" FCFLAGS="-fPIC -m64 -tp=sandybridge" CC=mpicc CXX=mpic++ FC=mpif90 ./configure --enable-threadsafe --enable-fortran --enable-parallel --prefix=$LIBS --enable-unsupported
make
make test
make install

4) Build HeFFTe (https://bitbucket.org/icl/heffte/)

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
    -D CMAKE_INSTALL_PREFIX=/home/giorgos/libs/heffte/ \
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

5) build ALIAKMON-GPU with HeFFTe CUDA backend

cd aliakmon
make -f Makefile.pgi.aris

6) build ALIAKMON with heFFTe FFTW backend

cd aliakmon
comment out $(CUDAFFLAGS) in Makefile.heffte
make -f Makefile.heffte

7) Build ALIAKMON with Intel MKL FFT (this gets the most out of the CPUs on ARIS)

cd aliakmon
make -f Makefile.intel.aris

8) run ALIAKMON-GPU or ALIAKMON

cd aliakmon/
mkdir test_case
cp aliakmon.nml test_case
cd test_case
edit aliakon.nml (see comments inside this file)
mpirun -x OMP_NUM_THREADS=n -x UCX_MEMTYPE_CACHE=n -np n ../aliakmon.exe

9) As a check, open the vtk files with paraview

