export LDFLAGS='-lrt -lpthread -L /home/giorgos/opt/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin/'
cmake \
    -D CMAKE_C_COMPILER=`which pgcc` \
    -D CMAKE_CXX_COMPILER=`which pgc++`\
    -D MPI_CXX_COMPILER=`which mpic++` \
    -D CMAKE_CXX_FLAGS='-tp=native -gpu=cc75,cuda11.4 -acc -Mmkl'\
    -D MPI_C_COMPILER=`which mpicc` \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_SHARED_LIBS=ON     \
    -D CMAKE_INSTALL_PREFIX=/home/giorgos/libs/heffte/ \
    -D Heffte_ENABLE_STOCK=ON \
    -D Heffte_ENABLE_AVX=ON \
    -D Heffte_ENABLE_FFTW=ON \
    -D FFTW_ROOT=/home/giorgos/libs/fftw/ \
    -D Heffte_ENABLE_CUDA=ON \
    -D CUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-11.4/ \
    -D Heffte_ENABLE_FORTRAN=ON \
    -D CUDA_NVCC_FLAGS='-ccbin nvc++ --gpu-architecture compute_75 --gpu-code compute_75' \
    -D Heffte_DISABLE_GPU_AWARE_MPI=OFF \
    -D Heffte_ENABLE_MKL=ON \
    -D MKL_ROOT=$MKLROOT \
    -D Heffte_MKL_THREAD_LIBS='/home/giorgos/opt/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin/libiomp5.so' \
    -D Heffte_ENABLE_DOXYGEN=ON \
    ../heffte
