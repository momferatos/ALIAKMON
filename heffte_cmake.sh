export LDFLAGS='-lrt -lpthread -L /home/giorgos/opt/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin/'
cmake \
    -D CMAKE_C_COMPILER=`which pgcc` \
    -D CMAKE_CXX_COMPILER=`which pgc++`\
    -D MPI_CXX_COMPILER=`which mpic++` \
    -D CMAKE_CXX_FLAGS="-fast -fastsse -tp=native -Mmkl" \
    -D MPI_C_COMPILER=`which mpicc` \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_SHARED_LIBS=ON     \
    -D CMAKE_INSTALL_PREFIX=${HOME}/libs/heffte-NVHPC/ \
    -D Heffte_ENABLE_STOCK=ON \
    -D Heffte_ENABLE_AVX=ON \
    -D Heffte_ENABLE_FFTW=ON \
    -D FFTW_ROOT=${HOME}/libs/fftw-NVHPC/ \
    -D Heffte_ENABLE_CUDA=ON \
    -D CUDA_TOOLKIT_ROOT_DIR=$CUDA_ROOT \
    -D Heffte_ENABLE_FORTRAN=ON \
    -D CUDA_NVCC_FLAGS="-ccbin nvc++ --allow-unsupported-compiler" \
    -D Heffte_DISABLE_GPU_AWARE_MPI=OFF \
    -D Heffte_ENABLE_MKL=ON \
    -D MKL_ROOT=$MKLROOT \
    -D Heffte_MKL_THREAD_LIBS='${HOME}/opt/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin/libiomp5.so' \
    -D Heffte_ENABLE_DOXYGEN=OFF \
    ../heffte
