export LDFLAGS='-lrt -lpthread -L /home/giorgos/opt/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin/'
cmake \
    -D CMAKE_C_COMPILER=`which mpicc` \
    -D CMAKE_CXX_COMPILER=`which mpic++`\
    -D MPI_CXX_COMPILER=`which mpic++` \
    -D CMAKE_CXX_FLAGS='-Ofast -march=native -mtune=native' \
    -D MPI_C_COMPILER=`which mpicc` \
    -D CMAKE_BUILD_TYPE=Release \
    -D BUILD_SHARED_LIBS=ON     \
    -D CMAKE_INSTALL_PREFIX=${HOME}/libs/heffte-gnu/ \
    -D Heffte_ENABLE_STOCK=ON \
    -D Heffte_ENABLE_AVX=OFF \
    -D Heffte_ENABLE_FFTW=OFF \
    -D Heffte_ENABLE_FORTRAN=ON \
    -D Heffte_DISABLE_GPU_AWARE_MPI=OFF \
    -D Heffte_ENABLE_DOXYGEN=OFF \
    ../heffte
