# ============================================================================
# ALIAKMON — portable CMake build (drop-in replacement for the Makefile)
#
# Usage:
#   cmake -B build -DALIAKMON_BACKEND=fftw -DHeffte_ROOT=/path/to/heffte
#   cmake --build build -j
#
# IMPORTANT: the MPI Fortran module must match the Fortran compiler (a gfortran
# build cannot read an nvfortran-built mpi.mod and vice versa). The simplest way
# is to point CMake at the MPI compiler wrappers for your toolchain, e.g.:
#   cmake -B build -DCMAKE_Fortran_COMPILER=mpifort \
#         -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ...
#
# Common options (defaults in []):
#   -DALIAKMON_BACKEND=stock|fftw|mkl|cufft   heFFTe FFT backend     [fftw]
#   -DALIAKMON_DOUBLE=ON|OFF                   double precision       [OFF]
#   -DALIAKMON_OPENMP=ON|OFF                   OpenMP threading       [OFF]
#   -DALIAKMON_MPI=ON|OFF                      distributed (MPI)      [ON]
#   -DCMAKE_BUILD_TYPE=Release|Debug|...                             [Release]
#
# Dependency hints (if not in a standard location), any of:
#   -DHeffte_ROOT=...   -DHDF5_ROOT=...   -DMKL_ROOT=...
#   or the matching environment variables.
#
# RPATH: the built executable records the full paths of the libraries it links
# (heFFTe, HDF5, ...), so it runs without setting LD_LIBRARY_PATH.
# ============================================================================
cmake_minimum_required(VERSION 3.18)
project(aliakmon VERSION 1.0 LANGUAGES Fortran CXX)

# ----------------------------------------------------------------------------
# Options
# ----------------------------------------------------------------------------
set(ALIAKMON_BACKEND "fftw" CACHE STRING "heFFTe FFT backend: stock, fftw, mkl, cufft")
set_property(CACHE ALIAKMON_BACKEND PROPERTY STRINGS stock fftw mkl cufft)
option(ALIAKMON_DOUBLE "Build in double precision" OFF)
option(ALIAKMON_OPENMP "Enable OpenMP threading"   OFF)
option(ALIAKMON_MPI    "Enable distributed MPI"    ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Make .f90 sources run through the C preprocessor (they use #ifdef _MPI_ etc.).
# Portable across gfortran (-cpp) and other compilers.
set(CMAKE_Fortran_PREPROCESS ON)

# Use RPATH from the build tree so linked shared libs are found at runtime.
set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# ----------------------------------------------------------------------------
# Sources  (CMake scans `module`/`use` to order Fortran compilation itself —
# no manual ordering needed, unlike the Makefile's OBJS list).
# ----------------------------------------------------------------------------
set(ALIAKMON_FSOURCES
  parameters.f90
  data.f90
  fft_heffte.f90
  numerics.f90
  hdf5.f90
  validation.f90
  initial_conditions.f90
  input_output.f90
  aliakmon.f90)

set(ALIAKMON_CXXSOURCES heffte_init.cpp)
if(ALIAKMON_BACKEND STREQUAL "cufft")
  list(APPEND ALIAKMON_CXXSOURCES fft_cuda.cpp)
endif()

add_executable(aliakmon ${ALIAKMON_FSOURCES} ${ALIAKMON_CXXSOURCES})
set_target_properties(aliakmon PROPERTIES
  OUTPUT_NAME "aliakmon.${ALIAKMON_BACKEND}"
  SUFFIX ".exe")

# Use DT_RPATH (old dtags) instead of DT_RUNPATH so the executable's RPATH also
# resolves the TRANSITIVE dependencies of the libraries it links (e.g. heFFTe's
# own dependency on MKL) — RUNPATH only covers direct dependencies.
if(NOT APPLE AND NOT WIN32)
  target_link_options(aliakmon PRIVATE "-Wl,--disable-new-dtags")
endif()

# ----------------------------------------------------------------------------
# Preprocessor defines (mirror the config.<backend> files)
# ----------------------------------------------------------------------------
target_compile_definitions(aliakmon PRIVATE
  $<$<BOOL:${ALIAKMON_MPI}>:_MPI_>
  $<$<BOOL:${ALIAKMON_DOUBLE}>:_DOUBLE_>
  $<$<BOOL:${ALIAKMON_OPENMP}>:_OPENMP_>)

if(ALIAKMON_BACKEND STREQUAL "mkl")
  target_compile_definitions(aliakmon PRIVATE _MKL_)
elseif(ALIAKMON_BACKEND STREQUAL "fftw")
  target_compile_definitions(aliakmon PRIVATE _FFTW_)
elseif(ALIAKMON_BACKEND STREQUAL "cufft")
  target_compile_definitions(aliakmon PRIVATE _CUFFT_)
endif()  # stock => no backend define (the #else branch in fft_heffte.f90)

# ----------------------------------------------------------------------------
# MPI
# ----------------------------------------------------------------------------
if(ALIAKMON_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran CXX)
  target_link_libraries(aliakmon PRIVATE MPI::MPI_Fortran MPI::MPI_CXX)
endif()

# ----------------------------------------------------------------------------
# HDF5 (C + Fortran bindings; Fortran for hdf5.mod / hdf5_fortran).
# The code does collective MPI-IO, so prefer the parallel build when MPI is on
# (override with -DHDF5_ROOT=<parallel hdf5 prefix> if the wrong one is picked).
# ----------------------------------------------------------------------------
if(ALIAKMON_MPI)
  set(HDF5_PREFER_PARALLEL ON)
endif()
find_package(HDF5 REQUIRED COMPONENTS Fortran C HL)
if(ALIAKMON_MPI AND NOT HDF5_IS_PARALLEL)
  message(WARNING "Found a SERIAL HDF5 but MPI is enabled; parallel output "
                  "needs a parallel HDF5. Set -DHDF5_ROOT=<parallel hdf5 prefix>.")
endif()
target_include_directories(aliakmon PRIVATE ${HDF5_Fortran_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
target_link_libraries(aliakmon PRIVATE
  ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES}
  ${HDF5_C_LIBRARIES} ${HDF5_HL_LIBRARIES})

# ----------------------------------------------------------------------------
# heFFTe: main C++ library + the backend-specific Fortran interface lib.
# The Fortran modules (heffte_<backend>.mod) live in heFFTe's include dir.
# Prefer the installed CMake package; fall back to manual find_*.
# ----------------------------------------------------------------------------
find_package(Heffte CONFIG QUIET HINTS ${Heffte_ROOT} $ENV{Heffte_ROOT} $ENV{HEFFTE_ROOT})

# Locate the include dir holding the Fortran .mod files (and heffte.h).
find_path(HEFFTE_MOD_DIR
  NAMES heffte_stock.mod heffte.h
  HINTS ${Heffte_ROOT} $ENV{Heffte_ROOT} $ENV{HEFFTE_ROOT}
  PATH_SUFFIXES include)

# Core library + the always-needed stock Fortran interface.
find_library(HEFFTE_LIB        NAMES heffte
  HINTS ${Heffte_ROOT} $ENV{Heffte_ROOT} $ENV{HEFFTE_ROOT} PATH_SUFFIXES lib lib64)
find_library(HEFFTE_F_STOCK    NAMES hefftestockfortran
  HINTS ${Heffte_ROOT} $ENV{Heffte_ROOT} $ENV{HEFFTE_ROOT} PATH_SUFFIXES lib lib64)

if(NOT HEFFTE_MOD_DIR OR NOT HEFFTE_LIB OR NOT HEFFTE_F_STOCK)
  message(FATAL_ERROR
    "heFFTe not found. Pass -DHeffte_ROOT=<install prefix> "
    "(expected <prefix>/include/heffte_stock.mod and <prefix>/lib/libheffte*).")
endif()

# Fortran interface libs to link: the base stock lib is always needed; non-stock
# backends add their own. (Don't pre-`set()` HEFFTE_F_BACKEND — a normal var
# would shadow find_library's cache result and read back empty.)
set(ALIAKMON_HEFFTE_FLIBS ${HEFFTE_F_STOCK})
if(NOT ALIAKMON_BACKEND STREQUAL "stock")
  find_library(HEFFTE_F_BACKEND NAMES heffte${ALIAKMON_BACKEND}fortran
    HINTS ${Heffte_ROOT} $ENV{Heffte_ROOT} $ENV{HEFFTE_ROOT} PATH_SUFFIXES lib lib64)
  if(NOT HEFFTE_F_BACKEND)
    message(FATAL_ERROR
      "heFFTe Fortran interface 'heffte${ALIAKMON_BACKEND}fortran' not found. "
      "Was heFFTe built with -DHeffte_ENABLE_FORTRAN=ON and this backend enabled? "
      "Pass -DHeffte_ROOT=<install prefix>.")
  endif()
  list(PREPEND ALIAKMON_HEFFTE_FLIBS ${HEFFTE_F_BACKEND})
endif()

target_include_directories(aliakmon PRIVATE ${HEFFTE_MOD_DIR})
target_link_libraries(aliakmon PRIVATE ${ALIAKMON_HEFFTE_FLIBS} ${HEFFTE_LIB})
if(TARGET Heffte::Heffte)
  target_link_libraries(aliakmon PRIVATE Heffte::Heffte)
endif()

# ----------------------------------------------------------------------------
# Backend-specific extras
# ----------------------------------------------------------------------------
if(ALIAKMON_BACKEND STREQUAL "mkl")
  # heFFTe's MKL backend needs the Intel MKL runtime. We deliberately do NOT use
  # find_package(MKL): its ILP64 default adds -fdefault-integer-8, which makes
  # integer literals 8-byte and breaks resolution of heFFTe's integer(C_INT)
  # constructors ("Too many components in structure constructor").
  #
  # heFFTe's own libs declare the MKL dependency (libmkl_intel_ilp64) but use
  # DT_RUNPATH (their own dir), so the exe's RPATH cannot resolve their MKL
  # child. We therefore force MKL to be a DIRECT dependency of the exe
  # (--no-as-needed, since our Fortran code references no MKL symbols) so the
  # exe's own RPATH loads it; once loaded it satisfies heFFTe globally.
  find_path(MKL_RUNTIME_LIBDIR
    NAMES libmkl_intel_ilp64.so libmkl_core.so libmkl_rt.so
    HINTS ${MKL_ROOT} $ENV{MKLROOT} $ENV{MKL_ROOT} ${MKLROOT}
    PATH_SUFFIXES lib lib/intel64 lib64 intel64)
  if(MKL_RUNTIME_LIBDIR)
    find_library(MKL_ILP64_LIB NAMES mkl_intel_ilp64 HINTS ${MKL_RUNTIME_LIBDIR})
    target_link_libraries(aliakmon PRIVATE
      "-Wl,--no-as-needed" "${MKL_ILP64_LIB}" "-Wl,--as-needed")
    set_property(TARGET aliakmon APPEND PROPERTY BUILD_RPATH   "${MKL_RUNTIME_LIBDIR}")
    set_property(TARGET aliakmon APPEND PROPERTY INSTALL_RPATH "${MKL_RUNTIME_LIBDIR}")
    message(STATUS "ALIAKMON: MKL runtime RPATH -> ${MKL_RUNTIME_LIBDIR}")
  else()
    message(WARNING "MKL runtime lib dir not found; set MKLROOT or run with "
                    "LD_LIBRARY_PATH=<mkl>/lib. (Build still succeeds.)")
  endif()
elseif(ALIAKMON_BACKEND STREQUAL "cufft")
  # The cuFFT backend offloads the fields to the GPU via OpenACC (the !$acc
  # directives in numerics.f90 etc.) and hands device pointers to heFFTe's
  # cuFFT backend. It therefore REQUIRES the NVIDIA HPC Fortran compiler
  # (nvfortran) — gfortran ignores !$acc and the fields would stay on the host.
  # Configure with -DCMAKE_Fortran_COMPILER=mpifort (wrapping nvfortran) or
  # =nvfortran, and use the matching NVHPC MPI.
  if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    message(WARNING
      "cufft backend needs nvfortran (NVHPC) for OpenACC GPU offload, but the "
      "Fortran compiler is '${CMAKE_Fortran_COMPILER_ID}'. The build may succeed "
      "while silently NOT offloading to the GPU. Set "
      "-DCMAKE_Fortran_COMPILER=<nvfortran or its mpi wrapper>.")
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  # OpenACC GPU offload flags (NVHPC). -cuda enables CUDA interop so device
  # pointers can be passed to heFFTe; -Minfo=accel reports offload decisions.
  # Override the accel target with e.g. -DALIAKMON_ACC_FLAGS="-acc=gpu;-gpu=cc75".
  set(ALIAKMON_ACC_FLAGS "-acc=gpu;-cuda;-Minfo=accel"
      CACHE STRING "NVHPC OpenACC/GPU flags for the cuFFT backend")
  target_compile_options(aliakmon PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${ALIAKMON_ACC_FLAGS}>)
  target_link_options(aliakmon PRIVATE -acc=gpu -cuda)
  target_link_libraries(aliakmon PRIVATE CUDA::cufft CUDA::cudart CUDA::cuda_driver)
endif()

# ----------------------------------------------------------------------------
# OpenMP
# ----------------------------------------------------------------------------
if(ALIAKMON_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran CXX)
  target_link_libraries(aliakmon PRIVATE OpenMP::OpenMP_Fortran OpenMP::OpenMP_CXX)
endif()

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
message(STATUS "ALIAKMON: backend=${ALIAKMON_BACKEND} "
               "double=${ALIAKMON_DOUBLE} mpi=${ALIAKMON_MPI} omp=${ALIAKMON_OPENMP} "
               "-> aliakmon.${ALIAKMON_BACKEND}.exe")
