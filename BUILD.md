# Building ALIAKMON (CMake)

A portable CMake build (`CMakeLists.txt`) that replaces the hand-written
Makefile. It uses `find_package`/`find_library` for all dependencies (no
hardcoded paths), lets CMake order the Fortran modules automatically, and bakes
the dependency directories into the executable's RPATH so it runs without
`LD_LIBRARY_PATH`.

## Quick start

```bash
# gfortran + system OpenMPI, FFTW backend:
cmake -B build \
  -DCMAKE_Fortran_COMPILER=mpifort \
  -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
  -DALIAKMON_BACKEND=fftw \
  -DHeffte_ROOT=/path/to/heffte
cmake --build build -j
# -> build/aliakmon.fftw.exe
```

## Options

| Option | Values | Default | Meaning |
|--------|--------|---------|---------|
| `ALIAKMON_BACKEND` | `stock` `fftw` `mkl` `cufft` | `fftw` | heFFTe FFT backend (sets `_MKL_`/`_FFTW_`/`_CUFFT_`; `stock` = none) |
| `ALIAKMON_DOUBLE`  | `ON`/`OFF` | `OFF` | double precision (`-D _DOUBLE_`) |
| `ALIAKMON_OPENMP`  | `ON`/`OFF` | `OFF` | OpenMP threading (`-D _OPENMP_`) |
| `ALIAKMON_MPI`     | `ON`/`OFF` | `ON`  | distributed MPI (`-D _MPI_`) |
| `CMAKE_BUILD_TYPE` | `Release`/`Debug`/… | `Release` | optimisation level |

Dependency hints (cache var or environment): `Heffte_ROOT`, `HDF5_ROOT`,
`MKLROOT`, `CUDAToolkit_ROOT`.

## Backends

- **stock** — heFFTe's built-in reference FFT. Slow, but needs only
  `libheffte` + `libhefftestockfortran`. Good for a smoke test; **builds and
  runs end-to-end here.**
- **fftw** — needs `libhefftefftwfortran` and FFTW. Recommended for single-node.
- **mkl** — needs `libhefftemklfortran` and the Intel MKL runtime. See gotcha #2.
- **cufft** — GPU backend. Needs CUDA + `libhefftecufftfortran`, and
  **requires the NVHPC `nvfortran` compiler**: the fields
  are moved to the GPU with OpenACC (`!$acc`) directives and device pointers are
  passed to heFFTe's cuFFT backend, so it is built with `-acc=gpu -cuda`
  (overridable via `-DALIAKMON_ACC_FLAGS`). gfortran cannot do this offload.
  Configure with the NVHPC MPI wrappers, e.g.
  `-DCMAKE_Fortran_COMPILER=mpifort` (wrapping nvfortran), and ensure the GPU's
  compute capability is among those heFFTe/CUDA were built for (set
  `-DALIAKMON_ACC_FLAGS="-acc=gpu;-gpu=cc75"` for a GTX 1650, for example).

## RPATH (no `LD_LIBRARY_PATH` needed)

The executable records the directories of heFFTe, HDF5, and (for `mkl`) the MKL
runtime in its RPATH, using **DT_RPATH** (`--disable-new-dtags`) so it also
covers transitive dependencies. For the MKL backend, MKL is force-linked as a
**direct** dependency (`--no-as-needed`) because heFFTe's own libs use
DT_RUNPATH and would otherwise hide their MKL child from the exe's RPATH.

## Gotchas (environment, not the build)

1. **MPI must match the Fortran compiler.** A gfortran build cannot read an
   nvfortran-built `mpi.mod` (and vice versa) — you get
   *"… mpi.mod … is not a GNU Fortran module file."* `find_package(MPI)` may
   pick the wrong MPI if several are installed. Fix: point CMake at the matching
   MPI wrappers, e.g. `-DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc
   -DCMAKE_CXX_COMPILER=mpicxx`.

2. **The MKL runtime must match the one heFFTe was built against.** If heFFTe's
   MKL backend was built against an older MKL, running against a much newer MKL
   (e.g. oneAPI 2026.0) can abort at FFT init (SIGABRT) even though everything
   links and loads. Fix: rebuild heFFTe against the current MKL, or point
   `MKLROOT` at the MKL heFFTe was built with.

   - heFFTe's MKL Fortran interface uses `integer(C_INT)`. Do **not** build with
     MKL's ILP64 interface — it adds `-fdefault-integer-8`, which makes integer
     literals 8-byte and breaks the `heffte_fft3d_r2c_mkl(...)` constructor
     (*"Too many components in structure constructor"*). The CMakeLists avoids
     this by not using `find_package(MKL)` for the link.

## Status (validated here)

| Backend | Configure | Build | Run |
|---------|-----------|-------|-----|
| stock   | ✅ | ✅ | ✅ (KE/ε correct) |
| mkl     | ✅ | ✅ | loads w/o `LD_LIBRARY_PATH`; aborts at FFT init (gotcha #2) |

Backends `fftw`/`cufft` were not built here (Fortran interface libs / CUDA not
present in this environment) but are wired in the CMakeLists.
