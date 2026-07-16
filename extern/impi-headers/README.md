# Intel MPI headers (vendored, compile-time only)

These are the C and Fortran MPI headers from Intel MPI Library, vendored
so the MUMPS build can compile against a real MPI ABI without requiring
every contributor to install the runtime locally.

## Provenance

Extracted from `pip install impi-devel==2021.17.2`:

- `mpi.h`  ← `<venv>/include/mpi.h`   (C API, 285 KB)
- `mpif.h` ← `<venv>/include/mpif.h`  (Fortran F77 include, 23 KB)

Only these two headers are vendored. MUMPS 5.8.2 uses `INCLUDE 'mpif.h'`
in 231 files and never `USE mpi`, so Fortran `.mod` files are not needed
(and would be useless anyway — Intel `.mod` files are ifx/ifort-specific
and not readable by gfortran).

## What's NOT here

- **Runtime libraries** (`libmpi.so`, `libmpifort.so`, …) — not vendored.
  At link/run time, use the `impi-rt` Python package (installed as a
  transitive dep of `impi-devel`) or any system MPI (OpenMPI, MPICH).
  The ABI at the header level is stable across Intel MPI releases.
- **Compiler wrappers** (`mpicc`, `mpifort`) — not needed; the build
  passes `-I` and `-L` flags directly.
- **Fortran `.mod` files** — compiler-specific binaries, would not work
  with gfortran. Irrelevant since MUMPS doesn't `USE mpi`.
- **MPI I/O headers** (`mpio.h`, `mpiof.h`) — MUMPS doesn't call any
  `MPI_File_*` / `MPI_IO` routines.

## License

Intel Simplified Software License (October 2022). See `LICENSE.txt` in
this directory. The license permits redistribution without modification
provided the copyright notice and license text are preserved alongside.

Intel(R) MPI Library: Copyright (C) 2009 Intel Corporation.

## Refreshing

To pull a newer version:

```
python3 -m venv /tmp/impi && /tmp/impi/bin/pip install impi-devel
cp /tmp/impi/lib/python*/site-packages/../../include/mpi.h external/impi-headers/
cp /tmp/impi/lib/python*/site-packages/../../include/mpif.h external/impi-headers/
cp /tmp/impi/lib/python*/site-packages/../../share/doc/mpi/licensing/license.txt \
   external/impi-headers/LICENSE.txt
# update version in this README
```
