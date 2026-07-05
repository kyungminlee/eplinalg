# eplinalg examples

Standalone example programs that consume an **installed** eplinalg release via
`find_package`. They are *not* part of the eplinalg build tree — each is its own
CMake project you point at a release you built and installed.

## `mmsolve` — MPI sparse direct solver (MUMPS + Intel MKL)

Reads a sparse matrix and a right-hand side in [Matrix Market][mm] format,
solves `A x = b` with MUMPS, and writes the solution. Precision is selected by a
type prefix, and it runs under MPI.

```
mmsolve -t <s|c|d|z> [-v] <matrix.mtx> <rhs.mtx> <solution.mtx>

  -t   arithmetic:  s = single real   c = single complex
                    d = double real   z = double complex
  -v   leave MUMPS diagnostics on (default: silent)
```

- **matrix** — Matrix Market `coordinate` (sparse), `real` or `complex`,
  `general` or `symmetric`. A `symmetric` matrix stores only its lower triangle
  and is solved with MUMPS `SYM=2` (which is *complex-symmetric* `A = Aᵀ` for
  complex data, not Hermitian). `hermitian` and `skew-symmetric` are rejected
  rather than silently mis-solved — expand them to `general` first.
- **rhs** — Matrix Market `array` (dense), `n × 1`, matching field.
- **solution** — written as a Matrix Market `array` file.

MPI uses MUMPS' centralized-input model: rank 0 reads the files and owns the
assembled matrix, RHS, and solution; every rank participates in the collective
factor/solve. Run with `mpirun -n N ./mmsolve …`.

### Why only s/c/d/z (and why MKL)

MUMPS ships four "genuine" arithmetics — `s c d z` — and those four **coexist in
one binary** (they share the single `mumps_common` runtime and expose distinct
`?mumps_c` entry points). Intel MKL provides a ScaLAPACK/BLACS/PBLAS backend for
exactly those four, so MKL is a drop-in parallel backend for the genuine family.

eplinalg's *extended* precisions (`e y` kind10, `q x` kind16, `m w` multifloat)
each come from a migrated MUMPS stack that re-emits the same `dmumps_*`/`zmumps_*`
symbols as the genuine bridge. Consequently:

- A migrated pair (e.g. `e y`) coexists internally, but **not** with the genuine
  family and **not** with another migrated pair — they collide on those symbols.
- ⇒ **One executable per family, minimum.** A single binary cannot cover all ten
  types, and MKL — which only provides s/c/d/z — cannot serve the extended
  families at all; those require the in-tree migrated ScaLAPACK/LAPACK/BLAS.

This example targets the genuine `s/c/d/z` family on MKL. An extended-family
build from the same sources (linking `eplinalg::<fam>mumps_full` + the in-tree
`eplinalg::<fam>scalapack` instead of MKL) is a straightforward follow-on.

### Building

You need: an installed eplinalg release (its `${LIB_PREFIX}mumps` package, which
bundles the genuine `dzmumps`/`scmumps` solvers), Intel MKL, and MPI.

```sh
cmake -S examples/mmsolve -B build/mmsolve \
      -DCMAKE_PREFIX_PATH="<eplinalg-install>;<mkl-root>/lib/cmake" \
      -DEPLINALG_MUMPS_PACKAGE=emumps \
      -DMKL_MPI=openmpi
cmake --build build/mmsolve
```

Two things **must match the release you installed**:

1. **Compiler family/version.** The release tags its libraries by compiler
   (e.g. `gfortran-15-openmpi-4.1`); `find_package` selects the matching
   pre-built target and errors if your `CMAKE_Fortran_COMPILER` doesn't match.
   Configure with the same compiler the release was built with.
2. **MKL BLACS flavor (`-DMKL_MPI=…`) must match the release's MPI.** An OpenMPI
   release needs `MKL_MPI=openmpi` (`libmkl_blacs_openmpi_lp64`); an Intel-MPI
   release needs `MKL_MPI=intelmpi`. A mismatch fails at link or run time.

`EPLINALG_MUMPS_PACKAGE` names the installed package that carries the genuine
solvers (kind10 → `emumps`, kind16 → `qmumps`, …); the genuine targets inside
are identical regardless of which one you use.

If MKL is installed to a non-standard prefix, add its runtime to the loader path
before running: `export LD_LIBRARY_PATH=<mkl-root>/lib:$LD_LIBRARY_PATH`.

### Running the bundled example

`mmsolve/data/` has a 5×5 real system (solution `[1 2 3 4 5]`) and a 3×3
complex-symmetric system (solution `[1+i, 2−i, −1+2i]`):

```sh
cd examples/mmsolve
mpirun -n 2 ../../build/mmsolve/mmsolve -t d data/A_real.mtx data/b_real.mtx x.mtx
mpirun -n 4 ../../build/mmsolve/mmsolve -t z data/A_cplx.mtx data/b_cplx.mtx xz.mtx
```

### Files

| file | role |
|------|------|
| `mmsolve.c`    | arg parsing, MPI orchestration, per-type MUMPS solve (macro-generated for s/c/d/z) |
| `mmio_min.h/.c`| minimal self-contained Matrix Market reader/writer (no external mmio dependency) |
| `CMakeLists.txt` | standalone consumer project |
| `data/`        | small real + complex test systems with known solutions |

[mm]: https://math.nist.gov/MatrixMarket/formats.html
