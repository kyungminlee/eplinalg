# eplinalg examples

Standalone example programs that consume an **installed** eplinalg release via
`find_package`. They are *not* part of the eplinalg build tree — each is its own
CMake project you point at a release you built and installed.

## `mmsolve` — MPI sparse direct solver (MUMPS + Intel MKL)

Reads a sparse matrix and a right-hand side in [Matrix Market][mm] format,
solves `A x = b` with MUMPS, and writes the solution. Precision is selected by a
type prefix, and it runs under MPI.

```
mmsolve -t <s|c|d|z|e|y|q|x|m|w> [-o <ordering>] [-v] <matrix.mtx> <rhs.mtx> <solution.mtx>

  -t   arithmetic:  s = single real         c = single complex
                    d = double real         z = double complex
                    e = long-double real    y = long-double complex
                    q = __float128 real     x = __float128 complex
                    m = double-double real  w = double-double complex
  -o   ordering: default | pord | scotch | metis | ptscotch
  -v   leave MUMPS diagnostics on (default: silent)
```

The genuine `s c d z` and the extended `e y` (kind10), `q x` (kind16), and
`m w` (multifloats) families **all live in this one binary** — see "How ten
arithmetics share one binary" below.

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

### How ten arithmetics share one binary (and why MKL)

MUMPS ships four "genuine" arithmetics — `s c d z` — that coexist in one binary:
they share the single `mumps_common` runtime and expose distinct `?mumps_c` entry
points. Intel MKL provides a ScaLAPACK/BLACS/PBLAS backend for exactly those four,
so MKL is a drop-in parallel backend for the genuine family.

eplinalg's *extended* precisions add six more `?mumps_c` entry points on the same
shared runtime: `e y` (kind10 long double), `q x` (kind16 `__float128`), and
`m w` (double-double multifloats), each backed by its own in-tree migrated
ScaLAPACK/LAPACK/BLAS (MKL provides none of these precisions). This binary links
**all of them** — the genuine `s c d z` on MKL and the extended families on their
in-tree stacks — so a single `mmsolve` covers all ten. The seam that makes this
safe is symbol layering:

- **Genuine typed** routines (`pdgetrf_`, `pzpotrf_`, …) → **MKL**.
- **Type-agnostic plumbing** (`blacs_gridinit_`, `descinit_`, `numroc_`,
  `Cblacs_gridinfo`) → **MKL** as well. There is exactly *one* copy of each in the
  process, and both the program and MKL's own `pdgetrf` internals use it.
- **Extended typed** routines (`pegetrf_`, `pqpotrf_`, `Cmgamx2d`, …) → the
  in-tree e/y, q/x, and m/w archives, which MKL cannot provide.

This works because MKL's `libmkl_blacs_openmpi_lp64` is the netlib **reference**
BLACS recompiled — it exports the same `BI_*` context internals — so the in-tree
long-double reference routines share MKL's exact BLACS context representation. No
ABI clash, and no `-Wl,--allow-multiple-definition`.

The one non-obvious requirement is **link order**: MKL must appear *ahead* of the
in-tree archives so it wins the genuine + type-agnostic symbols (the in-tree
reference ScaLAPACK also *defines* `pdgetrf_`, so if it came first it would
capture them). The `q x` and `m w` families follow the same MKL-first layering
against their own in-tree archives; their custom MPI reduction operators
additionally require the `quad_mpi_init()` / `multifloats_mpi_init()` calls
`mmsolve` makes right after `MPI_Init`.

### Building

You need: an installed eplinalg release (its `${LIB_PAIR_PREFIX}mumps` packages —
`eymumps`/`qxmumps`/`mwmumps` — which bundle the genuine `dzmumps`/`scmumps`
solvers *and* the extended typed solvers, plus the per-precision
ScaLAPACK/BLAS closures they pull), Intel MKL, and MPI. The consumer
`find_package`s that whole closure — see the comments at the top of
`mmsolve/CMakeLists.txt`.

```sh
cmake -S example/mmsolve -B build/mmsolve \
      -DCMAKE_PREFIX_PATH="<eplinalg-install>;<mkl-root>/lib/cmake" \
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

If MKL is installed to a non-standard prefix, add its runtime to the loader path
before running: `export LD_LIBRARY_PATH=<mkl-root>/lib:$LD_LIBRARY_PATH`.

### Running the bundled example

`mmsolve/data/` has a 5×5 real system (solution `[1 2 3 4 5]`) and a 3×3
complex-symmetric system (solution `[1+i, 2−i, −1+2i]`):

```sh
cd example/mmsolve
mpirun -n 2 ../../build/mmsolve/mmsolve -t d data/A_real.mtx data/b_real.mtx x.mtx
mpirun -n 4 ../../build/mmsolve/mmsolve -t z data/A_cplx.mtx data/b_cplx.mtx xz.mtx
# extended long-double (e/y) run through the same binary:
mpirun -n 2 ../../build/mmsolve/mmsolve -t e data/A_real.mtx data/b_real.mtx xe.mtx
mpirun -n 2 ../../build/mmsolve/mmsolve -t y data/A_cplx.mtx data/b_cplx.mtx xy.mtx
```

### Files

| file | role |
|------|------|
| `mmsolve.c`    | arg parsing, MPI orchestration, per-type MUMPS solve (macro-generated for s/c/d/z/e/y/q/x/m/w) |
| `mmio_min.h/.c`| minimal self-contained Matrix Market reader/writer (no external mmio dependency) |
| `CMakeLists.txt` | standalone consumer project |
| `data/`        | small real + complex test systems with known solutions |

[mm]: https://math.nist.gov/MatrixMarket/formats.html
