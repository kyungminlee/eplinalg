# Hand-written parallel BLAS overlay — design & phased plan

Date: 2026-05-13  
Branch: `parallel-blas`

This document specifies the plan for a hand-written, thread-parallel
overlay BLAS for the extended-precision targets (`kind10`, `kind16`,
`multifloats`). The migrated `${LIB_PREFIX}blas` archives produced by
the migrator are correctness-baselined (`1125/1125` differential
precision tests against a Netlib quad reference) but **serial and
unblocked** — they are 1:1 type-rewrites of Netlib's reference BLAS.
This project layers cache-blocked, OpenMP-parallel kernels on top
without touching the migrator pipeline.

## 1. Strategy

Build a **separate** static archive per target — `lib{e,q,m}blas_parallel.a`
— that supplies hand-written replacements for a small set of
high-leverage routines. The migrated archive continues to provide
everything else.

The overlay archive `lib{e,q,m}blas_parallel.a` is **always built**
(no opt-in gate). What `PARALLEL_BLAS` (default ON) controls is
whether the public CMake target `${LIB_PREFIX}blas` is rewired into a
**composite INTERFACE** that links the overlay *before* the migrated
archive. Concretely:

| target name                       | type      | role                                                |
|-----------------------------------|-----------|-----------------------------------------------------|
| `${LIB_PREFIX}blas_serial`        | STATIC    | the migrated archive (was `${LIB_PREFIX}blas`)      |
| `${LIB_PREFIX}blas_parallel`      | STATIC    | the hand-written C overlay (this project)           |
| `${LIB_PREFIX}blas`               | INTERFACE | composite: links `_parallel` then `_serial`         |
| `${LIB_PREFIX}blas_migrated`      | STATIC    | `objcopy`-renamed `_serial` (symbols → `*_migrated_`); test-only |

Effect: every existing `target_link_libraries(... ${LIB_PREFIX}blas)`
— in `${LIB_PREFIX}lapack`, `${LIB_PREFIX}scalapack`,
`${LIB_PREFIX}mumps`, downstream apps, the installed package — picks
up overlay symbols automatically. Linker resolves the overlay's
`egemm_` first; everything not in the overlay falls through to the
serial archive. No per-consumer link-line edits anywhere.

When `PARALLEL_BLAS=OFF`, `${LIB_PREFIX}blas` remains a plain STATIC
target pointing at the migrated archive — current behavior unchanged.
The overlay archive still builds and is exercised by its own
consistency / fuzz / bench tests, but no downstream consumer links
it. This mode exists for A/B'ing end-to-end performance.

The overlay source lives at `src/parallel_blas/<target>/` in the
project repo and is staged into the unified build tree by
`migrator stage`.

### Language

**C99-and-up, not Fortran.** Reasons:

- Better auto-vectorization where the type permits.
- Idiomatic `restrict`, alignment attributes, `__attribute__((hot))`.
- `#pragma omp parallel for` more directly than Fortran's `!$omp`.
- Easier control of packing buffer layout and prefetch.

Type mapping:

| target        | C type                  | notes                                |
|---------------|-------------------------|--------------------------------------|
| `kind10`      | `long double`           | x86-64 80-bit, ABI matches `REAL(KIND=10)` |
| `kind16`      | `__float128` (quadmath) | every op is a libquadmath call       |
| `multifloats` | multifloats C API       | upstream provides `float64x2_t` etc. |

### Fortran ↔ C boundary

Hand-written kernels are exported with the gfortran symbol mangling
that the migrated archive uses (lowercase + trailing `_`). The overlay
symbol *is* the BLAS entry point — no `iso_c_binding` shim.

- Scalar args by pointer (Fortran reference semantics).
- Arrays as raw pointers, column-major.
- Character args: pointer + hidden trailing `size_t` length per
  gfortran ABI. A common helper `static int parse_trans(const char *,
  size_t)` normalizes to `'N' | 'T' | 'C'`.

### Primary numerical oracle: the migrated archive

The hand-written parallel kernels are tested **against the migrated
BLAS at the same precision** (not against the Netlib quad reference).

- Same precision both sides → only divergence comes from reduction
  order (blocking, threading).
- **Tolerance: 10 ulp relative error** (user-set). Bit-exact across
  thread counts is explicitly NOT required.
- The existing `tests/blas/level{1,2,3}` precision suite (vs quad
  reference) is kept as a secondary regression check — runs
  automatically against whatever resolves `egemm_`, so flipping
  `PARALLEL_BLAS=ON` re-tests the overlay at the existing tolerance.

To allow both implementations to coexist in one test binary, the
migrated archive is renamed via `objcopy --redefine-syms` to
`*_migrated_` (`egemm_ → egemm_migrated_`, etc.). The test program
calls both `egemm_` (overlay) and `egemm_migrated_` (renamed migrated)
and diffs.

## 2. Routine priority

| phase   | routines                                           | rationale |
|---------|----------------------------------------------------|-----------|
| P2.a    | `egemm` / `ygemm`                                  | Only L3 routine that matters for downstream LAPACK/ScaLAPACK perf. Goto's blocked-and-packed algorithm. |
| P2.b    | `etrsm` / `ytrsm`                                  | Block-triangular solve, reduces to gemm on trailing block. |
| P2.c    | `esymm`/`yhemm`, `esyrk`/`yherk`, `esyr2k`/`yher2k`, `etrmm` | All reduce to gemm. |
| P3      | L2 `egemv` / `ygemv`                               | Memory-bound; OMP across rows, blocked for L2. |
| P4      | L1 `eaxpy` / `edot` / `enrm2`                      | Likely skip — OMP overhead exceeds arithmetic for typical sizes. Only if profile says so. |

Each phase only lands once fuzz + bench are green on the previous one.

## 3. Build integration

- Source tree: `src/parallel_blas/<target>/{egemm.c, etrsm.c, ...}` plus
  `src/parallel_blas/<target>/CMakeLists.txt` that defines the static
  archive `${LIB_PREFIX}blas_parallel`.
- Top-level `src/parallel_blas/CMakeLists.txt` gates on
  `option(PARALLEL_BLAS)` and dispatches to the per-target subdir.
- `cmake/CMakeLists.txt` adds `add_subdirectory(parallel_blas)` if the
  staged tree contains it.
- `migrator stage` copies `src/parallel_blas/` into the staging dir
  (one-line addition to `cmd_stage` in `src/migrator/__main__.py`).
- Compile flags: `-O3 -fopenmp -march=native -ffp-contract=fast`
  (env-overridable per kernel for the bench autotune step).

## 4. Test infrastructure

Two test layers, both always-on whenever `PARALLEL_BLAS=ON`.

### 4.1 Consistency suite — **new, primary**

`tests/blas_parallel/consistency/` — one driver per routine.

- Calls both `egemm_` (overlay) and `egemm_migrated_` (renamed migrated).
- Compares result matrices/vectors in `REAL(KIND=10)` (same precision
  both sides).
- Tolerance: `max_rel_err ≤ 10 * eps10` per element (kind10), with
  `eps10 ≈ 1.08e-19`. Equivalent thresholds for kind16, multifloats.
- Fails with the seed printed so failures are reproducible.

### 4.2 Precision regression — existing `tests/blas/level{1,2,3}`

Unchanged. Runs against whatever resolves the BLAS symbols. With
`PARALLEL_BLAS=ON` the overlay gets exercised; tolerance is the
existing `c * K * eps_quad` scheme. Catches any accuracy regression
the parallel kernel introduces relative to the Netlib quad oracle.

### 4.3 Fuzz harness — **new**

`tests/blas_parallel/fuzz/`. Per routine, parameterized:

- Random `M, N, K` mixing edges {0, 1, 2, MR-1, MR, MR+1, 2·MR, ...}
  with mid (32–256) and tail (≥KC, ≥MC, log-distributed).
- Random `LDA, LDB, LDC ≥ required` (padding stresses stride code).
- All TRANS/UPLO/DIAG combos.
- Random `alpha, beta` including {0, 1, -1, denormal, near-overflow}.
- Random `INCX, INCY` ∈ {1, 2, -1, -3} for L1/L2.
- Seed from `BLAS_FUZZ_SEED` env (default: time-based, logged).
  Fail-with-seed for repro.
- Iteration count: `BLAS_FUZZ_CASES` (default 500 fast, 50000 nightly
  via `ctest -L nightly`).
- Oracle selectable via `BLAS_FUZZ_ORACLE = migrated|quad|both`
  (default `both`).

### 4.4 Threading torture

Fuzz re-run at `OMP_NUM_THREADS ∈ {1, 2, 4, 8, all}`. Pass if every
thread count is within 10 ulp of every other (not bit-exact).

## 5. Benchmark harness — **new**

`tests/blas_parallel/bench/`. Per routine:

- Geometric size sweep (e.g. GEMM: 64, 128, 256, 512, 1024, 2048).
- Warm-up + N timed iterations; report median GFLOP/s, min, max.
- Block sizes `MC, KC, NC` taken from env (`EBLAS_MC`, etc.) so a
  Python launcher can grid-search.
- JSON output under `bench_reports/`, schema mirroring
  `precision_reports/`.
- A `scripts/bench_report.py` aggregator that plots GFLOP/s vs N and
  speedup vs threads, and compares to the serial migrated baseline
  (`PARALLEL_BLAS=OFF` build).

## 6. Kernel design (kind10, GEMM)

Goto's blocked-and-packed scheme:

- Outer blocks: `NC` (~1024 cols of B/C), `KC` (~256), `MC` (~256).
  Tunable via env in bench.
- **Pack** A panel (MC×KC, into MR-row stripes) and B block (KC×NC,
  into NR-column stripes) into 64-byte-aligned scratch
  (`posix_memalign` / `aligned_alloc`), one per thread.
- **Micro-kernel** MR×NR = 4×4 for real, 4×2 for complex. Accumulator
  block kept in scalar registers (no SIMD intrinsics — gfortran/gcc
  have no 80-bit vector ISA on x86-64). `static inline`,
  `__attribute__((hot))`.
- `restrict`-qualify all inner pointers.
- **OpenMP**: `#pragma omp parallel for schedule(static)` on the
  NC-loop. Per-thread packing scratch indexed by
  `omp_get_thread_num()`. Reduction order is *not* preserved across
  thread counts (allowed by the 10-ulp tolerance).

For **kind16**: same skeleton; every inner FMA is a libquadmath call.
No vectorization. Blocking still pays for cache + TLB reasons.

For **multifloats**: inner kernel becomes an error-free transform
sequence. Consider exposing a fused `dot_acc` from upstream
multifloats if available. OMP outside as before.

## 7. Phasing

| step | deliverable                                                                            | gate |
|------|----------------------------------------------------------------------------------------|------|
| **0**  | C build plumbing: `src/parallel_blas/` skeleton, `PARALLEL_BLAS` option, passthrough stub for `egemm_` | builds clean; existing 1125 unaffected when `PARALLEL_BLAS=OFF` |
| **0b** | `objcopy --redefine-syms` step that emits `lib${LIB_PREFIX}blas_migrated.a`; first test exe links both | passthrough call lands in migrated impl |
| **0c** | `cmd_stage` copies `src/parallel_blas/` into staging                                   | `migrator stage` produces a self-contained tree |
| **0d** | Composite rewiring of `${LIB_PREFIX}blas` (rename existing target → `_serial`, add INTERFACE shim) under `PARALLEL_BLAS=ON` | `${LIB_PREFIX}lapack` links unchanged, downstream consumers transparently pick up the overlay |
| 1     | Fuzz harness for `egemm` (runs against the migrated archive only — validates the harness itself) | finds zero false positives on migrated-vs-migrated |
| 2     | Bench harness for `egemm` (baseline numbers committed)                                 | reproducible GFLOP/s for migrated baseline |
| 3     | Serial blocked + packed `egemm` in C (no OMP yet)                                      | fuzz green; bench ≥ 5× serial Netlib at N=1024 |
| 4     | Add OMP to `egemm`                                                                     | fuzz green at all thread counts; ~linear scaling to N_core/2 |
| 5     | Block-size autotune sweep; commit tuned `MC/KC/NC`                                     | aggregator picks winners; merge to `develop` |
| 6     | Repeat for `ygemm` (complex)                                                           | fuzz + bench |
| 7     | `etrsm` / `etrmm` / `esymm` / `esyrk` family                                           | fuzz + bench per routine |
| 8     | Port to `kind16`                                                                       | reuse harness 1:1 |
| 9     | Port to `multifloats`                                                                  | reuse harness 1:1 |

Phases 1–4 stay on branch `parallel-blas`; merge to `develop` only
after step 5 lands and is green.

## 8. Tolerance & determinism

- Per-test pass criterion: `max(|r_overlay - r_migrated| / |r_migrated|)
  ≤ 10 * eps_target` (with `eps_target = epsilon(1.0_K)`).
- Bit-exact across thread counts NOT required — frees the OMP
  partitioning to choose any schedule.
- Determinism within a fixed `OMP_NUM_THREADS` and fixed schedule is
  guaranteed by the algorithm (no atomics, no dynamic schedule for the
  reduction loop).

## 9. Open questions / explicit defaults

| item                  | default                                                           |
|-----------------------|-------------------------------------------------------------------|
| branch                | `parallel-blas`; merge to `develop` after step 5                  |
| bench iterations      | 5 timed + 1 warm-up; fuzz fast=500, nightly=50000                 |
| determinism           | 10 ulp rel err; not bit-exact                                     |
| overlay packaging     | separate `lib*blas_parallel.a`; opt-in via `PARALLEL_BLAS=ON`     |
| Fortran ABI           | gfortran lowercase + trailing-`_` mangling, no `iso_c_binding` shim |
| symbol-aliasing trick | `objcopy --redefine-syms`; fallback `-Wl,--wrap` if it misbehaves |

## 10. References

- Project memory: [project-parallel-blas-design](../../../.claude/projects/-home-kyungminlee-Code-fortran-migrator/memory/project_parallel_blas_design.md)
- Goto, K., van de Geijn, R. — *Anatomy of High-Performance Matrix
  Multiplication* (ACM TOMS, 2008). Source of the
  MC/KC/NC + pack-A / pack-B layout.
- BLIS framework — practical reference implementation of the same
  algorithm, in C, for arbitrary scalar types.
