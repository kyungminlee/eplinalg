# epopenblas — OpenBLAS D/Z port to kind10 (REAL/COMPLEX(KIND=10))

Direct ports of OpenBLAS double-precision (D) real and double-complex (Z)
sources, retyped to 80-bit long double / `_Complex long double`. Builds a
separate static archive `${LIB_PREFIX}blas_epopenblas` (not wired into the
public `${LIB_PREFIX}blas` composite). Tested standalone against the
migrated baseline via `tests/blas_epopenblas/`.

Naming map (project convention, matches `blas/src/`):
- `d*`  (double real)             → `e*`
- `z*`  (double complex)          → `y*`
- `dz*` (real-from-complex)       → `ey*`
- `iz*` (integer-from-complex)    → `iy*`
- `zd*` (complex scaled by real)  → `ye*`

OpenBLAS source roots:
- `external/openblas/interface/<name>.c` — dispatch / argument check / shape
- `external/openblas/kernel/generic/<name>.c` — portable C kernel (no asm)
- `external/openblas/kernel/x86_64/d<name>.c` — D-specific kernel (some have
  microkernel SIMD; ignore — kind10 has no SSE/AVX path)
- `external/openblas/driver/level3/level3.c` — GotoBLAS three-level blocking
  driver for GEMM family (templated by data type via macros)

Per row:
- **status** — `todo` / `wip` / `done`
- **smoke** — fuzz passes, bench builds: `n` / `y`
- **bench-omp1**, **bench-omp4** — most recent GFLOPS overlay vs migrated

**`faithful` column:**
- `yes` — algorithm + threading mirror OpenBLAS (interface + kernel + thread driver structure all ported; OMP replaces blas_queue but matches the partitioning/reduction shape).
- `partial` — algorithm matches OpenBLAS but the threading uses a simpler OMP partition rather than the custom load-balanced width formula from the `*_thread.c` driver. Result is equivalent at OMP=1/coarse-grain; fine-grain partition differs.
- `no` — diverges from OpenBLAS (e.g., a different algorithm such as Blue's nrm2, or serial when OpenBLAS uses a thread driver).
- `n/a` — scalar / O(1) / no SMP path exists in OpenBLAS for this routine.

## Level 1 — real

| target  | openblas source        | status | smoke | faithful | bench-omp1            | bench-omp4            | notes                              |
|---------|------------------------|--------|-------|----------|-----------------------|-----------------------|------------------------------------|
| eaxpy   | interface/axpy.c       | done   | y     | yes      | 1.40/1.31 (×1.07)     | 5.14/1.31 (×3.94)     | OMP block-partition (mirrors axpy_thread implicit pattern) |
| ecopy   | interface/copy.c       | done   | y     | yes      | 26.43/26.32 (×1.00)   | 92.11/25.32 (×3.64)   | OMP block-partition |
| eswap   | interface/swap.c       | done   | y     | yes      | 24.87/25.99 (×0.96)   | 90.43/25.59 (×3.53)   | OMP block-partition |
| escal   | interface/scal.c       | done   | y     | yes      | 0.82/0.80 (×1.01)     | 2.92/0.79 (×3.69)     | OMP block-partition |
| edot    | interface/dot.c        | done   | y     | yes      | 2.33/2.19 (×1.06)     | 7.38/1.95 (×3.78)     | seq 5-way unroll, 1 acc (fuzz tol) |
| enrm2   | interface/nrm2.c       | done   | y     | **no**   | n/a (Blue scl)        | n/a (Blue scl)        | Blue's algorithm — diverges from OpenBLAS's scaled-SSQ |
| easum   | interface/asum.c       | done   | y     | yes      | 1.97/1.42 (×1.39)     | 5.77/1.42 (×4.07)     | fabsl + 6-acc                      |
| erot    | interface/rot.c        | done   | y     | yes      | 1.96/1.83 (×1.07)     | 6.36/1.77 (×3.60)     | OMP block-partition                |
| erotg   | interface/rotg.c       | done   | y     | n/a      | 0.50/0.46 (×1.09)     | 0.12/0.12 (×0.97)     | scalar; OpenBLAS has no SMP path   |
| erotm   | interface/rotm.c       | done   | y     | yes      | 1.31/1.32 (×1.00)     | 4.27/1.32 (×3.25)     | OMP block-partition, unswitched per dflag |
| erotmg  | interface/rotmg.c      | done   | y     | n/a      | 0.48/0.59 (×0.82)     | 0.68/0.53 (×1.29)     | scalar; no SMP                     |
| ieamax  | interface/imax.c       | done   | y     | yes      | 0.60/0.53 (×1.12)     | 2.65/0.52 (×5.09)     |                                    |

## Level 1 — complex

| target  | openblas source        | status | smoke | faithful | bench-omp1            | bench-omp4            | notes                              |
|---------|------------------------|--------|-------|----------|-----------------------|-----------------------|------------------------------------|
| yaxpy   | interface/axpy.c       | done   | y     | yes      | 2.00/2.00 (×1.00)     | 7.78/2.01 (×3.87)     |                                    |
| ycopy   | interface/copy.c       | done   | y     | yes      | 24.43/24.64 (×0.99)   | 93.82/25.51 (×3.68)   |                                    |
| yswap   | interface/swap.c       | done   | y     | yes      | 25.00/24.37 (×1.03)   | 76.73/23.71 (×3.24)   |                                    |
| yscal   | interface/scal.c       | done   | y     | yes      | 2.00/1.87 (×1.07)     | 7.11/1.97 (×3.61)     | manual re/im, xi*ar first (-1 fxch)|
| yescal  | interface/scal.c       | done   | y     | yes      | 2.18/2.03 (×1.07)     | 8.01/2.11 (×3.79)     | explicit 2-fmul pair walk          |
| ydotu   | interface/dot.c        | done   | y     | yes      | 3.16/3.11 (×1.02)     | 10.97/3.03 (×3.62)    |                                    |
| ydotc   | interface/dot.c        | done   | y     | yes      | 2.79/2.77 (×1.01)     | 10.72/2.79 (×3.85)    |                                    |
| eynrm2  | interface/nrm2.c       | done   | y     | **no**   | n/a (Blue scl)        | n/a (Blue scl)        | Blue's algorithm — diverges from OpenBLAS |
| eyasum  | interface/asum.c       | done   | y     | yes      | 1.24/0.78 (×1.59)     | 4.66/0.78 (×6.00)     |                                    |
| yerot   | interface/rot.c (zrot) | done   | y     | yes      | 1.96/1.01 (×1.94)     | 6.79/0.99 (×6.84)     | OMP block-partition                |
| yrotg   | interface/rotg.c       | done   | y     | n/a      | 0.05/0.05 (×0.96)     | 0.05/0.05 (×0.97)     | scalar; no SMP                     |
| iyamax  | interface/imax.c       | done   | y     | yes      | 1.15/0.95 (×1.22)     | 4.06/0.93 (×4.37)     |                                    |
| ecabs1  | kernel/generic/cabs.c  | done   | y     | n/a      | n/a (scalar)          | n/a (scalar)          | scalar O(1)                        |

## Level 2 — real  (bench at N=512, kernel-only, unit strides)

| target  | openblas source        | status | smoke | faithful | bench-omp1            | bench-omp4            | notes |
|---------|------------------------|--------|-------|----------|-----------------------|-----------------------|-------|
| egemv   | interface/gemv.c       | done   | y     | yes      | 1.69/1.88 (×0.90)     | 5.67/1.83 (×3.09)     | M-split + gemv_thread split_x fallback (m < 4·nthreads, m·n large) |
| eger    | interface/ger.c        | done   | y     | yes      | 1.08/1.10 (×0.99)     | 4.72/1.20 (×3.94)     | N-block partition + ger_thread x-buffer-copy for incx≠1 |
| esymv   | interface/symv.c       | done   | y     | yes      | 1.73/1.96 (×0.88)     | 6.59/1.84 (×3.57)     | symv_thread sqrt-partition (mask=3 min-4) + AXPY-chain reduce |
| esyr    | interface/syr.c        | done   | y     | yes      | 1.01/1.00 (×1.01)     | 2.76/0.99 (×2.79)     | syr_thread sqrt-partition (mask=7 min-16); UPPER reverse-mapped |
| esyr2   | interface/syr2.c       | done   | y     | yes      | 1.85/1.87 (×0.99)     | 6.28/1.67 (×3.76)     | syr2_thread sqrt-partition (mask=7 min-16); UPPER reverse-mapped |
| espr    | interface/spr.c        | done   | y     | yes      | 1.20/1.20 (×1.00)     | 3.95/1.05 (×3.74)     | spr_thread sqrt-partition (mask=7 min-16); packed col offsets |
| espr2   | interface/spr2.c       | done   | y     | yes      | 0.95/1.85 (×0.51)     | 6.31/1.68 (×3.75)     | spr2_thread sqrt-partition (mask=7 min-16); packed col offsets |
| espmv   | interface/spmv.c       | done   | y     | yes      | 1.69/1.87 (×0.90)     | 6.73/1.87 (×3.60)     | spmv_thread sqrt-partition + AXPY-chain reduce; packed col offsets |
| esbmv   | interface/sbmv.c       | done   | y     | yes      | 1.71/1.69 (×1.01)     | 1.71/1.69 (×1.01)     | sbmv_thread port: load-balanced cols + y_priv reduce + alpha-AXPY |
| egbmv   | interface/gbmv.c       | done   | y     | yes      | 1.00/1.16 (×0.86)     | 1.02/1.16 (×0.88)     | gbmv_thread port: col-partition + y_priv reduce + alpha-AXPY |
| etrmv   | interface/trmv.c       | done   | y     | yes      | 1.20/1.03 (×1.17)     | 1.21/1.05 (×1.15)     | trmv_thread port: sqrt-partition + DTB_ENTRIES tile + per-thread slot + AXPY-reduce |
| etrsv   | interface/trsv.c       | done   | y     | yes      | 1.19/1.22 (×0.97)     | 1.20/1.20 (×1.00)     | serial — matches OpenBLAS (no trsv_thread.c) |
| etbmv   | interface/tbmv.c       | done   | y     | yes      | 0.98/1.12 (×0.88)     | 0.98/1.06 (×0.92)     | tbmv_thread port: sqrt/even partition + private slots + full AXPY-reduce |
| etbsv   | interface/tbsv.c       | done   | y     | yes      | 0.92/1.00 (×0.92)     | 0.97/1.01 (×0.96)     | serial — matches OpenBLAS (no tbsv_thread.c) |
| etpmv   | interface/tpmv.c       | done   | y     | yes      | 1.02/1.05 (×0.97)     | 1.02/1.04 (×0.98)     | tpmv_thread port: sqrt-partition + per-thread slot + AXPY-reduce |
| etpsv   | interface/tpsv.c       | done   | y     | yes      | 1.18/1.14 (×1.03)     | 1.20/1.20 (×1.00)     | serial — matches OpenBLAS (no tpsv_thread.c) |

## Level 2 — complex (general / sym / triangular / band / packed)

| target  | openblas source        | status | smoke | faithful | bench-omp1            | bench-omp4            | notes |
|---------|------------------------|--------|-------|----------|-----------------------|-----------------------|-------|
| ygemv   | interface/gemv.c       | done   | y     | yes      | 1.77/2.67 (×0.66)     | 6.76/2.61 (×2.59)     | M-split + gemv_thread split_x fallback (m < 4·nthreads, m·n large) |
| ygeru   | interface/ger.c (zgeru)| done   | y     | yes      | 1.89/1.90 (×1.00)     | 7.50/1.89 (×3.97)     | N-block partition + ger_thread x-buffer-copy for incx≠1 |
| ygerc   | interface/ger.c (zgerc)| done   | y     | yes      | 1.84/1.80 (×1.02)     | 7.44/1.89 (×3.92)     | N-block partition + conjg y; same shape as ygeru |
| ygbmv   | interface/zgbmv.c      | done   | y     | yes      | 1.94/1.95 (×1.00)     | 1.85/1.95 (×0.95)     | gbmv_thread port: col-partition + y_priv reduce + alpha-AXPY (N/T/C) |
| ytrmv   | interface/trmv.c       | done   | y     | yes      | 1.96/1.97 (×0.99)     | 1.97/1.89 (×1.04)     | trmv_thread port: sqrt-partition + DTB_ENTRIES tile + per-thread slot + AXPY-reduce (N/T/C) |
| ytrsv   | interface/trsv.c       | done   | y     | yes      | 0.55/0.53 (×1.03)     | 0.54/0.55 (×0.98)     | serial — matches OpenBLAS |
| ytbmv   | interface/tbmv.c       | done   | y     | yes      | 1.96/1.87 (×1.05)     | 1.88/1.86 (×1.01)     | tbmv_thread port: sqrt/even partition + private slots + full AXPY-reduce (N/T/C) |
| ytbsv   | interface/tbsv.c       | done   | y     | yes      | 0.06/0.06 (×1.01)     | 0.06/0.06 (×0.99)     | serial — matches OpenBLAS |
| ytpmv   | interface/tpmv.c       | done   | y     | yes      | 1.99/2.04 (×0.97)     | 2.01/2.02 (×0.99)     | tpmv_thread port: sqrt-partition + per-thread slot + AXPY-reduce (N/T/C) |
| ytpsv   | interface/tpsv.c       | done   | y     | yes      | 0.56/0.58 (×0.97)     | 0.56/0.57 (×0.98)     | serial — matches OpenBLAS |

## Level 2 — complex Hermitian

| target  | openblas source        | status | smoke | faithful | bench-omp1            | bench-omp4            | notes |
|---------|------------------------|--------|-------|----------|-----------------------|-----------------------|-------|
| yhemv   | interface/zhemv.c      | done   | y     | yes      | 2.66/1.99 (×1.34)     | 10.36/2.33 (×4.45)    | symv_thread+HEMV sqrt-partition + AXPY-chain reduce; conj A^H reflection |
| yhbmv   | interface/zhbmv.c      | done   | y     | yes      | 2.56/2.01 (×1.27)     | 2.58/2.01 (×1.28)     | sbmv_thread+HEMV port: DOTC band + REAL diagonal |
| yher    | interface/zher.c       | done   | y     | yes      | 1.69/1.69 (×1.00)     | 6.89/1.87 (×3.68)     | syr_thread+HEMV sqrt-partition (mask=7 min-16); alpha REAL |
| yher2   | interface/zher2.c      | done   | y     | yes      | 2.46/2.23 (×1.10)     | 9.69/2.76 (×3.51)     | syr2_thread+HEMV sqrt-partition (mask=7 min-16) |
| yhpmv   | interface/zhpmv.c      | done   | y     | yes      | 2.66/1.98 (×1.34)     | 5.55/1.80 (×3.09)     | spmv_thread+HEMV sqrt-partition + AXPY-chain reduce; packed |
| yhpr    | interface/zhpr.c       | done   | y     | yes      | 1.95/1.95 (×1.00)     | 7.09/1.86 (×3.81)     | spr_thread+HEMV sqrt-partition (mask=7 min-16); packed |
| yhpr2   | interface/zhpr2.c      | done   | y     | yes      | 2.79/2.46 (×1.13)     | 10.52/2.92 (×3.61)    | spr2_thread+HEMV sqrt-partition (mask=7 min-16); packed |

## Level 3 — real

| target  | openblas source        | status | smoke | faithful | bench-omp1 | bench-omp4 | notes |
|---------|------------------------|--------|-------|----------|------------|------------|-------|
| egemm   | interface/gemm.c       | todo   | n     |          |            |            |       |
| egemmtr | interface/gemmt.c      | todo   | n     |          |            |            |       |
| esymm   | interface/symm.c       | todo   | n     |          |            |            |       |
| esyrk   | interface/syrk.c       | todo   | n     |          |            |            |       |
| esyr2k  | interface/syr2k.c      | todo   | n     |          |            |            |       |
| etrmm   | interface/trmm.c       | todo   | n     |          |            |            |       |
| etrsm   | interface/trsm.c       | todo   | n     |          |            |            |       |

## Level 3 — complex (general / sym / triangular)

| target  | openblas source        | status | smoke | faithful | bench-omp1 | bench-omp4 | notes |
|---------|------------------------|--------|-------|----------|------------|------------|-------|
| ygemm   | interface/gemm.c       | todo   | n     |          |            |            |       |
| ygemmtr | interface/gemmt.c      | todo   | n     |          |            |            |       |
| ysymm   | interface/symm.c       | todo   | n     |          |            |            |       |
| ysyrk   | interface/syrk.c       | todo   | n     |          |            |            |       |
| ysyr2k  | interface/syr2k.c      | todo   | n     |          |            |            |       |
| ytrmm   | interface/trmm.c       | todo   | n     |          |            |            |       |
| ytrsm   | interface/trsm.c       | todo   | n     |          |            |            |       |

## Level 3 — complex Hermitian

| target  | openblas source        | status | smoke | faithful | bench-omp1 | bench-omp4 | notes |
|---------|------------------------|--------|-------|----------|------------|------------|-------|
| yhemm   | interface/zhemm.c?     | todo   | n     |          |            |            |       |
| yherk   | interface/zherk.c?     | todo   | n     |          |            |            |       |
| yher2k  | interface/zher2k.c?    | todo   | n     |          |            |            |       |

(Some Hermitian L3 sources live in `driver/level3/` rather than `interface/`
because the dispatcher hands directly to the GEMM driver with a flag —
verify the actual entry points when porting.)

---

Total: 75 routines (matches `blas/src/` migrated set).

How to update this file as ports land: change `status` to `wip` when starting,
`done` when fuzz passes; flip `smoke` to `y` once the CMake list is wired and
the test driver builds; fill `bench-omp1`/`bench-omp4` with the most recent
median GFLOPS, formatted `ov / mig (×speedup)`, e.g. `12.3 / 8.1 (×1.52)`.
Use `notes` to record kernel-shape decisions worth remembering (algorithm
variant, what was tried and reverted, why a faithful straight port was
swapped out): e.g. `fabsl + 6-acc (was ternary+4-acc)` or `Blue's algorithm
(Algorithm 978)`.

Bench numbers are **kernel-only**: the perf harness measures each BLAS
call individually via `perf_now_s()` and excludes the per-iter memcpy
reset from the timed window. Previously the SCAL/COPY/SWAP family
appeared to cap at ~2× OMP=4 speedup, but that was a serial-memcpy
Amdahl artifact in the harness, not the kernel. The peak GFLOPS columns
above now reflect actual kernel scaling (4 threads on this 6-core box).

Note for the **rot family** (erot, erotm, yerot): the earlier
"DRAM-BW bound" diagnosis (2026-05-22) was a measurement artifact —
when re-benched in a clean stage-e build on the same 6-core
DDR4-2666 box, OMP=4 scales 3.2–3.6× on top of OMP=1 for erot
and erotm at N=65536, and yerot stacks the overlay's ~×1.9 OMP=1
margin with another ×3.5 from the 4 cores. The OMP wiring matters
even on a consumer box, not just for the HPC deployment target.
