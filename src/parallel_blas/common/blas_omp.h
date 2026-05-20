/* blas_omp.h — shared OpenMP runtime helpers for the parallel BLAS overlay.
 *
 * Caches `omp_get_max_threads()` at the first overlay call so we don't pay
 * a libgomp ICV lookup on every dispatch. The cache lives in blas_omp.c
 * as a SINGLE static — earlier this was a `static inline` body in the
 * header, which gave each translation unit its own copy. With per-TU
 * caches, the kernel whose first call arrived at OMP_NUM_THREADS=1 (e.g.
 * xgemv during a perf bench's cur_1 phase) would lock itself to 1 even
 * after the caller switched to omp_set_num_threads(N>1) for the next
 * timing, silently disabling parallelism in that kernel.
 *
 * BLAS callers virtually never change OMP_NUM_THREADS mid-run in
 * production, and the migrated reference path doesn't honour mid-run
 * changes either, so a one-shot cache is correct for real use.
 */
#ifndef PARALLEL_BLAS_BLAS_OMP_H
#define PARALLEL_BLAS_BLAS_OMP_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _OPENMP
#include <omp.h>
extern int blas_omp_max_threads(void);
#else
static inline int blas_omp_max_threads(void) { return 1; }
#endif

#ifdef __cplusplus
}
#endif

#endif /* PARALLEL_BLAS_BLAS_OMP_H */
