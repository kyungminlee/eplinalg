/* blas_omp.h — shared OpenMP runtime helpers for the parallel BLAS overlay.
 *
 * The overlays were paying an `omp_get_max_threads()` libgomp call on every
 * BLAS invocation just to decide whether to fork. At small problem sizes
 * (N=64 L2 matrices, ~µs of FP work) this added 10-15% overhead. This
 * header caches the value after first use.
 */
#ifndef PARALLEL_BLAS_BLAS_OMP_H
#define PARALLEL_BLAS_BLAS_OMP_H

#ifdef _OPENMP
#include <omp.h>

/* Cached omp_get_max_threads(). Returns 1 in non-OMP builds. The cache
 * is process-global; OpenMP allows OMP_NUM_THREADS to change at runtime
 * but BLAS callers virtually never do, and the migrated reference path
 * doesn't honour mid-run changes either. */
static inline int blas_omp_max_threads(void)
{
    static int cached = 0;  /* 0 = uninitialized */
    int v = __atomic_load_n(&cached, __ATOMIC_RELAXED);
    if (__builtin_expect(v == 0, 0)) {
        v = omp_get_max_threads();
        if (v < 1) v = 1;
        __atomic_store_n(&cached, v, __ATOMIC_RELAXED);
    }
    return v;
}
#else
static inline int blas_omp_max_threads(void) { return 1; }
#endif

#endif /* PARALLEL_BLAS_BLAS_OMP_H */
