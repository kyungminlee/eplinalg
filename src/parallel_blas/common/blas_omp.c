/* blas_omp.c — single definition of the cached omp_get_max_threads
 * wrapper. See blas_omp.h for context on why the cache is centralized
 * rather than per-translation-unit static. */

#include "blas_omp.h"

#ifdef _OPENMP
int blas_omp_max_threads(void)
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
#endif
