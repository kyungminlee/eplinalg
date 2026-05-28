/*
 * BI_eMPI_sum -- user-defined MPI reduction op for the real-precision
 * BLACS gsum2d combine. Mirrors the upstream ``BI_zMPI_sum.c`` pattern,
 * which the complex ``zgsum2d_`` already uses, but for the real path.
 *
 * The d-prefix BLACS gsum2d_ default-topology branch was historically
 * the only BLACS reduction that called the stock ``MPI_SUM``. On
 * extended-precision targets that resolves against ``MPI_REAL16`` /
 * ``MPI_LONG_DOUBLE``, neither of which Intel MPI 2021.17 implements
 * ``MPI_SUM`` for. Switching to a user-op closes the gap by performing
 * the per-element addition in user code (native ``+=`` on the target
 * real type after migration), bypassing the broken builtin entirely.
 *
 * See the B01 description in recipes/blacs.yaml for full context.
 */
#include "Bdef.h"
void BI_eMPI_sum(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_evvsum(Int, char *, char *);
   BI_evvsum(*N, inout, in);
}
