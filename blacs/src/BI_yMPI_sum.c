#include "Bdef.h"
void BI_yMPI_sum(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_yvvsum(Int, char *, char *);
   BI_yvvsum(*N, inout, in);
}
