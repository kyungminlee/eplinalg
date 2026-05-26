#include "Bdef.h"
void BI_yMPI_amx2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_yvvamx2(Int, char *, char *);
   BI_yvvamx2(*N, inout, in);
}
