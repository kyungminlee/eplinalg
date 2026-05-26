#include "Bdef.h"
void BI_eMPI_amx2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_evvamx2(Int, char *, char *);
   BI_evvamx2(*N, inout, in);
}
