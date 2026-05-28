#include "Bdef.h"
void BI_eMPI_amn2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_evvamn2(Int, char *, char *);
   BI_evvamn2(*N, inout, in);
}
