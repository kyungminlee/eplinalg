#include "Bdef.h"

void BI_eMPI_amx(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_evvamx(Int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_evvamx(BI_AuxBuff.Len, inout, in);
}
