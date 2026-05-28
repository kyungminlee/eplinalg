#include "Bdef.h"

void BI_eMPI_amn(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_evvamn(Int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_evvamn(BI_AuxBuff.Len, inout, in);
}
