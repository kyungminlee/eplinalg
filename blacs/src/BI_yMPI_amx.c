#include "Bdef.h"

void BI_yMPI_amx(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_yvvamx(Int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_yvvamx(BI_AuxBuff.Len, inout, in);
}
