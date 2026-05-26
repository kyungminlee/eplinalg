#include "Bdef.h"

void BI_yMPI_amn(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_yvvamn(Int, char *, char *);
   extern BLACBUFF BI_AuxBuff;

   BI_yvvamn(BI_AuxBuff.Len, inout, in);
}
