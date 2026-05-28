#include "Bdef.h"
void BI_yMPI_amn2(void *in, void *inout, MpiInt *N, MPI_Datatype *dtype)
{
   void BI_yvvamn2(Int, char *, char *);
   BI_yvvamn2(*N, inout, in);
}
