#include "Bdef.h"
void BI_evvsum(Int N, char *vec1, char *vec2)
{
   EREAL *v1=(EREAL*)vec1, *v2=(EREAL*)vec2;
   Int k;
   for (k=0; k < N; k++) v1[k] += v2[k];
}
