#include "Bdef.h"
void BI_evvamx2(Int N, char *vec1, char *vec2)
{
   Int k;
   EREAL *v1=(EREAL*)vec1, *v2=(EREAL*)vec2;
   EREAL diff;

   for (k=0; k != N; k++)
   {
      diff = Rabs(v1[k]) - Rabs(v2[k]);
      if (diff < 0) v1[k] = v2[k];
      else if (diff == 0) if (v1[k] < v2[k]) v1[k] = v2[k];
   }
}
