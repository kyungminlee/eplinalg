#include "Bdef.h"
void BI_evvamx(Int N, char *vec1, char *vec2)
{
   EREAL *v1=(EREAL*)vec1, *v2=(EREAL*)vec2;
   EREAL diff;
   BI_DistType *dist1, *dist2;
   Int i, k;

   k = N * sizeof(EREAL);
   i = k % sizeof(BI_DistType);
   if (i) k += sizeof(BI_DistType) - i;
   dist1 = (BI_DistType *) &vec1[k];
   dist2 = (BI_DistType *) &vec2[k];

   for (k=0; k < N; k++)
   {
      diff = Rabs(v1[k]) - Rabs(v2[k]);
      if (diff < 0)
      {
         v1[k] = v2[k];
         dist1[k] = dist2[k];
      }
      else if (diff == 0)
      {
         if (dist1[k] > dist2[k])
         {
            v1[k] = v2[k];
            dist1[k] = dist2[k];
         }
      }
   }
}
