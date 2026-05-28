#include "Bdef.h"

#if (INTFACE == C_CALL)
EREAL Cewalltime00(void)
#else
F_DOUBLE_FUNC ewalltime00_(void)
#endif
{
   return(MPI_Wtime());
}
