#include "Bdef.h"

#if (INTFACE == C_CALL)
EREAL Cecputime00(void)
#else
F_DOUBLE_FUNC ecputime00_(void)
#endif
{
   return(-1.0);
}
