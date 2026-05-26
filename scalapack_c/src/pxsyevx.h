

/*
 * These macros define how C routines will be called.  ADD_ assumes that
 * they will be called by fortran, which expects C routines to have an
 * underscore postfixed to the name (Suns, and the Intel expect this).
 * NOCHANGE indicates that fortran will be calling, and that it expects
 * the name called by fortran to be identical to that compiled by the C
 * (RS6K's do this).  UPCASE says it expects C routines called by fortran
 * to be in all upcase (CRAY wants this). 
 */

#define ADD_       0
#define NOCHANGE   1
#define UPCASE     2
#define C_CALL     3

#ifdef UpCase
#define F77_CALL_C UPCASE
#endif

#ifdef NoChange
#define F77_CALL_C NOCHANGE
#endif

#ifdef Add_
#define F77_CALL_C ADD_
#endif

#ifndef F77_CALL_C
#define F77_CALL_C ADD_
#endif

#if (F77_CALL_C == ADD_)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine
 * No redefinition necessary to have following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call pdgemm(...)           void pdgemm_(...)
 *           call pegemm(...)           void pegemm_(...)
 *
 * This is the default.
 */

#endif

#if (F77_CALL_C == UPCASE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call pdgemm(...)           void PDGEMM(...)
 *           call pegemm(...)           void PEGEMM(...)
 */
                                                            /* TOOLS */
#define pdlasnbt_           PDLASNBT
#define pelasnbt_           PELASNBT
#define pdlachkieee_        PDLACHKIEEE
#define pelachkieee_        PELACHKIEEE
#define pdlaiectl_          PDLAIECTL
#define pelaiectl_          PELAIECTL
#define pdlaiectb_          PDLAIECTB
#define pelaiectb_          PELAIECTB

#define pslasnbt_           PSLASNBT
#define pelasnbt_           PELASNBT
#define pslachkieee_        PSLACHKIEEE
#define pelachkieee_        PELACHKIEEE
#define pslaiect_           PSLAIECT

#endif

#if (F77_CALL_C == NOCHANGE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * for following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call pdgemm(...)           void pdgemm(...)
 *           call pegemm(...)           void pegemm(...)
 */
                                                            /* TOOLS */
#define pdlasnbt_           pdlasnbt
#define pelasnbt_           pelasnbt
#define pdlachkieee_        pdlachkieee
#define pelachkieee_        pelachkieee
#define pdlaiectl_          pdlaiectl
#define pelaiectl_          pelaiectl
#define pdlaiectb_          pdlaiectb
#define pelaiectb_          pelaiectb

#define pslasnbt_           pslasnbt
#define pelasnbt_           pelasnbt
#define pslachkieee_        pslachkieee
#define pelachkieee_        pelachkieee
#define pslaiect_           pslaiect
#endif
