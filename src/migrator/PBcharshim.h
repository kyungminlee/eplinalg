#ifndef EP_PBCHARSHIM_H
#define EP_PBCHARSHIM_H
/* ---------------------------------------------------------------------
*
*  PBcharshim.h -- hidden-CHARACTER-length ABI bridge for the migrated
*                  gfortran extended-precision BLAS facades.
*
*  The type-generic PBLAS drivers reach the BLAS through PBTYP_T function
*  pointers (GEMM_T, HEMV_T, ...) whose C prototypes use
*
*        typedef char * F_CHAR_T;   #define C2F_CHAR(a) (a)
*
*  i.e. they pass NO hidden Fortran CHARACTER-length arguments and reserve
*  no stack for them.  A migrated gfortran facade (ygemm_, ytrsm_, ...) is a
*  genuine Fortran routine: per the SysV AMD64 ABI it may freely use its
*  incoming hidden-length slots as scratch, and gfortran's codegen for these
*  routines does exactly that (e.g. `mov %edx,0x1a0(%rsp)`).  Because the C
*  caller never provisioned those slots, the writes clobber the caller's live
*  locals -> wild pointers -> SIGSEGV / garbage MPI ranks.  Reference/MKL BLAS
*  merely never happen to touch those slots.
*
*  A trampoline with the driver's exact no-length signature that forwards to
*  the facade with the hidden lengths appended fixes this: the facade's
*  scratch stores land in the trampoline's own outgoing-argument area, and the
*  CHARACTER lengths (always 1 here -- every argument is CHARACTER*1) are
*  correct.  Confined to the extended typesets; the MKL s/c/d/z path is
*  untouched.
*
*  Include AFTER pblas.h (needs F_CHAR_T, Int, F_VOID_FCT and the *_T
*  typedefs).  Requires the length-carrying (non-Cray) F_CHAR_T; the Cray
*  _fcd path already carries lengths and must not use these shims.
*
*  ---------------------------------------------------------------------
*/
#include <stddef.h>

/*
*  Length-carrying function-pointer types: identical to the pblas.h callback
*  typedefs but with one trailing size_t per hidden CHARACTER-length argument.
*  Shapes are shared across callbacks with the same layout (e.g. GEMV/AGEMV,
*  SYMV/ASYMV/HEMV/AHEMV, TRMV/TRSV, SYR/HER, ...).
*/
typedef F_VOID_FCT (*EPL_TZPAD_T)   ( F_CHAR_T, F_CHAR_T, Int*, Int*, Int*,
                                      char*, char*, char*, Int*,
                                      size_t, size_t );
typedef F_VOID_FCT (*EPL_TZPADCPY_T)( F_CHAR_T, F_CHAR_T, Int*, Int*, Int*,
                                      char*, Int*, char*, Int*,
                                      size_t, size_t );
typedef F_VOID_FCT (*EPL_TZSCAL_T)  ( F_CHAR_T, Int*, Int*, Int*, char*,
                                      char*, Int*, size_t );
typedef F_VOID_FCT (*EPL_GEMV_T)    ( F_CHAR_T, Int*, Int*, char*, char*,
                                      Int*, char*, Int*, char*, char*, Int*,
                                      size_t );
typedef F_VOID_FCT (*EPL_SYMV_T)    ( F_CHAR_T, Int*, char*, char*, Int*,
                                      char*, Int*, char*, char*, Int*,
                                      size_t );
typedef F_VOID_FCT (*EPL_TRMV_T)    ( F_CHAR_T, F_CHAR_T, F_CHAR_T, Int*,
                                      char*, Int*, char*, Int*,
                                      size_t, size_t, size_t );
typedef F_VOID_FCT (*EPL_ATRMV_T)   ( F_CHAR_T, F_CHAR_T, F_CHAR_T, Int*,
                                      char*, char*, Int*, char*, Int*,
                                      char*, char*, Int*,
                                      size_t, size_t, size_t );
typedef F_VOID_FCT (*EPL_SYR_T)     ( F_CHAR_T, Int*, char*, char*, Int*,
                                      char*, Int*, size_t );
typedef F_VOID_FCT (*EPL_SYR2_T)    ( F_CHAR_T, Int*, char*, char*, Int*,
                                      char*, Int*, char*, Int*, size_t );
typedef F_VOID_FCT (*EPL_GEMM_T)    ( F_CHAR_T, F_CHAR_T, Int*, Int*, Int*,
                                      char*, char*, Int*, char*, Int*,
                                      char*, char*, Int*, size_t, size_t );
typedef F_VOID_FCT (*EPL_SYMM_T)    ( F_CHAR_T, F_CHAR_T, Int*, Int*, char*,
                                      char*, Int*, char*, Int*, char*, char*,
                                      Int*, size_t, size_t );
typedef F_VOID_FCT (*EPL_SYRK_T)    ( F_CHAR_T, F_CHAR_T, Int*, Int*, char*,
                                      char*, Int*, char*, char*, Int*,
                                      size_t, size_t );
typedef F_VOID_FCT (*EPL_SYR2K_T)   ( F_CHAR_T, F_CHAR_T, Int*, Int*, char*,
                                      char*, Int*, char*, Int*, char*, char*,
                                      Int*, size_t, size_t );
typedef F_VOID_FCT (*EPL_TRMM_T)    ( F_CHAR_T, F_CHAR_T, F_CHAR_T, F_CHAR_T,
                                      Int*, Int*, char*, char*, Int*, char*,
                                      Int*, size_t, size_t, size_t, size_t );

/*
*  Per-shape trampoline generators.  EP_MK_<SHAPE>(tag, fn) defines a static
*  function PB_shim_##tag with the driver's no-length signature that forwards
*  to the Fortran facade `fn` with the hidden CHARACTER lengths (all 1)
*  appended.  Assign it to the PBTYP_T field cast to the field's own typedef.
*/
#define EP_MK_TZPAD(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, F_CHAR_T A2, Int *A3, Int *A4, \
   Int *A5, char *A6, char *A7, char *A8, Int *A9 ) \
{ ((EPL_TZPAD_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8,A9, (size_t)1,(size_t)1 ); }

#define EP_MK_TZPADCPY(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, F_CHAR_T A2, Int *A3, Int *A4, \
   Int *A5, char *A6, Int *A7, char *A8, Int *A9 ) \
{ ((EPL_TZPADCPY_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8,A9, (size_t)1,(size_t)1 ); }

#define EP_MK_TZSCAL(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, Int *A2, Int *A3, Int *A4, char *A5, \
   char *A6, Int *A7 ) \
{ ((EPL_TZSCAL_T)(fn))( A1,A2,A3,A4,A5,A6,A7, (size_t)1 ); }

#define EP_MK_GEMV(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, Int *A2, Int *A3, char *A4, char *A5, \
   Int *A6, char *A7, Int *A8, char *A9, char *A10, Int *A11 ) \
{ ((EPL_GEMV_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11, (size_t)1 ); }

#define EP_MK_SYMV(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, Int *A2, char *A3, char *A4, Int *A5, \
   char *A6, Int *A7, char *A8, char *A9, Int *A10 ) \
{ ((EPL_SYMV_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8,A9,A10, (size_t)1 ); }

#define EP_MK_TRMV(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, F_CHAR_T A2, F_CHAR_T A3, Int *A4, \
   char *A5, Int *A6, char *A7, Int *A8 ) \
{ ((EPL_TRMV_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8, \
                      (size_t)1,(size_t)1,(size_t)1 ); }

#define EP_MK_ATRMV(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, F_CHAR_T A2, F_CHAR_T A3, Int *A4, \
   char *A5, char *A6, Int *A7, char *A8, Int *A9, char *A10, char *A11, \
   Int *A12 ) \
{ ((EPL_ATRMV_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12, \
                       (size_t)1,(size_t)1,(size_t)1 ); }

#define EP_MK_SYR(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, Int *A2, char *A3, char *A4, Int *A5, \
   char *A6, Int *A7 ) \
{ ((EPL_SYR_T)(fn))( A1,A2,A3,A4,A5,A6,A7, (size_t)1 ); }

#define EP_MK_SYR2(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, Int *A2, char *A3, char *A4, Int *A5, \
   char *A6, Int *A7, char *A8, Int *A9 ) \
{ ((EPL_SYR2_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8,A9, (size_t)1 ); }

#define EP_MK_GEMM(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, F_CHAR_T A2, Int *A3, Int *A4, \
   Int *A5, char *A6, char *A7, Int *A8, char *A9, Int *A10, char *A11, \
   char *A12, Int *A13 ) \
{ ((EPL_GEMM_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13, \
                      (size_t)1,(size_t)1 ); }

#define EP_MK_SYMM(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, F_CHAR_T A2, Int *A3, Int *A4, \
   char *A5, char *A6, Int *A7, char *A8, Int *A9, char *A10, char *A11, \
   Int *A12 ) \
{ ((EPL_SYMM_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12, \
                      (size_t)1,(size_t)1 ); }

#define EP_MK_SYRK(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, F_CHAR_T A2, Int *A3, Int *A4, \
   char *A5, char *A6, Int *A7, char *A8, char *A9, Int *A10 ) \
{ ((EPL_SYRK_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8,A9,A10, \
                      (size_t)1,(size_t)1 ); }

#define EP_MK_SYR2K(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, F_CHAR_T A2, Int *A3, Int *A4, \
   char *A5, char *A6, Int *A7, char *A8, Int *A9, char *A10, char *A11, \
   Int *A12 ) \
{ ((EPL_SYR2K_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12, \
                       (size_t)1,(size_t)1 ); }

#define EP_MK_TRMM(tag,fn) \
static void PB_shim_##tag( F_CHAR_T A1, F_CHAR_T A2, F_CHAR_T A3, F_CHAR_T A4, \
   Int *A5, Int *A6, char *A7, char *A8, Int *A9, char *A10, Int *A11 ) \
{ ((EPL_TRMM_T)(fn))( A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11, \
                      (size_t)1,(size_t)1,(size_t)1,(size_t)1 ); }

#define EP_SHIM(tag) PB_shim_##tag

#endif /* EP_PBCHARSHIM_H */
