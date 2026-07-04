#include "Bdef.h"

/* kind16 (REAL(KIND=16) / __float128) lazy-init hook. The migrated MUMPS
 * reduces through custom MPI ops (MPI_QQ_SUM / MPI_XX_* on the standard
 * MPI_REAL16 / MPI_COMPLEX32 datatypes) because Intel MPI's builtin
 * MPI_SUM/MAX/MIN have no 16-byte-real kernel. Those ops are created at
 * runtime by quad_mpi_init() (quad_mpi.c) once MPI is up; BLACS's
 * blacs_pinfo is the first thing MUMPS calls after MPI_Init, so this is
 * where we register them — mirroring the multifloats bridge's hook.
 * Declared extern here rather than via a header to keep the plain-C
 * BLACS build free of any quad-specific include path. */
extern void quad_mpi_init(void);

#if (INTFACE == C_CALL)
void Cblacs_pinfo(Int *mypnum, Int *nprocs)
#else
F_VOID_FUNC blacs_pinfo_(Int *mypnum, Int *nprocs)
#endif
{
   Int ierr;
   extern Int BI_Iam, BI_Np;
   MpiInt flag, Iam = BI_Iam, Np = BI_Np;
   MpiInt argc=0;
   char **argv=NULL;
   if (BI_COMM_WORLD == NULL)
   {
      MPI_Initialized(&flag);

      if (!flag)
         ierr = MPI_Init(&argc,&argv);

      BI_COMM_WORLD = (Int *) malloc(sizeof(Int));
      *BI_COMM_WORLD = MPI_Comm_c2f(MPI_COMM_WORLD);
   }
   MPI_Comm_size(MPI_COMM_WORLD, &Np);
   MPI_Comm_rank(MPI_COMM_WORLD, &Iam);
   *mypnum = BI_Iam = Iam;
   *nprocs = BI_Np  = Np;

   quad_mpi_init();
}
