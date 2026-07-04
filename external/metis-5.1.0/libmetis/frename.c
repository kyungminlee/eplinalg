/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * Frename.c
 * 
 * THis file contains some renaming routines to deal with different Fortran compilers
 *
 * Started 9/15/97
 * George
 *
 */


#include "metislib.h"

#define FRENAME(name, dargs, cargs, name1, name2, name3, name4)   \
  int name1 dargs { return name cargs; }                          \
  int name2 dargs { return name cargs; }                          \
  int name3 dargs { return name cargs; }                          \
  int name4 dargs { return name cargs; }


FRENAME(
    METIS_MUMPS_PartGraphRecursive, 
    (idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
     idx_t *vsize, idx_t *adjwgt, idx_t *nparts, real_t *tpwgts, 
     real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part),
    (nvtxs, ncon, xadj, adjncy, vwgt, 
     vsize, adjwgt, nparts, tpwgts, 
     ubvec, options, edgecut, part),
    METIS_MUMPS_PARTGRAPHRECURSIVE, 
    metis_mumps_partgraphrecursive, 
    metis_mumps_partgraphrecursive_, 
    metis_mumps_partgraphrecursive__
) 
    

FRENAME(
    METIS_MUMPS_PartGraphKway,
    (idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, 
     idx_t *vsize, idx_t *adjwgt, idx_t *nparts, real_t *tpwgts, 
     real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part),
    (nvtxs, ncon, xadj, adjncy, vwgt, 
     vsize, adjwgt, nparts, tpwgts, 
     ubvec, options, edgecut, part),
    METIS_MUMPS_PARTGRAPHKWAY,
    metis_mumps_partgraphkway,
    metis_mumps_partgraphkway_,
    metis_mumps_partgraphkway__
)

FRENAME(
  METIS_MUMPS_MeshToDual,
  (idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *ncommon, idx_t *numflag, 
   idx_t **r_xadj, idx_t **r_adjncy),
  (ne, nn, eptr, eind, ncommon, numflag, r_xadj, r_adjncy),
  METIS_MUMPS_MESHTODUAL,
  metis_mumps_meshtodual,
  metis_mumps_meshtodual_,
  metis_mumps_meshtodual__
)


FRENAME(
  METIS_MUMPS_MeshToNodal,
  (idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *numflag, idx_t **r_xadj, 
   idx_t **r_adjncy),
  (ne, nn, eptr, eind, numflag, r_xadj, r_adjncy),
  METIS_MUMPS_MESHTONODAL,
  metis_mumps_meshtonodal,
  metis_mumps_meshtonodal_,
  metis_mumps_meshtonodal__
)
  

FRENAME(
  METIS_MUMPS_PartMeshNodal,
  (idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *vwgt, idx_t *vsize, 
   idx_t *nparts, real_t *tpwgts, idx_t *options, idx_t *objval, idx_t *epart, 
   idx_t *npart),
  (ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgts, options, objval, epart, npart),
  METIS_MUMPS_PARTMESHNODAL,
  metis_mumps_partmeshnodal,
  metis_mumps_partmeshnodal_,
  metis_mumps_partmeshnodal__
)


FRENAME(
  METIS_MUMPS_PartMeshDual,
  (idx_t *ne, idx_t *nn, idx_t *eptr, idx_t *eind, idx_t *vwgt, idx_t *vsize, 
   idx_t *ncommon, idx_t *nparts, real_t *tpwgts, idx_t *options, idx_t *objval, 
   idx_t *epart, idx_t *npart),
  (ne, nn, eptr, eind, vwgt, vsize, ncommon, nparts, tpwgts, options, objval, epart, npart),
  METIS_MUMPS_PARTMESHDUAL,
  metis_mumps_partmeshdual,
  metis_mumps_partmeshdual_,
  metis_mumps_partmeshdual__
)


FRENAME(
  METIS_MUMPS_NodeND,
  (idx_t *nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *options, idx_t *perm, 
   idx_t *iperm),
  (nvtxs, xadj, adjncy, vwgt, options, perm, iperm),
  METIS_MUMPS_NODEND,
  metis_mumps_nodend,
  metis_mumps_nodend_,
  metis_mumps_nodend__
)


FRENAME(
  METIS_MUMPS_Free,
  (void *ptr),
  (ptr),
  METIS_MUMPS_FREE,
  metis_mumps_free,
  metis_mumps_free_,
  metis_mumps_free__
)


FRENAME(
  METIS_MUMPS_SetDefaultOptions,
  (idx_t *options),
  (options),
  METIS_MUMPS_SETDEFAULTOPTIONS,
  metis_mumps_setdefaultoptions,
  metis_mumps_setdefaultoptions_,
  metis_mumps_setdefaultoptions__
)
    


