/* scotch_rename_pt_mumps.h — distributed (PT-Scotch) companion to
** scotch_rename_mumps.h.
**
** The vendored scotch_rename_mumps.h is generated from the SEQUENTIAL
** libscotch.a and therefore carries no distributed (dgraph/dorder)
** entries. libptscotch is compiled with the same
** -DSCOTCH_RENAME -DSCOTCH_NAME_SUFFIX=_mumps namespacing as the
** sequential library, so its public distributed symbols also carry the
** _mumps suffix (e.g. SCOTCH_dgraphInit -> SCOTCH_dgraphInit_mumps).
**
** MUMPS's C bridge (mumps_scotch.c) calls exactly one distributed
** Scotch entry point by its bare name — SCOTCH_dgraphInit (its
** MUMPS_DGRAPHINIT shim). Force-include this header (C-only) alongside
** scotch_rename_mumps.h on the -Dptscotch build so that bare call binds
** to the suffixed archive symbol, mirroring how the sequential bare
** SCOTCH_* calls are remapped. Every other distributed entry point MUMPS
** uses is reached through the Fortran SCOTCHFDGRAPH* interface, remapped
** separately via CPP defines in the cmake _scotch_pt_f_renames list.
*/

#ifndef SCOTCH_RENAME_PT_MUMPS_H
#define SCOTCH_RENAME_PT_MUMPS_H

#define SCOTCH_dgraphInit               SCOTCH_dgraphInit_mumps

#endif /* SCOTCH_RENAME_PT_MUMPS_H */
