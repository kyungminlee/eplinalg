# Coexisting with MKL (LP64)

The extended-precision stacks (kind10 `ey`, kind16 `qx`, multifloats `mw`)
ship with a **privatized** copy of the family-independent BLACS / PBLAS /
ScaLAPACK engine: every family-independent Netlib name (302 of them —
`blacs_gridinit`, `descinit`, `numroc`, the `BI_*` internals, …) is emitted
at source level with an `ep_` prefix, in the `libep*_common` archives
(`libepblacs_common`, `libeppblas_common`, `libepscalapack_common`,
`libepptzblas_common`). None of those names exist in MKL, and the extended
archives reference only the `ep_` twins, so linking a full MKL
(`mkl_scalapack_lp64` + `mkl_blacs_*_lp64` + core) into the same executable
can never collapse state with the extended engine — in **either** link
order.

The genuine kind4/kind8 families (s/d/c/z, including the shipped `dz`/`sc`
MUMPS) intentionally keep the pristine Netlib names, so at link time they
bind to MKL's BLACS/ScaLAPACK when MKL is present. That is the supported
configuration: MKL serves double precision, the `ep_` engine serves the
extended families, both in one process.

## Consumer contract

**MUMPS path** (recommended: extended precision via `mmumps_c`, `qmumps_c`,
… or the Fortran equivalents):

- Call the family's MPI setup once, after `MPI_Init` and before the first
  MUMPS call — this is the pre-existing extended-MUMPS contract, unchanged
  by the privatization:
  - multifloats (`m`/`w`): `multifloats_mpi_init()` — registers the custom
    MPI datatypes and reduce ops.
  - kind16 (`q`/`x`): `quad_mpi_init()` — registers the `MPI_QQ`/`MPI_XX`
    reduce families (Intel MPI has no 16-byte-real reduce kernel).
  - kind10 (`e`/`y`): nothing — long-double reductions use builtin ops.
- Nothing else. The extended MUMPS internally creates its grids through the
  `ep_` engine.

**Direct extended ScaLAPACK/PBLAS use**: create grids and descriptors
through the privatized entry points — `ep_blacs_pinfo` / `ep_blacs_get` /
`ep_blacs_gridinit` / `ep_blacs_gridinfo` / `ep_blacs_gridexit`,
`ep_descinit`, `ep_numroc`, … (C bindings: `ep_Cblacs_*`). The
family-prefixed computational routines keep their usual names (`pmgemm`,
`pqgetrf`, …). Never pass an MKL BLACS context handle into an extended
routine, or an `ep_` context into MKL — the two engines are separate
universes with private state; a context is only meaningful inside the
engine that created it.

## Link line

- Both orders (eplinalg before or after the MKL group) are supported and
  verified; group the eplinalg archives with
  `-Wl,--start-group … -Wl,--end-group`.
- Mixed dz+extended MUMPS consumers need the ordering closure archives
  `libesmumps_mumps.a` and `libmetis_mumps.a`.
- Toolchains defaulting to `--as-needed` (most Linux distros) drop MKL
  shared libraries placed before the objects that reference them; wrap the
  MKL group in `-Wl,--push-state,--no-as-needed … -Wl,--pop-state`.

## What is *not* supported

Calling the generic BLACS names (`blacs_gridinit`, `Cblacs_get`, …) and
handing the resulting context to extended-family transfer or compute
routines. With MKL linked, the generic names are MKL's; before the
privatization this "worked" only when fragile link-order accidents let the
bundled engine win — the same accident that corrupted MKL's own state
(`BLACBUFF` layout mismatch) and crashed double-precision solves. Use the
`ep_` names for the extended universe.
