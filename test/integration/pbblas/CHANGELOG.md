# tests/pbblas — CHANGELOG

Resolved items, reverse-chronological. Open work lives in `TODO.md`.

## 2026-05-06 — `scalapack_tools` recipe

NUMROC / ICEIL / ILCM are now produced by the new
`recipes/scalapack_tools.yaml` (a tiny migrator output that carries
just those three precision-independent integer helpers from
`external/scalapack-2.2.3/TOOLS/`). Both `${LIB_PREFIX}pbblas` and
`${LIB_PREFIX}scalapack` PUBLIC-link `${LIB_PREFIX}scalapack_tools`,
so downstream test executables that link either parent archive pick
the symbols up via the normal dependency graph — no test-side shim,
no duplicate-symbol risk when an executable links both.

Replaces the hand-rolled `common/scalapack_tools.f90` shim (deleted
in this commit). The shim's documented `IDENOM == 0 → 0` /
`A==0 ∨ B==0 → 0` guards in ICEIL / ILCM no longer apply: PBBLAS
never exercises those input ranges, and the new library uses the
upstream-verbatim sources.

## 2026-05-02 — `pbdtrnv` / `pbztrnv` covered

`test_pbdtrnv.f90` / `test_pbztrnv.f90` (baseline XDIST='C', NZ=0,
2×2 grid) plus four families of follow-ups, all passing on every target
(kind10 / kind16 / multifloats):

- `test_pb[dz]trnv_xdistr.f90` — XDIST='R' (row-vector input →
  column-vector output) on the default 2×2 grid.
- `test_pb[dz]trnv_replicated.f90` — IXCOL=-1 and IYROW=-1
  simultaneously; X populated on every process column, Y verified on
  row 0 against the reference and cross-checked byte-equal against row
  `nprow-1`.
- `test_pb[dz]trnv_nzoffset.f90` — NZ > 0 (NB=4, NZ=2, N=12).
  Source-row / source-column local layouts honor the
  `numroc(NN) - NZ` shrink for the offset extended NN=N+NZ space;
  non-source processes hold full NB-sized blocks.
- `test_pb[dz]trnv_lcm.f90` — non-square grid via the new
  `grid_init_shape` helper. `pbdtrnv_lcm` runs 4×1 (LCMQ=4) with
  XDIST='C'; `pbztrnv_lcm` runs 1×4 (LCMP=4) with XDIST='R'. Both
  require a connected 4-rank MPI world: skipped (with a "skipped"
  pseudo-case) when the launcher spawns unconnected worlds. Verified
  on the linux-impi preset with `FI_PROVIDER=tcp` +
  `I_MPI_HYDRA_BOOTSTRAP=fork`.

## 2026-05-02 — `pbdtran` / `pbztran` covered

Baseline `test_pbdtran.f90` / `test_pbztran.f90` cover ADIST='C' with
IACOL=0 and ADIST='R' with IAROW=0. Replicated-A added via
`test_pb[dz]tran_replicated.f90` — both ADIST='C' with IACOL=-1 (every
process column holds A) and ADIST='R' with IAROW=-1 (every process row
holds A). These exercise the separate "all column processors have a
copy" branch in pbdtran.f (line ~352 onward) which uses `PBDTRGET` /
`PBDTRSRT` rather than the source-node-sends pattern.

2×3 / 3×2 grids with `MIN(LCMQ,CEIL(M,NB)) > 1` covered via
`test_pbdtran_lcm6.f90` (2×3, ADIST='C', M=12, N=NB=4 ⇒ LCMQ=2,
CEIL(M,NB)=3, MIN=2) and `test_pbztran_lcm6.f90` (3×2, ADIST='R',
M=NB=4, N=12 ⇒ LCMP=2, CEIL(N,NB)=3, MIN=2; TRANS='C'). These need a
connected 6-rank MPI world; the harness exposes a
`PBBLAS_TEST_NPROC_LCM6=6` knob for filenames matching `*_lcm6` and the
tests self-skip otherwise. Verified on the linux-impi preset for all
three targets — kind16 and kind10 (max_rel_err ≤ 3.3e-20) and
multifloats (max_rel_err ≤ 6.6e-33).

## 2026-04-30 — Multi-target coverage

All three target shims (`target_kind10/`, `target_kind16/`,
`target_multifloats/`) ship in this subtree and build cleanly via
`pbblas_test_target` on the corresponding stagings. Runtime exercise
still requires real MPI under `--preset=linux-impi`.
