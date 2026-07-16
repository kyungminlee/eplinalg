# Debugging

Hard-won, verified knowledge. Everything here was diagnosed in this
repo; when a symptom matches, start from the listed cause before
re-deriving it.

## Known pitfalls

- **Intel MPI user-op reduce**: Intel MPI's shm reduce can hand a
  commutative user-op callback a buffer aligned to 8-but-not-16 bytes
  above the short-message cutoff; aligned-SSE callbacks (`__float128`)
  fault (#GP → SIGSEGV with `si_addr=0`). The quad/multifloats reduce
  ops therefore register with `commute=0` under an Intel MPI
  compile-time guard (`commute=1` elsewhere — OpenMPI verified
  unaffected). Don't "optimize" it back without re-testing under
  Intel MPI.
- **MPI init contract**: extended collectives/MUMPS require
  `multifloats_mpi_init()` / `quad_mpi_init()` after `MPI_Init`
  (kind10 needs nothing). Missing init is the usual cause of
  "Intel MPI crashes in a quad reduction".
- **Shared-library links**: every produced `.so` needs
  `-Wl,--no-define-common` (Intel MPI COMMON sentinels fork per-`.so`
  otherwise → PT-Scotch parallel analysis dies, `INFOG(1)=-9980`) and
  `-Wl,-z,now` (lazy PLT resolution corrupts live vector-register FP
  state → flaky wrong double-double limbs; `LD_BIND_NOW=1` is the
  consumer-side rescue for an already-built `.so`). Details and the one
  exception: [../user/mkl-coexistence.md](../user/mkl-coexistence.md).
- **Scotch vs `_FORTIFY_SOURCE`**: fortify false-aborts Scotch's based
  arrays ("buffer overflow detected" at `ICNTL(7)=3`); the vendored
  build compiles scotch/esmumps with `-D_FORTIFY_SOURCE=0`.
- **Never rename shipped symbols**: MKL coexistence is solved by the
  source-level `ep_` privatization of the family-independent engine,
  not by link tricks; if a symbol exists in MKL, MKL must provide it.
- **ctest oversubscription**: high `-j` can deadlock MPI tests forever;
  use `-j2 --timeout 300` ([test.md](test.md)).

## Where to look

- **Upstream-source bugs and workarounds**:
  [../upstream-bugs/](../upstream-bugs/README.md) — methodology and
  per-library catalogue. Workarounds are applied as recipe patches;
  `external/` is never edited.
- **Minimal reproducers**: `repro/` holds standalone reproducers for
  the nastiest past bugs (Intel MPI alignment fault, PLT FP corruption,
  …) — useful as templates when isolating a new one.
- **np≥2-only MUMPS failures**: check the implicit API constraints
  first ([../user/api/mumps.md](../user/api/mumps.md)) — host-only
  sparse RHS and the distributed-solution `INFO(23)` slice are the
  classic false alarms.
- **Convergence regressions**: `diverge`/`converge` output pinpoints
  the transform pass; see [codegen.md](codegen.md).

## Technique notes

- Flakiness that varies run-to-run but is deterministic for a given
  address-space layout points at load-order/relocation effects (the
  lazy-PLT FP corruption looked exactly like this): disabling ASLR
  makes it reproducible while isolating.
- When a crash implicates a migrated routine, diff it against its
  Netlib original and its co-family sibling — a divergence between the
  two halves is usually the bug ([convergence.md](convergence.md)).
