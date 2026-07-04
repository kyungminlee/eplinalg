# Implicit MUMPS API constraints (np ≥ 2)

Some MUMPS input/output conventions are **invisible at `np = 1`** and only
manifest once the solve runs on two or more MPI ranks. A single-rank run is
simultaneously the host and the only worker, so any host-vs-slave asymmetry in
the API collapses and a mis-driven test still passes. The moment a slave rank
exists, the same code either error-stops with `INFO(1) = -22` or reads past the
end of an array and segfaults.

These are **not** migration defects and **not** upstream bugs — they are
documented (if easily-missed) requirements of the MUMPS API that the migrated
`?mumps` archives honour exactly as Netlib does. They are recorded here because
the `tests/mumps` harness tripped over each one when it moved from `mpiexec -n 1`
to a genuine multi-rank default (see `tests/mumps/TODO.md`, and the drivers
`fortran/test_?mumps_icntl_io.f90`).

The reference line numbers below are into the vendored source
`external/MUMPS_5.8.2/src/?sol_driver.F` (the `?` is `d`/`z`/… per arithmetic).

---

## 1. Centralized sparse RHS (`ICNTL(20)=1`) is a host-only input

**Symptom (np ≥ 2):** `INFOG(1) = -22`, `INFO(2) = 7` (the code for the `RHS`
family of arrays).

**Rule:** when a centralized sparse right-hand side is active
(`ICNTL(20)=1`, i.e. `KEEP(248)==1` internally), the arrays

- `id%RHS`
- `id%RHS_SPARSE`
- `id%IRHS_SPARSE`
- `id%IRHS_PTR`

must be **associated only on the host rank (`myid == 0`)** and must be left
**null on every slave rank**. The solve driver explicitly rejects any non-master
rank that has them associated (`?sol_driver.F`, roughly lines 1122–1147). The
host-side dense-RHS validator (`?MUMPS_CHECK_DENSE_RHS`, in `?mumps_driver.F`)
independently requires them present and correctly sized *on the host*, so the
two checks together pin down exactly one legal layout.

The coefficient matrix (`IRN` / `JCN` / `A`) is unaffected — centralized
assembly tolerates a full copy on every rank, so it may stay allocated
everywhere as in the dense cases.

**Correct pattern:**

```fortran
id%ICNTL(20) = 1
id%N   = n
id%NNZ = int(nz, kind=8)
allocate(id%IRN(nz));  id%IRN = irn      ! matrix: on every rank
allocate(id%JCN(nz));  id%JCN = jcn
allocate(id%A(nz));    id%A   = ...
id%NRHS = 1
if (myid == 0) then                      ! sparse RHS + solution: HOST ONLY
    id%NZ_RHS = n
    allocate(id%RHS_SPARSE(n));  id%RHS_SPARSE = ...
    allocate(id%IRHS_SPARSE(n))
    allocate(id%IRHS_PTR(2))
    id%IRHS_PTR = [1, n + 1]
    allocate(id%RHS(n))                  ! solution lands here, centralized
end if
id%JOB = 6
call ?mumps(id)
if (myid == 0) then
    ! read solution from id%RHS ...
    deallocate(id%RHS_SPARSE, id%IRHS_SPARSE, id%IRHS_PTR)
end if
```

The same rule applies to a centralized **dense** RHS (`ICNTL(20)=0`): `id%RHS`
is a host-only array. It is easy to miss there too, but the dense path is more
commonly driven host-only by habit, so the sparse path is where tests trip.

---

## 2. Distributed solution (`ICNTL(21)=1`) slice size is `INFO(23)`, not `NSOL_loc`

**Symptom (np ≥ 2):** SIGSEGV in the post-solve result loop (reading past the
end of `SOL_loc` / `ISOL_loc`).

**Rule:** with `ICNTL(21)=1` the solution is returned **distributed**: each rank
receives its own contiguous slice, and the number of local components on this
rank is `id%INFO(23)` (the count of pivots eliminated on this processor, per the
User's Guide). The global row index of local component `k` is `id%ISOL_loc(k)`;
the union of all ranks' slices is the full solution.

The per-rank loop bound is therefore `INFO(23)` — **not**:

- `id%NSOL_loc` — this field is **never written by the solve**; reading it as a
  loop bound uses an uninitialized / stale value and walks off the array.
- `id%N` — that is the global size, valid only when `np = 1` (where the lone
  rank owns all `N`).

The user must pre-allocate `SOL_loc` / `ISOL_loc` to at least the largest slice
any rank will receive, and set `id%LSOL_loc` to that capacity **before** the
solve. Using `N` as the capacity is a safe upper bound.

**Correct pattern:**

```fortran
id%ICNTL(21) = 1
id%LSOL_loc  = n                 ! capacity (upper bound on any rank's slice)
allocate(id%SOL_loc(n))
allocate(id%ISOL_loc(n))
id%JOB = 6
call ?mumps(id)
do k = 1, id%INFO(23)            ! <-- INFO(23), the documented local slice size
    gi = id%ISOL_loc(k)          ! global row index of this local component
    ! use id%SOL_loc(k) as x(gi) ...
end do
```

Because the slice is genuinely partitioned at `np ≥ 2` (e.g. at `np = 2` one
rank may own 0 components and the other all `N`), any correctness check that
compares against a global reference must reduce across ranks — see the scalar
`MPI_MAX` reduction of the local max-error in `test_?mumps_icntl_io.f90`.

---

## Operational note: serialize the ctest suite

Unrelated to the API but in the same spirit of "only bites at `np ≥ 2`": the
mumps ctest suite must be run **`-j1`**. Concurrent `mpiexec -n N` launches race
on the shared OpenMPI session directory and fail spuriously. The suite is
registered with a resource lock for this reason; when driving `ctest` by hand,
pass `-j1` (and export `OMPI_MCA_rmaps_base_oversubscribe=1`
`OMPI_MCA_mpi_yield_when_idle=1` when the rank count exceeds the core budget).
