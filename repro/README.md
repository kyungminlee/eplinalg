# Upstream LAPACK bug repro — `sgelss.f` workspace query inconsistency

Reproduces and documents the latent inconsistency between `sgelss.f`
and `dgelss.f` in their path-2a workspace estimates, catalogued as
**L10** in `doc/UPSTREAM_BUGS.md`.

## The bug

`sgelss.f` line 324 (LAPACK 3.12.1) computes the optimal LWORK for
its inner `SGELQF` call via the ILAENV block-size formula:

```fortran
MAXWRK = M + M*ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )
```

`dgelss.f` lines 308–310 and 328 query `DGELQF` directly with
`LWORK=-1` and use the returned size:

```fortran
CALL DGELQF( M, N, A, LDA, DUM(1), DUM(1), -1, INFO )
LWORK_DGELQF = INT( DUM(1) )
...
MAXWRK = M + LWORK_DGELQF
```

For the reference Netlib `SGELQF`, the optimal LWORK is exactly
`M*NB`, so the two formulas coincide. The reference build never
crashes from this inconsistency. The bug is latent.

The mismatch becomes consequential if any vendor (Intel MKL,
OpenBLAS, AOCL, …) ships an `SGELQF` whose optimal LWORK includes
overhead beyond `M*NB` — workspace for the triangular factor T,
lookahead buffers, or future algorithmic improvements. In that case
`SGELSS` path 2a would silently under-allocate, and the inner
`SGELQF` call would fail or corrupt memory.

## Files

- `sgelss_workspace.f90` — reproducer. Queries `SGELSS`, `SGELQF`,
  `DGELSS`, `DGELQF` independently, prints the workspace each
  reports, and confirms the ILAENV formula and the direct query
  coincide on Netlib reference.
- `patches/sgelss_workspace.patch` — proposed upstream fix, against
  Netlib LAPACK 3.12.1. Mirrors the D-half pattern: declares
  `LWORK_SGELQF`, adds the direct `SGELQF` query, replaces the
  ILAENV formula with `M + LWORK_SGELQF`. Same structure as
  `dgelss.f`'s analogous block.

## Build and run

```bash
gfortran -O2 sgelss_workspace.f90 -llapack -lblas -o sgelss_workspace
./sgelss_workspace
```

Expected output on a Netlib build:

```
sgelss workspace reproducer:
  M=64  N=512  NRHS=1

=== sgelss.f path 2a SGELQF term ===
  formula (line 324):  M + M*ILAENV(1,SGELQF) = 64 + 2048
                                    total =  2112
  direct  (would-be):  M + LWORK_SGELQF      = 64 + 2048
                                    total =  2112

=== dgelss.f path 2a DGELQF term (already correct) ===
  formula (line 328):  M + LWORK_DGELQF      = 64 + 2048
                                    total =  2112

  SGELSS reported optimal LWORK = 8448
  DGELSS reported optimal LWORK = 8448

=== Verdict ===
  Reference SGELQF: ILAENV formula and direct
  query coincide (both = M*NB). Bug is LATENT
  in the reference build — no visible crash.
```

## Applying the upstream fix

```bash
cd $LAPACK_ROOT
git apply path/to/patches/sgelss_workspace.patch
```

The patch is a 4-line semantic change (one new variable declaration,
one new query call, one line removed for the ILAENV formula
replaced) plus comments / line-numbering adjustments. After applying,
`sgelss.f`'s path-2a workspace estimate matches the direct-query
pattern that the rest of the LAPACK family already uses.
