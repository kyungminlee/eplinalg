# LAPACK differential precision tests

Differential precision tests for the migrated LAPACK libraries
(`qxlapack` / `eylapack` / `mwlapack`), following the scheme described in
`tests/README.md`. Each test runs a routine through the migrated
extended-precision implementation and through a REAL(KIND=16) reference
(`reflapack_quad`), then compares the outputs in KIND=16 and writes a
per-test JSON report.

Routines covered (15):

| group          | routines                                                |
|----------------|---------------------------------------------------------|
| linear_solve   | dgesv, dgetrf, dgetrs, dpotrf, dpotrs, zgesv            |
| factorization  | dgeqrf, dorgqr, zgeqrf                                  |
| eigenvalue     | dsyev, zheev                                            |
| svd            | dgesvd                                                  |
| auxiliary      | dlange, dlacpy, dlaset                                  |

For `dsyev`, `zheev`, and `dgesvd` the eigenvector / singular-vector
columns have a precision-sensitive sign convention (an internal LAPACK
test on a near-zero value can flip across precisions). The tests
therefore compare only the canonical unique outputs (eigenvalues W,
singular values S) directly, and validate the vectors through the
sign-immune invariants `||A·Z − Z·diag(W)||/||A||` and
`||A − U·diag(S)·Vᵀ||/||A||`.

