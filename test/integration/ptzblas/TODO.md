# tests/ptzblas — TODO / follow-ups

## Potential migrator follow-ups (NOT in this subtree)

None observed during the kind16 build — every test passed on the
first iteration after the upstream-semantics fixes documented in
`common/ref_quad_aux.f90` were added (CABS1 vs Euclidean magnitude
for the complex abs-matvec family; non-circular row/column shift
semantics for `q[r|c]shft`; A-is-M×N / B-is-N×M shape for `qmmtadd`).
If ctest later reports a bug that looks like a migrator artifact (e.g.
`cmplx16` aliasing analogous to the pzher2k issue noted in
`tests/pblas/TODO.md`), record it here with a reference test program
left in place and an explanatory header comment.
