# Compact summary: OMP=1 parallel-blas-vs-migrated speedup ranges

Per (precision, routine), parallel-blas-overlay vs migrated-Fortran speedup distribution across all measured (combo × size) cases.
Speedup = `t_migrated / t_parallel-blas` (>1 = parallel-blas wins).

| prec | routine | combos × sizes | min spd | median | max spd | <1.0× cases |
|------|---------|----------------|---------|--------|---------|-------------|
| kind10 | egemm | 36 | 0.82× | 1.42× | 2.99× | 6/36 |
| kind10 | esymm | 4 | 0.90× | 1.12× | 1.27× | 1/4 |
| kind10 | esyrk | 4 | 1.00× | 1.65× | 2.29× | 0/4 |
| kind10 | etrmm | 4 | 0.90× | 1.79× | 2.47× | 1/4 |
| kind10 | etrsm | 4 | 1.10× | 1.71× | 2.30× | 0/4 |
| kind10 | ygemm | 36 | 0.59× | 0.76× | 1.10× | 32/36 |
| kind10 | yhemm | 4 | 0.83× | 0.91× | 1.13× | 3/4 |
| kind10 | yherk | 4 | 1.02× | 1.09× | 1.19× | 0/4 |
| kind10 | ysymm | 4 | 1.02× | 1.17× | 1.19× | 0/4 |
| kind10 | ysyrk | 4 | 1.00× | 1.18× | 1.23× | 0/4 |
| kind10 | ytrmm | 4 | 1.03× | 1.10× | 1.12× | 0/4 |
| kind10 | ytrsm | 4 | 1.10× | 1.17× | 1.18× | 0/4 |
| kind16 | qgemm | 18 | 0.97× | 1.00× | 1.07× | 8/18 |
| kind16 | qsymm | 2 | 1.00× | 1.00× | 1.00× | 0/2 |
| kind16 | qsyrk | 2 | 1.02× | 1.06× | 1.11× | 0/2 |
| kind16 | qtrmm | 2 | 1.00× | 1.06× | 1.12× | 0/2 |
| kind16 | qtrsm | 2 | 0.99× | 1.04× | 1.10× | 1/2 |
| kind16 | xgemm | 18 | 0.99× | 1.00× | 1.04× | 2/18 |
| kind16 | xhemm | 2 | 1.00× | 1.01× | 1.02× | 0/2 |
| kind16 | xherk | 2 | 0.97× | 0.98× | 1.00× | 1/2 |
| kind16 | xsymm | 2 | 1.00× | 1.00× | 1.00× | 0/2 |
| kind16 | xsyrk | 2 | 0.98× | 0.98× | 0.99× | 2/2 |
| kind16 | xtrmm | 2 | 0.99× | 0.99× | 1.00× | 1/2 |
| kind16 | xtrsm | 2 | 0.99× | 0.99× | 1.00× | 1/2 |
| multifloats | mgemm | 36 | 20.19× | 29.26× | 30.58× | 0/36 |
| multifloats | msymm | 4 | 4.21× | 9.02× | 16.37× | 0/4 |
| multifloats | msyrk | 4 | 4.03× | 8.79× | 16.15× | 0/4 |
| multifloats | mtrmm | 4 | 4.13× | 9.19× | 16.33× | 0/4 |
| multifloats | mtrsm | 4 | 3.84× | 24.09× | 28.12× | 0/4 |
| multifloats | wgemm | 36 | 13.21× | 18.46× | 22.21× | 0/36 |
| multifloats | whemm | 4 | 2.12× | 4.92× | 9.33× | 0/4 |
| multifloats | wherk | 4 | 2.43× | 5.63× | 9.80× | 0/4 |
| multifloats | wsymm | 4 | 2.14× | 5.08× | 9.39× | 0/4 |
| multifloats | wsyrk | 4 | 2.23× | 5.00× | 9.43× | 0/4 |
| multifloats | wtrmm | 4 | 2.31× | 5.19× | 9.52× | 0/4 |
| multifloats | wtrsm | 4 | 2.24× | 15.57× | 17.44× | 0/4 |

## Worst regressions (parallel-blas slower than migrated)

| speedup | prec | routine | combo | size |
|---------|------|---------|-------|------|
| 0.59× | kind10 | ygemm | TN | 128 |
| 0.60× | kind10 | ygemm | TN | 64 |
| 0.61× | kind10 | ygemm | TN | 256 |
| 0.61× | kind10 | ygemm | TT | 64 |
| 0.62× | kind10 | ygemm | TN | 512 |
| 0.62× | kind10 | ygemm | TT | 128 |
| 0.65× | kind10 | ygemm | CN | 64 |
| 0.65× | kind10 | ygemm | CT | 128 |
| 0.65× | kind10 | ygemm | TC | 64 |
| 0.65× | kind10 | ygemm | TC | 128 |
| 0.66× | kind10 | ygemm | CN | 256 |
| 0.67× | kind10 | ygemm | CN | 128 |
| 0.67× | kind10 | ygemm | CN | 512 |
| 0.67× | kind10 | ygemm | CT | 64 |
| 0.69× | kind10 | ygemm | CC | 128 |
| 0.70× | kind10 | ygemm | CC | 64 |
| 0.73× | kind10 | ygemm | NN | 64 |
| 0.75× | kind10 | ygemm | TT | 256 |
| 0.77× | kind10 | ygemm | CC | 256 |
| 0.78× | kind10 | ygemm | CT | 256 |
| 0.81× | kind10 | ygemm | TC | 256 |
| 0.82× | kind10 | egemm | NN | 64 |
| 0.82× | kind10 | ygemm | NN | 128 |
| 0.82× | kind10 | ygemm | NT | 256 |
| 0.83× | kind10 | ygemm | NC | 64 |
| 0.83× | kind10 | ygemm | NC | 128 |
| 0.83× | kind10 | ygemm | NC | 256 |
| 0.83× | kind10 | ygemm | NN | 256 |
| 0.83× | kind10 | ygemm | NN | 512 |
| 0.83× | kind10 | ygemm | NT | 64 |