# Full OMP=1 parallel-blas overlay vs migrated Fortran benchmark

Per-routine speedup at OMP=1. Run with iters=3, warmup=1, sizes=64/128/256/512.

Scope: parallel-blas overlay only. (This report predates the epopenblas overlay; see `reports/cmp5/` for epopenblas vs parallel-blas comparisons.)

Speedup column: `t_migrated / t_parallel-blas` (>1 = parallel-blas wins).

## gemm

| prec | routine | trans | size | parallel-blas GF | migrated GF | speedup |
|------|---------|-------|------|-----------------|-------------|---------|
| kind10 | egemm | CC | 64 | 2.1308 | 2.0919 | 1.02× |
| kind10 | egemm | CC | 128 | 2.1766 | 1.9899 | 1.09× |
| kind10 | egemm | CC | 256 | 2.2158 | 1.3226 | 1.68× |
| kind10 | egemm | CC | 512 | 2.2336 | 0.7521 | 2.97× |
| kind10 | egemm | CN | 64 | 2.0964 | 2.3391 | 0.90× |
| kind10 | egemm | CN | 128 | 2.2385 | 2.2430 | 1.00× |
| kind10 | egemm | CN | 256 | 2.2695 | 2.3095 | 0.98× |
| kind10 | egemm | CN | 512 | 2.1532 | 2.3674 | 0.91× |
| kind10 | egemm | CT | 64 | 2.2130 | 1.9307 | 1.15× |
| kind10 | egemm | CT | 128 | 2.1788 | 1.9938 | 1.09× |
| kind10 | egemm | CT | 256 | 2.2137 | 1.3204 | 1.68× |
| kind10 | egemm | CT | 512 | 2.2346 | 0.7484 | 2.99× |
| kind10 | egemm | NC | 64 | 2.0502 | 1.2997 | 1.58× |
| kind10 | egemm | NC | 128 | 2.1400 | 1.3351 | 1.60× |
| kind10 | egemm | NC | 256 | 2.1856 | 1.3803 | 1.58× |
| kind10 | egemm | NC | 512 | 2.2150 | 1.3535 | 1.64× |
| kind10 | egemm | NN | 64 | 1.9420 | 2.3693 | 0.82× |
| kind10 | egemm | NN | 128 | 1.8175 | 1.3747 | 1.32× |
| kind10 | egemm | NN | 256 | 2.1840 | 1.3958 | 1.56× |
| kind10 | egemm | NN | 512 | 2.1988 | 1.4106 | 1.56× |
| kind10 | egemm | NT | 64 | 2.0792 | 1.3779 | 1.51× |
| kind10 | egemm | NT | 128 | 2.1418 | 1.3558 | 1.58× |
| kind10 | egemm | NT | 256 | 2.1969 | 1.3762 | 1.60× |
| kind10 | egemm | NT | 512 | 2.1995 | 1.3548 | 1.62× |
| kind10 | egemm | TC | 64 | 2.0255 | 1.7920 | 1.13× |
| kind10 | egemm | TC | 128 | 2.1717 | 1.9907 | 1.09× |
| kind10 | egemm | TC | 256 | 2.2201 | 1.2982 | 1.71× |
| kind10 | egemm | TC | 512 | 2.1739 | 0.7327 | 2.97× |
| kind10 | egemm | TN | 64 | 2.1142 | 2.3272 | 0.91× |
| kind10 | egemm | TN | 128 | 2.2345 | 2.2066 | 1.01× |
| kind10 | egemm | TN | 256 | 2.2302 | 2.2313 | 1.00× |
| kind10 | egemm | TN | 512 | 2.2424 | 2.3600 | 0.95× |
| kind10 | egemm | TT | 64 | 2.1519 | 2.0920 | 1.03× |
| kind10 | egemm | TT | 128 | 2.1660 | 1.9861 | 1.09× |
| kind10 | egemm | TT | 256 | 2.2112 | 1.2529 | 1.76× |
| kind10 | egemm | TT | 512 | 2.1224 | 0.7292 | 2.91× |
| kind10 | ygemm | CC | 64 | 1.5826 | 2.2587 | 0.70× |
| kind10 | ygemm | CC | 128 | 1.5853 | 2.2913 | 0.69× |
| kind10 | ygemm | CC | 256 | 1.6236 | 2.1048 | 0.77× |
| kind10 | ygemm | CC | 512 | 1.6320 | 1.4935 | 1.09× |
| kind10 | ygemm | CN | 64 | 1.5871 | 2.4371 | 0.65× |
| kind10 | ygemm | CN | 128 | 1.6188 | 2.4195 | 0.67× |
| kind10 | ygemm | CN | 256 | 1.6252 | 2.4780 | 0.66× |
| kind10 | ygemm | CN | 512 | 1.6324 | 2.4334 | 0.67× |
| kind10 | ygemm | CT | 64 | 1.5867 | 2.3804 | 0.67× |
| kind10 | ygemm | CT | 128 | 1.5918 | 2.4360 | 0.65× |
| kind10 | ygemm | CT | 256 | 1.6178 | 2.0736 | 0.78× |
| kind10 | ygemm | CT | 512 | 1.6215 | 1.4773 | 1.10× |
| kind10 | ygemm | NC | 64 | 1.5873 | 1.9069 | 0.83× |
| kind10 | ygemm | NC | 128 | 1.6223 | 1.9499 | 0.83× |
| kind10 | ygemm | NC | 256 | 1.6175 | 1.9374 | 0.83× |
| kind10 | ygemm | NC | 512 | 1.6317 | 1.9171 | 0.85× |
| kind10 | ygemm | NN | 64 | 1.3480 | 1.8482 | 0.73× |
| kind10 | ygemm | NN | 128 | 1.6110 | 1.9543 | 0.82× |
| kind10 | ygemm | NN | 256 | 1.6267 | 1.9693 | 0.83× |
| kind10 | ygemm | NN | 512 | 1.6318 | 1.9578 | 0.83× |
| kind10 | ygemm | NT | 64 | 1.5847 | 1.9095 | 0.83× |
| kind10 | ygemm | NT | 128 | 1.6118 | 1.9203 | 0.84× |
| kind10 | ygemm | NT | 256 | 1.6274 | 1.9741 | 0.82× |
| kind10 | ygemm | NT | 512 | 1.6314 | 1.9237 | 0.85× |
| kind10 | ygemm | TC | 64 | 1.5831 | 2.4421 | 0.65× |
| kind10 | ygemm | TC | 128 | 1.6005 | 2.4647 | 0.65× |
| kind10 | ygemm | TC | 256 | 1.6214 | 2.0063 | 0.81× |
| kind10 | ygemm | TC | 512 | 1.6310 | 1.5319 | 1.06× |
| kind10 | ygemm | TN | 64 | 1.5861 | 2.6566 | 0.60× |
| kind10 | ygemm | TN | 128 | 1.6062 | 2.7009 | 0.59× |
| kind10 | ygemm | TN | 256 | 1.6129 | 2.6518 | 0.61× |
| kind10 | ygemm | TN | 512 | 1.6213 | 2.5992 | 0.62× |
| kind10 | ygemm | TT | 64 | 1.5695 | 2.5703 | 0.61× |
| kind10 | ygemm | TT | 128 | 1.6224 | 2.6363 | 0.62× |
| kind10 | ygemm | TT | 256 | 1.6153 | 2.1679 | 0.75× |
| kind10 | ygemm | TT | 512 | 1.6297 | 1.6026 | 1.02× |
| kind16 | qgemm | CC | 128 | 0.0593 | 0.0602 | 0.98× |
| kind16 | qgemm | CC | 256 | 0.0599 | 0.0615 | 0.97× |
| kind16 | qgemm | CN | 128 | 0.0601 | 0.0603 | 1.00× |
| kind16 | qgemm | CN | 256 | 0.0617 | 0.0605 | 1.02× |
| kind16 | qgemm | CT | 128 | 0.0600 | 0.0601 | 1.00× |
| kind16 | qgemm | CT | 256 | 0.0604 | 0.0615 | 0.98× |
| kind16 | qgemm | NC | 128 | 0.0562 | 0.0547 | 1.03× |
| kind16 | qgemm | NC | 256 | 0.0579 | 0.0565 | 1.03× |
| kind16 | qgemm | NN | 128 | 0.0555 | 0.0520 | 1.07× |
| kind16 | qgemm | NN | 256 | 0.0576 | 0.0547 | 1.05× |
| kind16 | qgemm | NT | 128 | 0.0561 | 0.0546 | 1.03× |
| kind16 | qgemm | NT | 256 | 0.0579 | 0.0566 | 1.02× |
| kind16 | qgemm | TC | 128 | 0.0599 | 0.0602 | 0.99× |
| kind16 | qgemm | TC | 256 | 0.0600 | 0.0612 | 0.98× |
| kind16 | qgemm | TN | 128 | 0.0592 | 0.0603 | 0.98× |
| kind16 | qgemm | TN | 256 | 0.0618 | 0.0619 | 1.00× |
| kind16 | qgemm | TT | 128 | 0.0600 | 0.0604 | 0.99× |
| kind16 | qgemm | TT | 256 | 0.0605 | 0.0615 | 0.98× |
| kind16 | xgemm | CC | 128 | 0.0574 | 0.0575 | 1.00× |
| kind16 | xgemm | CC | 256 | 0.0587 | 0.0589 | 1.00× |
| kind16 | xgemm | CN | 128 | 0.0579 | 0.0579 | 1.00× |
| kind16 | xgemm | CN | 256 | 0.0590 | 0.0591 | 1.00× |
| kind16 | xgemm | CT | 128 | 0.0577 | 0.0572 | 1.01× |
| kind16 | xgemm | CT | 256 | 0.0589 | 0.0590 | 1.00× |
| kind16 | xgemm | NC | 128 | 0.0512 | 0.0515 | 0.99× |
| kind16 | xgemm | NC | 256 | 0.0538 | 0.0535 | 1.01× |
| kind16 | xgemm | NN | 128 | 0.0498 | 0.0478 | 1.04× |
| kind16 | xgemm | NN | 256 | 0.0534 | 0.0534 | 1.00× |
| kind16 | xgemm | NT | 128 | 0.0522 | 0.0520 | 1.00× |
| kind16 | xgemm | NT | 256 | 0.0539 | 0.0534 | 1.01× |
| kind16 | xgemm | TC | 128 | 0.0576 | 0.0575 | 1.00× |
| kind16 | xgemm | TC | 256 | 0.0587 | 0.0588 | 1.00× |
| kind16 | xgemm | TN | 128 | 0.0579 | 0.0579 | 1.00× |
| kind16 | xgemm | TN | 256 | 0.0590 | 0.0591 | 1.00× |
| kind16 | xgemm | TT | 128 | 0.0573 | 0.0576 | 0.99× |
| kind16 | xgemm | TT | 256 | 0.0589 | 0.0590 | 1.00× |
| multifloats | mgemm | CC | 64 | 2.4905 | 0.1197 | 20.81× |
| multifloats | mgemm | CC | 128 | 3.3797 | 0.1213 | 27.87× |
| multifloats | mgemm | CC | 256 | 3.6852 | 0.1213 | 30.37× |
| multifloats | mgemm | CC | 512 | 3.7311 | 0.1226 | 30.43× |
| multifloats | mgemm | CN | 64 | 2.7733 | 0.1202 | 23.07× |
| multifloats | mgemm | CN | 128 | 3.3696 | 0.1218 | 27.66× |
| multifloats | mgemm | CN | 256 | 3.6729 | 0.1226 | 29.96× |
| multifloats | mgemm | CN | 512 | 3.7166 | 0.1231 | 30.19× |
| multifloats | mgemm | CT | 64 | 2.4462 | 0.1201 | 20.37× |
| multifloats | mgemm | CT | 128 | 3.4017 | 0.1217 | 27.95× |
| multifloats | mgemm | CT | 256 | 3.6888 | 0.1222 | 30.20× |
| multifloats | mgemm | CT | 512 | 3.7303 | 0.1227 | 30.41× |
| multifloats | mgemm | NC | 64 | 2.8515 | 0.1215 | 23.47× |
| multifloats | mgemm | NC | 128 | 3.6016 | 0.1226 | 29.37× |
| multifloats | mgemm | NC | 256 | 3.7488 | 0.1231 | 30.45× |
| multifloats | mgemm | NC | 512 | 3.7714 | 0.1233 | 30.58× |
| multifloats | mgemm | NN | 64 | 2.5164 | 0.1209 | 20.82× |
| multifloats | mgemm | NN | 128 | 3.2587 | 0.1218 | 26.76× |
| multifloats | mgemm | NN | 256 | 3.5779 | 0.1228 | 29.14× |
| multifloats | mgemm | NN | 512 | 3.7276 | 0.1229 | 30.33× |
| multifloats | mgemm | NT | 64 | 2.4523 | 0.1215 | 20.19× |
| multifloats | mgemm | NT | 128 | 3.4922 | 0.1225 | 28.51× |
| multifloats | mgemm | NT | 256 | 3.7066 | 0.1227 | 30.21× |
| multifloats | mgemm | NT | 512 | 3.7721 | 0.1234 | 30.58× |
| multifloats | mgemm | TC | 64 | 2.7860 | 0.1201 | 23.19× |
| multifloats | mgemm | TC | 128 | 3.3997 | 0.1216 | 27.97× |
| multifloats | mgemm | TC | 256 | 3.6760 | 0.1221 | 30.10× |
| multifloats | mgemm | TC | 512 | 3.7285 | 0.1227 | 30.38× |
| multifloats | mgemm | TN | 64 | 2.8300 | 0.1202 | 23.54× |
| multifloats | mgemm | TN | 128 | 3.4383 | 0.1215 | 28.30× |
| multifloats | mgemm | TN | 256 | 3.6848 | 0.1226 | 30.05× |
| multifloats | mgemm | TN | 512 | 3.7178 | 0.1227 | 30.31× |
| multifloats | mgemm | TT | 64 | 2.8384 | 0.1201 | 23.63× |
| multifloats | mgemm | TT | 128 | 3.4275 | 0.1216 | 28.18× |
| multifloats | mgemm | TT | 256 | 3.6754 | 0.1217 | 30.20× |
| multifloats | mgemm | TT | 512 | 3.7295 | 0.1227 | 30.41× |
| multifloats | wgemm | CC | 64 | 3.1261 | 0.1815 | 17.22× |
| multifloats | wgemm | CC | 128 | 3.6915 | 0.1791 | 20.61× |
| multifloats | wgemm | CC | 256 | 3.8331 | 0.1726 | 22.21× |
| multifloats | wgemm | CC | 512 | 2.7766 | 0.1563 | 17.77× |
| multifloats | wgemm | CN | 64 | 3.1612 | 0.1878 | 16.83× |
| multifloats | wgemm | CN | 128 | 3.6936 | 0.1895 | 19.49× |
| multifloats | wgemm | CN | 256 | 3.8347 | 0.1900 | 20.18× |
| multifloats | wgemm | CN | 512 | 3.8865 | 0.1895 | 20.51× |
| multifloats | wgemm | CT | 64 | 3.0519 | 0.1874 | 16.29× |
| multifloats | wgemm | CT | 128 | 3.7046 | 0.1906 | 19.44× |
| multifloats | wgemm | CT | 256 | 3.8395 | 0.1896 | 20.25× |
| multifloats | wgemm | CT | 512 | 3.8916 | 0.1841 | 21.14× |
| multifloats | wgemm | NC | 64 | 3.3094 | 0.2097 | 15.78× |
| multifloats | wgemm | NC | 128 | 3.7446 | 0.2100 | 17.83× |
| multifloats | wgemm | NC | 256 | 3.8547 | 0.2112 | 18.25× |
| multifloats | wgemm | NC | 512 | 2.7902 | 0.2112 | 13.21× |
| multifloats | wgemm | NN | 64 | 3.1672 | 0.2124 | 14.91× |
| multifloats | wgemm | NN | 128 | 3.5596 | 0.1921 | 18.53× |
| multifloats | wgemm | NN | 256 | 3.7455 | 0.2132 | 17.56× |
| multifloats | wgemm | NN | 512 | 3.8769 | 0.2152 | 18.01× |
| multifloats | wgemm | NT | 64 | 3.1316 | 0.2103 | 14.89× |
| multifloats | wgemm | NT | 128 | 3.7564 | 0.2110 | 17.80× |
| multifloats | wgemm | NT | 256 | 3.8584 | 0.2117 | 18.23× |
| multifloats | wgemm | NT | 512 | 3.8976 | 0.2120 | 18.39× |
| multifloats | wgemm | TC | 64 | 3.2070 | 0.1889 | 16.98× |
| multifloats | wgemm | TC | 128 | 3.7411 | 0.1896 | 19.73× |
| multifloats | wgemm | TC | 256 | 3.8465 | 0.1893 | 20.32× |
| multifloats | wgemm | TC | 512 | 3.8999 | 0.1859 | 20.98× |
| multifloats | wgemm | TN | 64 | 3.1334 | 0.1991 | 15.74× |
| multifloats | wgemm | TN | 128 | 3.7120 | 0.2000 | 18.56× |
| multifloats | wgemm | TN | 256 | 3.8434 | 0.2008 | 19.14× |
| multifloats | wgemm | TN | 512 | 3.8908 | 0.2017 | 19.29× |
| multifloats | wgemm | TT | 64 | 3.1651 | 0.1987 | 15.93× |
| multifloats | wgemm | TT | 128 | 3.7209 | 0.1990 | 18.70× |
| multifloats | wgemm | TT | 256 | 3.8356 | 0.1999 | 19.18× |
| multifloats | wgemm | TT | 512 | 3.8978 | 0.1954 | 19.94× |

## trsm

| prec | routine | size | parallel-blas GF | migrated GF | speedup |
|------|---------|------|-----------------|-------------|---------|
| kind10 | etrsm | 64 | 1.6885 | 1.5408 | 1.10× |
| kind10 | etrsm | 128 | 1.0248 | 0.7097 | 1.44× |
| kind10 | etrsm | 256 | 1.4356 | 0.7296 | 1.97× |
| kind10 | etrsm | 512 | 1.6927 | 0.7373 | 2.30× |
| kind10 | ytrsm | 64 | 1.4194 | 1.2881 | 1.10× |
| kind10 | ytrsm | 128 | 1.5185 | 1.3227 | 1.15× |
| kind10 | ytrsm | 256 | 1.5573 | 1.3229 | 1.18× |
| kind10 | ytrsm | 512 | 1.5808 | 1.3402 | 1.18× |
| kind16 | qtrsm | 128 | 0.0607 | 0.0552 | 1.10× |
| kind16 | qtrsm | 256 | 0.0612 | 0.0620 | 0.99× |
| kind16 | xtrsm | 128 | 0.0539 | 0.0539 | 1.00× |
| kind16 | xtrsm | 256 | 0.0553 | 0.0557 | 0.99× |
| multifloats | mtrsm | 64 | 0.4571 | 0.1192 | 3.84× |
| multifloats | mtrsm | 128 | 2.6670 | 0.1219 | 21.89× |
| multifloats | mtrsm | 256 | 3.0708 | 0.1168 | 26.29× |
| multifloats | mtrsm | 512 | 3.4074 | 0.1212 | 28.12× |
| multifloats | wtrsm | 64 | 0.4441 | 0.1979 | 2.24× |
| multifloats | wtrsm | 128 | 3.0646 | 0.2049 | 14.95× |
| multifloats | wtrsm | 256 | 3.3623 | 0.2076 | 16.20× |
| multifloats | wtrsm | 512 | 3.6886 | 0.2115 | 17.44× |

## trmm

| prec | routine | size | parallel-blas GF | migrated GF | speedup |
|------|---------|------|-----------------|-------------|---------|
| kind10 | etrmm | 64 | 1.5489 | 1.7251 | 0.90× |
| kind10 | etrmm | 128 | 1.0879 | 0.7100 | 1.53× |
| kind10 | etrmm | 256 | 1.4857 | 0.7294 | 2.04× |
| kind10 | etrmm | 512 | 1.8157 | 0.7353 | 2.47× |
| kind10 | ytrmm | 64 | 1.7933 | 1.7484 | 1.03× |
| kind10 | ytrmm | 128 | 1.4984 | 1.3892 | 1.08× |
| kind10 | ytrmm | 256 | 1.5448 | 1.3933 | 1.11× |
| kind10 | ytrmm | 512 | 1.5793 | 1.4075 | 1.12× |
| kind16 | qtrmm | 128 | 0.0618 | 0.0617 | 1.00× |
| kind16 | qtrmm | 256 | 0.0617 | 0.0553 | 1.12× |
| kind16 | xtrmm | 128 | 0.0534 | 0.0541 | 0.99× |
| kind16 | xtrmm | 256 | 0.0554 | 0.0556 | 1.00× |
| multifloats | mtrmm | 64 | 0.5000 | 0.1212 | 4.13× |
| multifloats | mtrmm | 128 | 0.8563 | 0.1219 | 7.03× |
| multifloats | mtrmm | 256 | 1.3905 | 0.1225 | 11.35× |
| multifloats | mtrmm | 512 | 1.9949 | 0.1222 | 16.33× |
| multifloats | wtrmm | 64 | 0.4811 | 0.2084 | 2.31× |
| multifloats | wtrmm | 128 | 0.8422 | 0.2115 | 3.98× |
| multifloats | wtrmm | 256 | 1.3560 | 0.2120 | 6.40× |
| multifloats | wtrmm | 512 | 2.0488 | 0.2153 | 9.52× |

## syrk

| prec | routine | size | parallel-blas GF | migrated GF | speedup |
|------|---------|------|-----------------|-------------|---------|
| kind10 | esyrk | 64 | 1.7491 | 1.7449 | 1.00× |
| kind10 | esyrk | 128 | 1.0262 | 0.7302 | 1.41× |
| kind10 | esyrk | 256 | 1.3942 | 0.7434 | 1.88× |
| kind10 | esyrk | 512 | 1.7183 | 0.7519 | 2.29× |
| kind10 | ysyrk | 64 | 1.7495 | 1.7562 | 1.00× |
| kind10 | ysyrk | 128 | 1.5569 | 1.3422 | 1.16× |
| kind10 | ysyrk | 256 | 1.6179 | 1.3487 | 1.20× |
| kind10 | ysyrk | 512 | 1.6654 | 1.3491 | 1.23× |
| kind16 | qsyrk | 128 | 0.0527 | 0.0475 | 1.11× |
| kind16 | qsyrk | 256 | 0.0553 | 0.0544 | 1.02× |
| kind16 | xsyrk | 128 | 0.0487 | 0.0498 | 0.98× |
| kind16 | xsyrk | 256 | 0.0516 | 0.0521 | 0.99× |
| multifloats | msyrk | 64 | 0.4533 | 0.1126 | 4.03× |
| multifloats | msyrk | 128 | 0.8045 | 0.1209 | 6.66× |
| multifloats | msyrk | 256 | 1.3395 | 0.1226 | 10.92× |
| multifloats | msyrk | 512 | 1.9947 | 0.1235 | 16.15× |
| multifloats | wsyrk | 64 | 0.4563 | 0.2049 | 2.23× |
| multifloats | wsyrk | 128 | 0.7540 | 0.2101 | 3.59× |
| multifloats | wsyrk | 256 | 1.3524 | 0.2113 | 6.40× |
| multifloats | wsyrk | 512 | 2.0193 | 0.2142 | 9.43× |

## herk

| prec | routine | size | parallel-blas GF | migrated GF | speedup |
|------|---------|------|-----------------|-------------|---------|
| kind10 | yherk | 64 | 1.7947 | 1.7614 | 1.02× |
| kind10 | yherk | 128 | 1.4652 | 1.4034 | 1.04× |
| kind10 | yherk | 256 | 1.5985 | 1.4059 | 1.14× |
| kind10 | yherk | 512 | 1.6825 | 1.4080 | 1.19× |
| kind16 | xherk | 128 | 0.0501 | 0.0501 | 1.00× |
| kind16 | xherk | 256 | 0.0504 | 0.0518 | 0.97× |
| multifloats | wherk | 64 | 0.4802 | 0.1973 | 2.43× |
| multifloats | wherk | 128 | 0.8339 | 0.2053 | 4.06× |
| multifloats | wherk | 256 | 1.3901 | 0.1929 | 7.21× |
| multifloats | wherk | 512 | 2.0720 | 0.2115 | 9.80× |

## symm

| prec | routine | size | parallel-blas GF | migrated GF | speedup |
|------|---------|------|-----------------|-------------|---------|
| kind10 | esymm | 64 | 0.8154 | 0.7204 | 1.13× |
| kind10 | esymm | 128 | 1.3093 | 1.4556 | 0.90× |
| kind10 | esymm | 256 | 1.6416 | 1.4647 | 1.12× |
| kind10 | esymm | 512 | 1.8650 | 1.4689 | 1.27× |
| kind10 | ysymm | 64 | 2.3335 | 1.9632 | 1.19× |
| kind10 | ysymm | 128 | 2.0682 | 1.7625 | 1.17× |
| kind10 | ysymm | 256 | 2.0757 | 1.7928 | 1.16× |
| kind10 | ysymm | 512 | 1.8250 | 1.7949 | 1.02× |
| kind16 | qsymm | 128 | 0.0556 | 0.0555 | 1.00× |
| kind16 | qsymm | 256 | 0.0573 | 0.0576 | 1.00× |
| kind16 | xsymm | 128 | 0.0520 | 0.0520 | 1.00× |
| kind16 | xsymm | 256 | 0.0551 | 0.0552 | 1.00× |
| multifloats | msymm | 64 | 0.5044 | 0.1200 | 4.21× |
| multifloats | msymm | 128 | 0.8406 | 0.1219 | 6.90× |
| multifloats | msymm | 256 | 1.3374 | 0.1201 | 11.13× |
| multifloats | msymm | 512 | 2.0143 | 0.1231 | 16.37× |
| multifloats | wsymm | 64 | 0.4344 | 0.2032 | 2.14× |
| multifloats | wsymm | 128 | 0.7711 | 0.2016 | 3.82× |
| multifloats | wsymm | 256 | 1.2793 | 0.2020 | 6.33× |
| multifloats | wsymm | 512 | 1.9719 | 0.2100 | 9.39× |

## hemm

| prec | routine | size | parallel-blas GF | migrated GF | speedup |
|------|---------|------|-----------------|-------------|---------|
| kind10 | yhemm | 64 | 2.2442 | 1.9902 | 1.13× |
| kind10 | yhemm | 128 | 1.9835 | 2.1866 | 0.91× |
| kind10 | yhemm | 256 | 2.0485 | 2.2504 | 0.91× |
| kind10 | yhemm | 512 | 1.8725 | 2.2578 | 0.83× |
| kind16 | xhemm | 128 | 0.0535 | 0.0525 | 1.02× |
| kind16 | xhemm | 256 | 0.0552 | 0.0551 | 1.00× |
| multifloats | whemm | 64 | 0.4212 | 0.1990 | 2.12× |
| multifloats | whemm | 128 | 0.7342 | 0.1970 | 3.73× |
| multifloats | whemm | 256 | 1.2353 | 0.2018 | 6.12× |
| multifloats | whemm | 512 | 1.8968 | 0.2032 | 9.33× |
