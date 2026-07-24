# eplinalg

Extended-precision numerical linear algebra: a Python migrator that
retargets the Netlib BLAS/LAPACK/ScaLAPACK/MUMPS stack to quad
(`REAL(16)`), 80-bit extended (`REAL(10)`), and double-double
(multifloats) arithmetic, plus the build system and runtime support
libraries for the produced archives.

```{toctree}
:maxdepth: 1
:caption: Using the libraries

user/index
```

```{toctree}
:maxdepth: 1
:caption: Developing eplinalg

dev/index
```

```{toctree}
:maxdepth: 1
:caption: Reference

changelog
```

Also in the repository: [upstream-bugs/](upstream-bugs/README.md)
(vendored-source bug catalogue, linked from the developer guide) and
`archive/` (historical reports, not rendered).
