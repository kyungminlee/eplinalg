"""Fortran migration engine, split into single-concern submodules.

Extracted verbatim from the former ``fortran_migrator.py`` monolith as part of
the migrator file-restructuring refactor. Layering is acyclic and bottom-up::

    lex  <-  per-line rewriter clusters  <-  form orchestrators

``fortran_migrator.py`` hosts only the top-level drivers that sequence
these passes and imports just the names those drivers use; callers reach
per-pass helpers through the owning submodule directly.
"""
