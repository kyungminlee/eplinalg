"""Fortran migration engine, split into single-concern submodules.

Extracted verbatim from the former ``fortran_migrator.py`` monolith as part of
the migrator file-restructuring refactor. Layering is acyclic and bottom-up::

    lex  <-  per-line rewriter clusters  <-  form orchestrators

``fortran_migrator.py`` remains the public facade and re-imports the moved
names, so every existing caller keeps working unchanged.
"""
