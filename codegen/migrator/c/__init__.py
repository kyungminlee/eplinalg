"""Library-specific C migration passes.

Mirrors the ``fortran/`` subpackage layout: the generic clone/rename
engine stays in ``migrator.c_migrator``; per-library knowledge (BLACS
directory handling and Bdef.h ABI patching, PBLAS header/idiom rules)
lives here, one module per library, with the shared clone primitives
in ``clone``.
"""
