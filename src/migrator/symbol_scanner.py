"""Symbol scanning from source files or compiled libraries.

Extracts SUBROUTINE/FUNCTION names from Fortran source, or function
names from C source, or defined symbols from compiled archives via nm.
"""

import re
import subprocess
from pathlib import Path

# Fortran patterns. ``MODULE`` is included so per-arithmetic module
# names (e.g. ``DMUMPS_STATIC_PTR_M`` / ``SMUMPS_STATIC_PTR_M``)
# collapse to their target-prefixed form (``QMUMPS_STATIC_PTR_M``)
# like subroutines. ``MODULE PROCEDURE`` and ``END MODULE`` are
# excluded — they reference existing names, not new definitions.
_FORTRAN_DEF_RE = re.compile(
    r'(?<!\w)(?<!END\s)(?<!end\s)'
    r'(?:SUBROUTINE|FUNCTION|MODULE(?!\s+(?:PROCEDURE|SUBROUTINE|FUNCTION)\b))'
    r'\s+([A-Za-z]\w*)',
    re.IGNORECASE,
)

# Derived-type definitions: ``TYPE FOO``, ``TYPE :: FOO``, ``TYPE,
# attr :: FOO``. Anchored so ``TYPE(FOO)`` (a type reference, not a
# definition) is rejected. Needed so precision-prefixed types like
# MUMPS's ``DMUMPS_STRUC`` — defined in dmumps_struc.h — are picked up
# and renamed alongside routines in the same family.
_FORTRAN_TYPE_DEF_RE = re.compile(
    r'^\s*TYPE\b(?!\s*\()(?:\s*,[^:]*)?(?:\s*::)?\s*'
    r'(?!(?:IS|DEFAULT)\b)([A-Za-z]\w*)\s*(?:!.*)?$',
    re.IGNORECASE,
)

# Built-in C return types recognized by the function-definition scanner.
# Recipes may add more via the ``c_return_types`` YAML key; the final
# regex is constructed by :func:`_build_c_func_re`.
_C_DEFAULT_RETURN_TYPES: tuple[str, ...] = (
    'void', 'int', 'float', 'double', 'Int', 'BVOID',
    'F_VOID_FUNC', 'F_INT_FUNC', 'F_DOUBLE_FUNC',
    'SCOMPLEX', 'DCOMPLEX', 'QCOMPLEX', 'MPI_Datatype',
    r'BLACBUFF\s*\*',
)


def _build_c_func_re(extra_return_types: tuple[str, ...] = ()) -> re.Pattern:
    """Compile the C function-definition regex for a given type set."""
    types = list(_C_DEFAULT_RETURN_TYPES) + list(extra_return_types)
    alt = '|'.join(types)
    return re.compile(
        r'^\s*(?:' + alt + r')\s+(\w+)\s*\(',
        re.MULTILINE,
    )


# #define alias name_ — exposes ``name_`` as the Fortran-callable entry
# (the preprocessor substitutes the alias used in the function
# definition to ``name_``). Used by ScaLAPACK's REDIST/SRC to export
# ``psgemr2d_`` et al. Also catches the reverse form ``#define name_
# other_alias`` in case a recipe uses it.
_C_DEFINE_FORTRAN_RE = re.compile(
    r'^\s*#\s*define\s+(\w+)\s+(\w+)\s*$', re.MULTILINE,
)


# Default compiled pattern (no recipe-specific extensions).
_C_FUNC_RE = _build_c_func_re()


def scan_fortran_source(src_dir: Path,
                        extensions: list[str] | None = None) -> set[str]:
    """Scan Fortran source files for SUBROUTINE/FUNCTION definitions."""
    if extensions is None:
        extensions = ['.f', '.f90', '.for', '.f95']

    names: set[str] = set()
    for f in sorted(src_dir.iterdir()):
        if f.suffix.lower() not in extensions:
            continue
        for line in f.read_text(errors='replace').splitlines():
            # Skip fixed-form comments
            if line and line[0] in ('C', 'c', '*', '!'):
                continue
            # Skip free-form comments
            stripped = line.lstrip()
            if stripped.startswith('!'):
                continue
            m = _FORTRAN_DEF_RE.search(line)
            if m:
                names.add(m.group(1).upper())
                continue
            m = _FORTRAN_TYPE_DEF_RE.match(line)
            if m:
                names.add(m.group(1).upper())
    return names


def scan_c_source(src_dir: Path,
                  extensions: list[str] | None = None,
                  extra_return_types: tuple[str, ...] = ()) -> set[str]:
    """Scan C source files for function definitions.

    ``extra_return_types`` extends the built-in set with recipe-specific
    return types (as regex fragments, e.g. ``r'PBTYP_T\\s*\\*'``).
    """
    if extensions is None:
        extensions = ['.c']

    pattern = (_build_c_func_re(tuple(extra_return_types))
               if extra_return_types else _C_FUNC_RE)

    names: set[str] = set()
    for f in sorted(src_dir.iterdir()):
        if f.suffix.lower() not in extensions:
            continue
        text = f.read_text(errors='replace')
        for m in pattern.finditer(text):
            sym = m.group(1)
            # Strip trailing underscore (Fortran name mangling in C wrappers)
            if sym.endswith('_'):
                sym = sym[:-1]
            if sym:
                names.add(sym.upper())
        # #define alias name_  — ``name_`` is the resulting Fortran
        # entry-point symbol. Either side may carry the trailing
        # underscore depending on how the file is written.
        for m in _C_DEFINE_FORTRAN_RE.finditer(text):
            for tok in (m.group(1), m.group(2)):
                if tok.endswith('_') and len(tok) > 1:
                    names.add(tok[:-1].upper())
    return names


def scan_library_archive(lib_path: Path) -> set[str]:
    """Extract defined symbols from a static or shared library using nm."""
    cmd = ['nm', '--defined-only', '-g', str(lib_path)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True,
                                timeout=30)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return set()

    nm_re = re.compile(r'\s*[0-9a-fA-F]*\s+[TtWw]\s+(\S+)')
    names: set[str] = set()
    for line in result.stdout.splitlines():
        m = nm_re.match(line)
        if m:
            sym = m.group(1)
            # Strip leading underscore (macOS) and trailing underscore (Fortran)
            if sym.startswith('_'):
                sym = sym[1:]
            if sym.endswith('_'):
                sym = sym[:-1]
            if sym:
                names.add(sym.upper())
    return names


def scan_symbols(source_dir: Path, language: str,
                 extensions: list[str] | None = None,
                 library_path: Path | None = None,
                 extra_c_return_types: tuple[str, ...] = ()) -> set[str]:
    """Scan for symbols using the appropriate method.

    If library_path is provided, uses nm. Otherwise scans source.
    ``extra_c_return_types`` is forwarded to :func:`scan_c_source` when
    ``language == 'c'``.
    """
    if library_path and library_path.exists():
        return scan_library_archive(library_path)

    if language == 'fortran':
        return scan_fortran_source(source_dir, extensions)
    elif language == 'c':
        return scan_c_source(source_dir, extensions, extra_c_return_types)
    else:
        raise ValueError(f'Unsupported language: {language}')
