"""GFortran parse tree interface via subprocess.

Invokes ``gfortran -fdump-fortran-original -fsyntax-only`` to obtain a
symbol-table + code dump, then extracts the same transformation-relevant
facts as :mod:`flang_parser`:

  - Type declarations (which type, which variables)
  - Subroutine/function definitions (name, return type)
  - Call sites and function references
  - External and intrinsic declarations
  - Literal constants (real, complex)
  - Character literal constants (for XERBLA strings)

The output reuses :class:`flang_parser.ParseTreeFacts` and its component
dataclasses so that the downstream migrator is parser-agnostic.
"""

import functools
import re
import shutil
import subprocess
from pathlib import Path

from .flang_parser import (
    CallSite, ParseTreeFacts, RoutineDef, TypeDecl, UseStmtInfo
)


# ---------------------------------------------------------------------------
# Locate gfortran
# ---------------------------------------------------------------------------

@functools.cache
def find_gfortran() -> str | None:
    """Find the gfortran executable on PATH. Cached — the migrator
    calls this once per source file when --parser gfortran is set."""
    for name in ('gfortran', 'gfortran-14', 'gfortran-13', 'gfortran-12'):
        path = shutil.which(name)
        if path:
            return path
    return None


# ---------------------------------------------------------------------------
# Run gfortran
# ---------------------------------------------------------------------------

def run_gfortran_parse_tree(source_path: Path,
                            gfortran_cmd: str | None = None) -> str | None:
    """Run gfortran to get a parse tree dump.

    Returns the dump text, or None if gfortran is unavailable or fails.
    """
    if gfortran_cmd is None:
        gfortran_cmd = find_gfortran()
    if gfortran_cmd is None:
        return None

    # gfortran doesn't have a direct equivalent to -fno-sema for dumping.
    # -fdump-fortran-original always runs after resolution.
    cmd = [gfortran_cmd, '-fdump-fortran-original', '-fsyntax-only',
           str(source_path)]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=30,
        )
        # gfortran prints the dump to stdout and may still return 0 even
        # with warnings; accept any output.
        if not result.stdout:
            return None
        return result.stdout
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None


# ---------------------------------------------------------------------------
# Type spec helpers
# ---------------------------------------------------------------------------

# Maps gfortran ``(TYPE kind)`` to our canonical type-spec names.
# Canonical mapping: (gfortran_type, byte_kind) → flang-style type_spec
_TYPE_MAP: dict[tuple[str, int], str] = {
    ('REAL', 4):    'Real',
    ('REAL', 8):    'DoublePrecision',
    ('REAL', 10):   'Real',
    ('REAL', 16):   'Real',
    ('COMPLEX', 4): 'Complex',
    ('COMPLEX', 8): 'Complex',      # COMPLEX*16 (8 bytes per component)
    ('COMPLEX', 16): 'Complex',
    ('INTEGER', 4): 'Integer',
    ('INTEGER', 8): 'Integer',
    ('LOGICAL', 4): 'Logical',
    ('LOGICAL', 8): 'Logical',
}


def _parse_type_spec(line: str) -> tuple[str | None, str | None]:
    """Parse ``type spec : (TYPE kind)`` → (canonical_type, kind_value|None).

    *kind_value* is set only when a ``*N`` or ``KIND=N`` qualifier is needed
    to disambiguate (e.g., ``COMPLEX*16`` → kind_value='16').
    """
    m = re.search(r'type spec\s*:\s*\((\w+)\s+(\d+)\)', line)
    if not m:
        return None, None
    gf_type = m.group(1)
    byte_kind = int(m.group(2))
    canonical = _TYPE_MAP.get((gf_type, byte_kind))
    if canonical is None:
        return None, None

    # Determine kind_value for COMPLEX*16 and explicit-kind declarations.
    kind_value: str | None = None
    if gf_type == 'COMPLEX' and byte_kind == 8:
        kind_value = '16'          # COMPLEX*16 (double complex)
    elif gf_type == 'REAL' and byte_kind not in (4, 8):
        kind_value = str(byte_kind)
    elif gf_type == 'COMPLEX' and byte_kind not in (4, 8):
        kind_value = str(byte_kind)

    return canonical, kind_value


# ---------------------------------------------------------------------------
# Main fact extractor
# ---------------------------------------------------------------------------

# Regex for symtree line: ``  symtree: 'name'  || symbol: 'name'``
_SYMTREE_RE = re.compile(
    r"^\s+symtree:\s+'(\w+)'\s*\|\|\s*symbol:\s+'(\w+)'")

# Regex for procedure name header: ``procedure name = name``.
# CONTAINS sub-programs appear indented (`  procedure name = sub_proc`)
# while top-level procedures are at column 0; allow optional leading
# whitespace so both reset in_code and update proc_name.
_PROC_NAME_RE = re.compile(r'^\s*procedure name\s*=\s*(\w+)', re.IGNORECASE)

# Regex for CALL in code section: ``CALL name (...)``
_CALL_RE = re.compile(r'\bCALL\s+(\w+)\s*\(')
_USE_STMT_RE = re.compile(
    r'^\s*USE(?:\s*,[^:]*)?\s*(?:::)?\s+(\w+)', re.IGNORECASE)
_USE_ASSOC_RE = re.compile(r'USE-ASSOC\(([\w]+)\)')
_ATTRS_RE = re.compile(r'attributes:\s*\(([^)]*)\)')

# Regex for function references: ``name[[ ... ]]``
_FUNCREF_RE = re.compile(r'\b(\w+)\[\[')

# Regex for real literals in code: ``digits.digitsEdigits_kind``
_REAL_LIT_RE = re.compile(
    r'(?<!\w)(\d+\.\d+(?:e[+-]?\d+)?_\d+)(?!\w)', re.IGNORECASE)


def parse_tree_facts(tree_text: str) -> ParseTreeFacts:
    """Extract transformation-relevant facts from gfortran dump text.

    The gfortran ``-fdump-fortran-original`` output has two sections per
    program unit:

    1. **Symbol table** – ``symtree`` entries with ``type spec`` and
       ``attributes`` for every symbol (variables, procedures, intrinsics).
    2. **Code section** – starting after the ``code:`` line, containing
       a linearised representation of all executable statements.

    We walk the symbol table to populate type declarations, routine
    definitions, external names, and intrinsic names.  We scan the code
    section for call sites, function references, character literals
    (XERBLA strings), and real literal constants.
    """
    facts = ParseTreeFacts()
    lines = tree_text.splitlines()

    proc_name: str | None = None  # current procedure
    in_code = False

    # --- Collected symbols (before we turn them into facts) ---
    # sym_name → {type_spec, kind_value, attrs_set, is_procedure_sym}
    symbols: dict[str, dict] = {}
    call_names: set[str] = set()
    func_ref_names: set[str] = set()
    # Module-name set mirroring facts.use_stmt_ranges for O(1) duplicate
    # detection. The previous ``any(u.module_name == mod_name ...)``
    # walk made USE-ASSOC harvesting O(M²) in #modules.
    seen_use_modules: set[str] = set()

    i = 0
    while i < len(lines):
        line = lines[i]

        # --- Procedure header ---
        pm = _PROC_NAME_RE.match(line)
        if pm:
            proc_name = pm.group(1).lower()
            in_code = False
            i += 1
            continue

        # --- Start of code section ---
        if line.strip() == 'code:':
            in_code = True
            i += 1
            continue

        if not in_code:
            # --- Use statement ---
            um = _USE_STMT_RE.match(line)
            if um:
                mod_name = um.group(1).upper()
                facts.use_stmt_ranges.append(
                    UseStmtInfo(mod_name))
                seen_use_modules.add(mod_name)

            # --- USE-ASSOC tag ---
            # Symbols imported from another module appear with
            # ``USE-ASSOC(modname)`` in their attributes line. The
            # `USE`-statement regex above only fires for explicit
            # ``USE`` lines, which gfortran's dump doesn't always
            # emit; harvesting USE-ASSOC tags from the symbol table
            # ensures the rename_map covers module names that the
            # caller imports a TYPE/symbol from.
            for am in _USE_ASSOC_RE.finditer(line):
                mod_name = am.group(1).upper()
                if mod_name not in seen_use_modules:
                    facts.use_stmt_ranges.append(
                        UseStmtInfo(mod_name))
                    seen_use_modules.add(mod_name)

            # --- Symbol table entry ---
            sm = _SYMTREE_RE.match(line)
            if sm:
                sym_name = sm.group(2).lower()
                # Read subsequent indented lines for type spec & attributes
                type_spec = None
                kind_value = None
                attrs_raw = ''
                j = i + 1
                while j < len(lines) and lines[j].startswith('    '):
                    sline = lines[j]
                    if 'type spec' in sline and type_spec is None:
                        type_spec, kind_value = _parse_type_spec(sline)
                    if 'attributes:' in sline:
                        attrs_raw = sline
                    # Harvest USE-ASSOC(modname) tags from any
                    # attribute line — these mark symbols imported
                    # from another module, which we need in
                    # use_stmt_ranges so callers' rename_map keeps
                    # the imported module name.
                    for am in _USE_ASSOC_RE.finditer(sline):
                        mod_name = am.group(1).upper()
                        if mod_name not in seen_use_modules:
                            facts.use_stmt_ranges.append(
                                UseStmtInfo(mod_name))
                            seen_use_modules.add(mod_name)
                    j += 1

                attrs = _parse_attributes(attrs_raw)
                # User-defined derived types appear twice in gfortran's
                # dump: once as ``symtree: 'Dmumps_intr_struc'`` (case-
                # preserved, attributes (DERIVED ...)``) and once as
                # ``symtree: 'dmumps_intr_struc'`` (autogenerated generic
                # interface, attributes (PROCEDURE FUNCTION ...)).
                # Both share the same lowercase ``symbol:`` key, so a
                # plain dict overwrite drops the type definition. Emit
                # the type entry immediately into routine_defs and skip
                # the dict mutation for it.
                if 'DERIVED' in attrs and not any(
                        a.startswith('USE-ASSOC') for a in attrs):
                    facts.routine_defs.append(
                        RoutineDef('type', sym_name.upper()))
                    i = j
                    continue
                symbols[sym_name] = {
                    'type_spec': type_spec,
                    'kind_value': kind_value,
                    'attrs': attrs,
                    # Record the procedure block this symbol was defined
                    # in. Without this, is_proc_sym below compares each
                    # symbol against the *last* proc_name seen in the
                    # file (the innermost CONTAINS procedure), so a
                    # top-level subroutine/function whose body has a
                    # CONTAINS subprogram never matches its own
                    # procedure block and gets dropped from routine_defs.
                    'proc_name': proc_name,
                }
                i = j
                continue
        else:
            # --- Code section: scan for calls, references, literals ---
            for cm in _CALL_RE.finditer(line):
                name = cm.group(1).lower()
                # Skip scope-qualified prefix (e.g. ``dgemm:``)
                if ':' in name:
                    continue
                call_names.add(name)

            for fm in _FUNCREF_RE.finditer(line):
                name = fm.group(1).lower()
                if ':' in name:
                    # e.g. ``dgemm:temp`` is a variable, not a function
                    continue
                func_ref_names.add(name)

            for rm in _REAL_LIT_RE.finditer(line):
                facts.real_literals.append(rm.group(1))

        i += 1

    # --- Convert symbol table into facts ---
    for sym_name, info in symbols.items():
        attrs = info['attrs']
        ts = info['type_spec']
        kv = info['kind_value']

        # Skip the procedure's own entry (e.g. ``symtree: 'dgemm'`` for
        # the dgemm subroutine) — we handle that via the ``procedure
        # name`` header plus attributes. Compare against the procedure
        # block this symbol was defined in (saved at parse time), not
        # the last proc_name seen in the file.
        sym_proc = info.get('proc_name')
        is_proc_sym = (sym_name == sym_proc)
        # Also include MODULE-PROC contained procedures (subroutines /
        # functions defined inside a MODULE block). Without this, a file
        # like dstatic_ptr_m.F whose proc_name is the MODULE name leaves
        # routine_defs empty, file_rename_map gets filtered to nothing,
        # and the contained subroutines (DMUMPS_GET_TMP_PTR etc.) and
        # the MODULE name itself never get renamed.
        is_module_proc = ('MODULE-PROC' in attrs)
        is_module_name = ('MODULE' in attrs and 'PROCEDURE' not in attrs
                          and 'SUBROUTINE' not in attrs
                          and 'FUNCTION' not in attrs)

        if 'SUBROUTINE' in attrs and (is_proc_sym or is_module_proc):
            facts.routine_defs.append(
                RoutineDef('subroutine', sym_name.upper()))
        elif 'FUNCTION' in attrs and (is_proc_sym or is_module_proc):
            facts.routine_defs.append(
                RoutineDef('function', sym_name.upper()))
        elif is_module_name:
            # Treat the module name itself like a routine_def so the
            # rename_map filter in _migrate_with_flang keeps the entry.
            facts.routine_defs.append(
                RoutineDef('module', sym_name.upper()))

        if 'EXTERNAL' in attrs and 'INTRINSIC' not in attrs:
            facts.external_names.append(sym_name.upper())

        # Variable declarations (non-procedure, non-intrinsic symbols)
        if ts is not None:
            if ts in ('Real', 'DoublePrecision', 'Complex', 'Integer', 'Logical'):
                facts.variable_types[sym_name.upper()] = ts

            if 'VARIABLE' in attrs:
                if ts in ('Real', 'DoublePrecision', 'Complex'):
                    facts.type_decls.append(
                        TypeDecl(ts, kv, [sym_name.upper()]))

    # --- Call sites from code section ---
    for name in sorted(call_names):
        facts.call_sites.append(CallSite(name.upper()))

    seen_call_names = {cs.name for cs in facts.call_sites}
    for name in sorted(func_ref_names):
        uname = name.upper()
        # Skip gfortran internal helper symbols (e.g. __max_i4, __convert_*)
        if uname.startswith('__'):
            continue
        if uname not in seen_call_names:
            facts.call_sites.append(CallSite(uname))
            seen_call_names.add(uname)

    return facts


def _parse_attributes(attrs_line: str) -> set[str]:
    """Extract attribute keywords from an ``attributes:`` line.

    Example input::

        attributes: (PROCEDURE EXTERNAL-PROC  EXTERNAL SUBROUTINE ARRAY-OUTER-DEPENDENCY)

    Returns ``{'PROCEDURE', 'EXTERNAL-PROC', 'EXTERNAL', 'SUBROUTINE', ...}``
    """
    m = _ATTRS_RE.search(attrs_line)
    if not m:
        return set()
    return set(m.group(1).split())


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def scan_file(source_path: Path,
              gfortran_cmd: str | None = None) -> ParseTreeFacts | None:
    """Scan a Fortran file using gfortran and return extracted facts.

    Returns None if gfortran is not available or fails.
    """
    tree_text = run_gfortran_parse_tree(source_path, gfortran_cmd)
    if tree_text is None:
        return None
    return parse_tree_facts(tree_text)
