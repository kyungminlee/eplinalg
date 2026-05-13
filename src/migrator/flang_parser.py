"""Flang parse tree interface via subprocess.

Invokes `flang-new -fc1 -fdebug-dump-parse-tree-no-sema` to get a
structured parse tree, then extracts transformation-relevant facts:
  - Type declarations (which type, which variables)
  - Subroutine/function definitions (name, return type)
  - Call sites and function references
  - External and intrinsic declarations
  - Literal constants (real, complex)
  - Character literal constants (for XERBLA strings)

This is used by the Fortran migrator as a "scanner oracle" — the parse
tree tells us WHAT needs transformation, then source-level regex applies
the actual byte-range replacements (preserving formatting).
"""

import functools
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class TypeDecl:
    """A type declaration found by Flang."""
    type_spec: str          # e.g., 'DoublePrecision', 'Real', 'Complex'
    kind_value: str | None  # e.g., '16' for COMPLEX*16, None for plain REAL
    names: list[str]        # declared entity names


@dataclass
class RoutineDef:
    """A subroutine or function definition."""
    kind: str               # 'subroutine' or 'function'
    name: str
    return_type: str | None # For functions: the type spec


@dataclass
class CallSite:
    """A CALL statement or function reference."""
    name: str               # called routine name
    is_call_stmt: bool      # True for CALL, False for function reference


@dataclass
class ParameterInfo:
    """Information about a PARAMETER statement."""
    name: str
    value: str
    type_spec: str          # 'FP' or 'integer' (simplified)
    line_number: int | None = None


@dataclass
class DataInfo:
    """Information about a DATA statement."""
    names: list[str]
    values: list[str]
    types: list[str] | None = None
    line_number: int | None = None


@dataclass
class ProcInfo:
    """Information about procedure boundaries."""
    name: str
    header_line: int
    first_exec_line: int | None
    end_line: int


@dataclass
class UseStmtInfo:
    """Information about a USE statement."""
    module_name: str
    line_range: tuple[int, int]
    only_list: list[str]


# The Flang parse-tree dump emits ``Name = 'IDENT'`` markers on
# almost every interesting line; the same regex is searched at nine
# call sites in this module. Hoist once.
_NAME_RE = re.compile(r"Name\s*=\s*'(\w+)'")


@dataclass
class ParseTreeFacts:
    """All facts extracted from a Flang parse tree."""
    type_decls: list[TypeDecl] = field(default_factory=list)
    routine_defs: list[RoutineDef] = field(default_factory=list)
    call_sites: list[CallSite] = field(default_factory=list)
    external_names: list[str] = field(default_factory=list)
    intrinsic_names: list[str] = field(default_factory=list)
    real_literals: list[str] = field(default_factory=list)
    char_literals: list[str] = field(default_factory=list)

    # New fields for Phase 1.5
    parameter_stmts: list[ParameterInfo] = field(default_factory=list)
    data_stmts: list[DataInfo] = field(default_factory=list)
    procedure_boundaries: list[ProcInfo] = field(default_factory=list)
    variable_types: dict[str, str] = field(default_factory=dict)
    use_stmt_ranges: list[UseStmtInfo] = field(default_factory=list)


@functools.cache
def find_flang() -> str | None:
    """Find the flang-new executable. Cached — shutil.which stats $PATH
    on every invocation, and the migrator calls this once per source
    file when --parser flang is set."""
    for name in ('flang-new', 'flang'):
        path = shutil.which(name)
        if path:
            return path
    return None


def run_flang_parse_tree(source_path: Path,
                         flang_cmd: str | None = None,
                         use_sema: bool = False) -> str | None:
    """Run Flang to get a parse tree dump.

    Returns the parse tree text, or None if Flang fails.
    """
    if flang_cmd is None:
        flang_cmd = find_flang()
    if flang_cmd is None:
        return None

    dump_flag = '-fdebug-dump-parse-tree' if use_sema else '-fdebug-dump-parse-tree-no-sema'
    cmd = [flang_cmd, '-fc1', dump_flag, str(source_path)]
    if use_sema:
        cmd.append('-fdebug-resolve-names')
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0:
            return None
        return result.stdout
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None


def parse_tree_facts(tree_text: str) -> ParseTreeFacts:
    """Extract transformation-relevant facts from Flang parse tree text.

    This is a line-by-line pattern matcher, not a full tree parser.
    It extracts the specific facts needed for type migration.
    """
    facts = ParseTreeFacts()

    lines = tree_text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]

        # --- Subroutine definition ---
        if 'SubroutineStmt' in line:
            name = _extract_first_name(lines, i + 1)
            if name:
                facts.routine_defs.append(
                    RoutineDef('subroutine', name, None))
                facts.procedure_boundaries.append(
                    ProcInfo(name, 0, None, 0))

        # --- Function definition ---
        elif 'FunctionStmt' in line:
            ret_type = None
            # Check for return type prefix on the same line
            if 'DoublePrecision' in line:
                ret_type = 'DoublePrecision'
            elif 'Complex' in line:
                ret_type = 'Complex'
            elif 'Real' in line and 'IntrinsicTypeSpec' in line:
                ret_type = 'Real'
            name = _extract_first_name(lines, i + 1)
            if name:
                facts.routine_defs.append(
                    RoutineDef('function', name, ret_type))
                facts.procedure_boundaries.append(
                    ProcInfo(name, 0, None, 0))

        # --- Type declaration ---
        elif 'TypeDeclarationStmt' in line:
            type_spec, kind_val, names = _parse_type_decl(lines, i)
            if type_spec and names:
                facts.type_decls.append(
                    TypeDecl(type_spec, kind_val, names))
                for name in names:
                    facts.variable_types[name] = type_spec

        # --- Use statement ---
        elif 'UseStmt' in line:
            m = _NAME_RE.search(line)
            if not m and i + 1 < len(lines):
                m = _NAME_RE.search(lines[i + 1])
            if m:
                mod_name = m.group(1).upper()
                facts.use_stmt_ranges.append(
                    UseStmtInfo(mod_name, (0, 0), []))

        # --- Parameter statement ---
        elif 'ParameterStmt' in line:
            p_name = _extract_first_name(lines, i + 1)
            if p_name:
                facts.parameter_stmts.append(
                    ParameterInfo(p_name, "unknown", "FP"))

        # --- Data statement ---
        elif 'DataStmt' in line:
            facts.data_stmts.append(DataInfo([], [], None))

        # --- Call statement ---
        elif 'CallStmt' in line:
            name = _extract_procedure_name(lines, i)
            if name:
                facts.call_sites.append(CallSite(name, True))

        # --- Function reference (not inside CallStmt) ---
        elif 'FunctionReference -> Call' in line and 'CallStmt' not in line:
            name = _extract_procedure_name(lines, i)
            if name:
                facts.call_sites.append(CallSite(name, False))

        # --- External declaration ---
        elif 'ExternalStmt' in line:
            m = _NAME_RE.search(line)
            if m:
                facts.external_names.append(m.group(1).upper())

        # --- Intrinsic declaration ---
        elif 'IntrinsicStmt' in line:
            m = _NAME_RE.search(line)
            if m:
                facts.intrinsic_names.append(m.group(1).upper())

        # --- Real literal constant ---
        elif 'RealLiteralConstant' in line:
            m = re.search(r"Real\s*=\s*'([^']+)'", line)
            if not m:
                # Check next line
                if i + 1 < len(lines):
                    m = re.search(r"Real\s*=\s*'([^']+)'", lines[i + 1])
            if m:
                facts.real_literals.append(m.group(1))

        # --- Character literal constant ---
        elif 'CharLiteralConstant' in line:
            m = re.search(r"string\s*=\s*'([^']*)'", line)
            if not m and i + 1 < len(lines):
                m = re.search(r"string\s*=\s*'([^']*)'", lines[i + 1])
            if m:
                facts.char_literals.append(m.group(1))

        i += 1

    # Also scan for all ProcedureDesignator names (function references
    # that appear as sub-expressions, like DCONJG, LSAME, etc.).
    # Use a set for O(1) duplicate detection — the previous
    # ``any(cs.name == name ...)`` scan was O(N) per ProcedureDesignator,
    # quadratic in total over the parse tree.
    seen_call_names = {cs.name for cs in facts.call_sites}
    for line in lines:
        if 'ProcedureDesignator -> Name' in line:
            m = _NAME_RE.search(line)
            if m:
                name = m.group(1).upper()
                if name not in seen_call_names:
                    facts.call_sites.append(CallSite(name, False))
                    seen_call_names.add(name)

    return facts


def _extract_first_name(lines: list[str], start: int) -> str | None:
    """Extract the first Name = '...' value starting from a line index."""
    for j in range(start, min(start + 5, len(lines))):
        m = _NAME_RE.search(lines[j])
        if m:
            return m.group(1).upper()
    return None


def _extract_procedure_name(lines: list[str], start: int) -> str | None:
    """Extract ProcedureDesignator -> Name from around a line index."""
    for j in range(start, min(start + 10, len(lines))):
        if 'ProcedureDesignator' in lines[j]:
            m = _NAME_RE.search(lines[j])
            if m:
                return m.group(1).upper()
    return None


def _parse_type_decl(lines: list[str], start: int) -> tuple:
    """Parse a TypeDeclarationStmt block.

    Returns (type_spec, kind_value, [entity_names]).
    """
    type_spec = None
    kind_val = None
    names = []

    # Determine the indentation level of the TypeDeclarationStmt
    base_indent = _indent_level(lines[start])

    for j in range(start, min(start + 30, len(lines))):
        line = lines[j]
        if j > start and _indent_level(line) <= base_indent:
            break  # Left the TypeDeclarationStmt block

        if 'DoublePrecision' in line:
            type_spec = 'DoublePrecision'
        elif 'IntrinsicTypeSpec -> Real' in line:
            type_spec = 'Real'
        elif 'IntrinsicTypeSpec -> Complex' in line:
            type_spec = 'Complex'
        elif 'IntrinsicTypeSpec -> IntegerTypeSpec' in line:
            type_spec = 'Integer'
        elif 'IntrinsicTypeSpec -> Logical' in line:
            type_spec = 'Logical'
        elif 'IntrinsicTypeSpec -> Character' in line:
            type_spec = 'Character'

        # Kind selector (e.g., COMPLEX*16)
        if 'StarSize' in line:
            m = re.search(r"uint64_t\s*=\s*'(\d+)'", line)
            if m:
                kind_val = m.group(1)
        elif 'KindSelector' in line and 'Scalar' in line:
            m = re.search(r"uint64_t\s*=\s*'(\d+)'", line)
            if m:
                kind_val = m.group(1)

        # Entity names
        if 'EntityDecl' in line:
            # Name is usually on the next line or same line
            m = _NAME_RE.search(line)
            if m:
                names.append(m.group(1).upper())
            elif j + 1 < len(lines):
                m = _NAME_RE.search(lines[j + 1])
                if m:
                    names.append(m.group(1).upper())

    return type_spec, kind_val, names


def _indent_level(line: str) -> int:
    """Count the indentation level based on '| ' prefixes."""
    count = 0
    i = 0
    while i < len(line) - 1:
        if line[i] == '|' and line[i + 1] == ' ':
            count += 1
            i += 2
        else:
            break
    return count


def scan_file(source_path: Path,
              flang_cmd: str | None = None,
              use_sema: bool = False) -> ParseTreeFacts | None:
    """Scan a Fortran file using Flang and return extracted facts.

    Returns None if Flang is not available or fails.
    """
    tree_text = run_flang_parse_tree(source_path, flang_cmd, use_sema)
    if tree_text is None:
        return None
    return parse_tree_facts(tree_text)
