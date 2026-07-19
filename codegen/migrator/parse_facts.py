"""Parser-agnostic fact model shared by the parse-tree scanners.

:mod:`flang_parser` and :mod:`gfortran_parser` both reduce a compiler
parse-tree dump to the same :class:`ParseTreeFacts` structure, so the
downstream migrator (:func:`fortran_migrator._migrate_with_facts`) is
parser-agnostic.
"""

from dataclasses import dataclass, field


@dataclass
class TypeDecl:
    """A type declaration found in the parse tree."""
    type_spec: str          # e.g., 'DoublePrecision', 'Real', 'Complex'
    kind_value: str | None  # e.g., '16' for COMPLEX*16, None for plain REAL
    names: list[str]        # declared entity names


@dataclass
class RoutineDef:
    """A subroutine, function, module, or derived-type definition."""
    kind: str               # 'subroutine', 'function', 'module' or 'type'
    name: str


@dataclass
class ParseTreeFacts:
    """All facts extracted from a compiler parse-tree dump.

    ``call_sites`` holds upper-cased called-routine / referenced-function
    names in discovery order; ``use_modules`` holds the upper-cased names
    of every module the file references via ``USE`` (or, for gfortran,
    via a ``USE-ASSOC(modname)`` attribute tag).
    """
    type_decls: list[TypeDecl] = field(default_factory=list)
    routine_defs: list[RoutineDef] = field(default_factory=list)
    call_sites: list[str] = field(default_factory=list)
    external_names: list[str] = field(default_factory=list)
    real_literals: list[str] = field(default_factory=list)
    variable_types: dict[str, str] = field(default_factory=dict)
    use_modules: set[str] = field(default_factory=set)
