"""Sanity checks on the intrinsic-rename tables.

These tables are large hand-curated lookups; the value of testing them
is structural — guard against accidental edits that break the contract
the rest of the migrator relies on (every value is a 2-tuple of
(generic_name, needs_kind_arg)).
"""

import re

from migrator.intrinsics import INTRINSIC_MAP, INTRINSIC_DECL_MAP


def test_intrinsic_map_value_shape():
    """Every entry must be (str, bool)."""
    for old, val in INTRINSIC_MAP.items():
        assert isinstance(val, tuple) and len(val) == 2, (
            f'{old!r} → {val!r}: must be (new_name, needs_kind_arg) tuple'
        )
        new_name, needs_kind = val
        assert isinstance(new_name, str) and new_name, (
            f'{old!r}: empty/non-string new name'
        )
        assert isinstance(needs_kind, bool), (
            f'{old!r}: needs_kind must be bool, got {type(needs_kind)}'
        )


def test_intrinsic_map_keys_uppercase_identifiers():
    """Source intrinsic names are matched case-insensitively in the
    migrator, but the table is canonicalised in upper-case so callers
    can do a direct dict lookup on ``name.upper()``."""
    ident = re.compile(r'^[A-Z][A-Z0-9_]*$')
    for old in INTRINSIC_MAP:
        assert ident.match(old), f'{old!r}: not an upper-case identifier'


def test_well_known_real_intrinsic_renames():
    """Spot-check a handful of the most-used renames so a
    refactor that drops an entry shows up immediately."""
    assert INTRINSIC_MAP['DABS'] == ('ABS', False)
    assert INTRINSIC_MAP['DSQRT'] == ('SQRT', False)
    assert INTRINSIC_MAP['DCONJG'] == ('CONJG', False)
    assert INTRINSIC_MAP['DIMAG'] == ('AIMAG', False)


def test_decl_map_values_uppercase():
    """The declaration-rename table maps old type-specific declared
    names (``DOUBLE PRECISION FUNCTION DABS``) to their generic
    counterpart; values must be upper-case identifiers."""
    ident = re.compile(r'^[A-Z][A-Z0-9_]*$')
    for old, new in INTRINSIC_DECL_MAP.items():
        assert ident.match(new), f'{old}→{new}: bad declaration target'
