"""Precision prefix classification via empirical pattern discovery.

Discovers precision structure by analyzing the full symbol set in two passes:

  Pass 1 (single precision): Replace S→#R, C→#C at candidate positions.
  Pass 2 (double precision): Replace D→#R, Z→#C at candidate positions.

Each replaced position carries a **type tag** (R=real, C=complex) so
that multi-position patterns like DZNRM2→#C#RNRM2 know position 0 is
complex and position 1 is real — and never collide with an unrelated
pattern where the two slots carry different types.

A family is formed by overlapping pass-1 and pass-2 results that share
the same tagged pattern AND the same slot positions. Pass-1 members
(S/C) come from single-precision source; pass-2 members (D/Z) come
from double-precision source. **Only the D/Z members drive renames**
— they are the extension sources (D→Q/E, Z→X/Y). S/C symbols keep
their original names and are copied through.

Because type tags are part of the grouping key, real and complex
families are always distinct:
  - (SGEMM, DGEMM) forms a real family with pattern "#R GEMM".
  - (CGEMM, ZGEMM) forms a complex family with pattern "#C GEMM".
  - (CSROT, ZDROT) forms a cross-type family with pattern
    "#C#R ROT" (pos 0 complex, pos 1 real).

This naturally handles any prefix convention without hardcoding:
  - Direct:      SGEMM/DGEMM, CGEMM/ZGEMM
  - I-prefix:    ISAMAX/IDAMAX, ICAMAX/IZAMAX
  - ScaLAPACK:   PSGESV/PDGESV, PCGESV/PZGESV
  - Cross-type:  CSROT/ZDROT
"""

from __future__ import annotations

from dataclasses import dataclass, field

# Character → type tag. Real and complex never merge: pass-1 and pass-2
# only share a tagged pattern when their tags agree position-by-position.
CHAR_TYPE: dict[str, str] = {
    'S': 'R', 'D': 'R',    # real
    'C': 'C', 'Z': 'C',    # complex
}

# Two passes: single-precision chars and double-precision chars
SINGLE_CHARS = {'S', 'C'}
DOUBLE_CHARS = {'D', 'Z'}


@dataclass
class PrecisionSlot:
    """One position in a symbol that carries a precision character."""
    position: int
    type_tag: str    # 'R' (real) or 'C' (complex)


@dataclass
class PrecisionFamily:
    """A group of symbols sharing the same tagged pattern."""
    pattern: str                     # e.g., '#RGEMM', '#C#RROT'
    members: dict[str, str]          # precision_key → symbol name
    slots: list[PrecisionSlot]       # positions + type tags


@dataclass
class SymbolClassification:
    """Result of analyzing a set of symbols for precision patterns."""
    families: list[PrecisionFamily] = field(default_factory=list)
    independent: set[str] = field(default_factory=set)
    _symbol_to_family: dict[str, PrecisionFamily] = field(
        default_factory=dict, repr=False)

    def get_family(self, name: str) -> PrecisionFamily | None:
        return self._symbol_to_family.get(name.upper())

    def target_name(self, name: str, target_mode) -> str | None:
        """Compute the target name for a symbol at the given mode.

        Every member of a family maps to the same target name: each
        slot's char is replaced by the type-tag-driven target (R→Q/E,
        C→X/Y). So SGEMM and DGEMM both map to QGEMM; CGEMM and ZGEMM
        both map to XGEMM. Collisions between co-family members are
        expected and checked for convergence at migration time.

        Cross-family collisions (two distinct upstream families
        producing the same target name) are caught by
        ``build_rename_map`` below and surface as a RuntimeError so
        the offending target prefix can be re-picked.
        """
        upper = name.upper()
        family = self._symbol_to_family.get(upper)
        if family is None:
            return None

        pmap = target_mode.prefix_map
        result = list(upper)
        for slot in family.slots:
            result[slot.position] = pmap[slot.type_tag]
        return ''.join(result)

    def build_rename_map(self, target_mode) -> dict[str, str]:
        """Build old_name → new_name mapping for every precision-dependent symbol."""
        rename: dict[str, str] = {}
        for sym in self._symbol_to_family:
            new_name = self.target_name(sym, target_mode)
            if new_name and new_name != sym:
                rename[sym] = new_name

        # Collision diagnostic 1: an orphaned (un-renamed) symbol whose
        # name equals another symbol's rename target will silently
        # shadow that target at link time. We hit this with multifloats'
        # DD/ZZ prefixes: BLAS DDOT renamed to DDDOT collided with
        # ScaLAPACK's orphaned DDDOT wrapper, corrupting pddpotf2.
        # Surface the conflict so the developer can investigate
        # (re-pick the prefix, force-rename, etc.) rather than letting
        # the migrator emit a silent miscompile.
        #
        # Collision diagnostic 2: two distinct source symbols both
        # rename to the same target. We hit this with multifloats'
        # D→T, Z→V where DVASUM and DZASUM both rename to TVASUM.
        # That's a hard merge; one of them needs an explicit
        # ``extra_renames`` entry in the owning library's recipe to
        # pin the source to a distinct target name, or the target's
        # ``prefixes:`` map in ``codegen/targets/<target>.yaml`` must move
        # the conflicting letter to an unused one.
        targets: dict[str, list[str]] = {}
        for src, tgt in rename.items():
            targets.setdefault(tgt, []).append(src)
        collisions = sorted(
            (tgt, sorted(srcs)) for tgt, srcs in targets.items()
            if tgt in self.independent
        )
        # Flag only CROSS-family collisions. Co-family pairs (SGEMM,
        # DGEMM) both rename to QGEMM by design — the same family
        # converges. The dangerous case is two DISTINCT families
        # (different upstream patterns) producing the same target,
        # which would silently merge two unrelated routines.
        cross_family: list[tuple[str, list[str]]] = []
        for tgt, srcs in targets.items():
            if len(srcs) <= 1:
                continue
            # Two sources are in the same family iff their classifier
            # entry is the same Python object. id() gives a stable
            # cross-family discriminator without needing the family
            # dataclass to be hashable.
            fam_ids = {id(self._symbol_to_family.get(s.upper())) for s in srcs}
            if len(fam_ids) > 1:
                cross_family.append((tgt, sorted(srcs)))
        cross_family.sort()
        if cross_family:
            lines = [
                f"  {tgt!r}: produced by {srcs}"
                for tgt, srcs in cross_family
            ]
            raise RuntimeError(
                "prefix-classifier cross-family rename collisions:\n"
                + "\n".join(lines)
                + "\n  Two distinct upstream families rename to the "
                "same target — the migrator would silently merge "
                "them, leading to nondeterministic link-time symbol "
                "selection and (typically) a SEGV from a calling "
                "convention mismatch. Add an ``extra_renames:`` "
                "entry in the owning library's recipe (codegen/recipes/"
                "<lib>.yaml) that pins one of the conflicting sources "
                "to a distinct target name, or remap the conflicting "
                "letter in ``prefixes:`` of codegen/targets/<target>.yaml to "
                "an unused letter (W and U are unused across the "
                "entire BLAS / LAPACK / ScaLAPACK / MUMPS source "
                "corpus)."
            )
        if collisions:
            import sys
            lines = [
                f"  {tgt!r}: rename target of {srcs}, "
                f"but also exists as an orphaned (un-renamed) symbol"
                for tgt, srcs in collisions
            ]
            print(
                "WARNING: prefix-classifier rename-map collisions:\n"
                + "\n".join(lines)
                + "\n  Orphan symbols and rename targets share names; the "
                "linker may silently pick the wrong definition. If the "
                "orphan is the intended implementation (e.g. an "
                "upstream-provided helper), this is benign. Otherwise, "
                "investigate the family-discovery picker or change the "
                "target prefix to avoid the overlap.",
                file=sys.stderr,
            )

        return rename


def _build_tagged_pattern(sym: str, positions: tuple[int, ...]) -> str:
    """Build a tagged pattern string for sym with the given positions replaced.

    Each replaced char becomes a 2-character marker ``#R`` or ``#C``
    which encodes the slot's type. The returned string is thus longer
    than the original symbol by (#positions) characters, and two
    symbols sharing this string necessarily share both slot positions
    AND per-position type tags.
    """
    out: list[str] = []
    pos_set = set(positions)
    for i, ch in enumerate(sym):
        if i in pos_set:
            out.append('#' + CHAR_TYPE[ch])
        else:
            out.append(ch)
    return ''.join(out)


def _find_families_single_pass(
    symbols: set[str],
    replacement_chars: set[str],
) -> dict[tuple[str, tuple[int, ...]], dict[str, tuple[str, ...]]]:
    """Run one pass of pattern discovery.

    Returns { (tagged_pattern, positions) → { symbol: type_tags } }.
    ``tagged_pattern`` encodes per-position R/C tags, so two entries
    with the same key are guaranteed to match slot-for-slot.
    """
    groups: dict[tuple[str, tuple[int, ...]], dict[str, tuple[str, ...]]] = {}

    for sym in symbols:
        candidate_positions = [
            i for i, ch in enumerate(sym) if ch in replacement_chars
        ]
        if not candidate_positions:
            continue

        # Only enumerate subsets that contain at most one real and at
        # most one complex slot — i.e. valid tag-combinations are
        # (R), (C), (R,C). Multiple R-slots or multiple C-slots in the
        # same pattern let unrelated routines masquerade as siblings
        # (e.g. DSPTRS vs DSPTRD getting grouped on their second R
        # slot alone), so we forbid them.
        r_positions = [p for p in candidate_positions
                       if CHAR_TYPE[sym[p]] == 'R']
        c_positions = [p for p in candidate_positions
                       if CHAR_TYPE[sym[p]] == 'C']

        subsets: list[tuple[int, ...]] = []
        for rp in [None] + r_positions:
            for cp in [None] + c_positions:
                combo = tuple(sorted(
                    p for p in (rp, cp) if p is not None
                ))
                if combo:
                    subsets.append(combo)

        # Also allow doubled same-type prefix at positions (0, 1).
        # ScaLAPACK uses ZZDOTC/CCDOTC (wrappers around ZDOTC/CDOTC)
        # where both leading characters carry the precision tag.
        # The single-position patterns can't match across halves
        # (e.g. #CZDOTC ≠ #CCDOTC), so the 2-position pattern
        # #C#CDOTC is the only way to form the family.  This is
        # safe because the overlap requirement (both pass-1 and
        # pass-2 members must exist) prevents false groupings —
        # e.g. BLACS's CC* routines have no ZZ* counterparts.
        if len(c_positions) >= 2 and 0 in c_positions and 1 in c_positions:
            subsets.append((0, 1))
        if len(r_positions) >= 2 and 0 in r_positions and 1 in r_positions:
            subsets.append((0, 1))

        for positions in subsets:
            pattern = _build_tagged_pattern(sym, positions)
            tags = tuple(CHAR_TYPE[sym[p]] for p in positions)
            key = (pattern, positions)
            groups.setdefault(key, {})[sym] = tags

    return groups


def classify_symbols(symbols: set[str]) -> SymbolClassification:
    """Analyze a set of symbols to discover precision families.

    A family is a (tagged_pattern, positions) group whose members come
    from BOTH pass 1 (S/C-sourced) AND pass 2 (D/Z-sourced): the two
    halves form the pairs that the user's algorithm calls for.
    """
    upper_symbols = {s.upper() for s in symbols}

    single_groups = _find_families_single_pass(upper_symbols, SINGLE_CHARS)
    double_groups = _find_families_single_pass(upper_symbols, DOUBLE_CHARS)

    # Keep only groups with at least one pass-1 member AND one pass-2
    # member — that's the "overlap/pair" requirement.
    candidates: dict[tuple[str, tuple[int, ...]],
                     tuple[dict[str, tuple[str, ...]],
                           dict[str, tuple[str, ...]]]] = {}
    for key, single_members in single_groups.items():
        double_members = double_groups.get(key)
        if double_members:
            candidates[key] = (single_members, double_members)

    # Sort: most combined members first, then fewest positions, then
    # prefer slots at low indices (LAPACK/BLAS convention: precision
    # char is a prefix). The min-position tiebreaker is what keeps
    # DSYTRS/DSYTRD from being grouped as a sibling pair via their
    # trailing S/D — the competing families (#RSYTRS@0, #RSYTRD@0)
    # win because their slot sits at position 0. The final
    # (positions, pattern) tiebreaker ensures deterministic output.
    def _size(item):
        (pattern, positions) = item[0]
        (single, double) = item[1]
        return (
            -(len(single) + len(double)),
            len(positions),
            min(positions) if positions else 0,
            positions,
            pattern,
        )

    sorted_candidates = sorted(candidates.items(), key=_size)

    assigned: set[str] = set()
    result = SymbolClassification()

    for (pattern, positions), (single, double) in sorted_candidates:
        # Skip members already claimed by a larger/tighter family
        unassigned_single = {s: t for s, t in single.items()
                             if s not in assigned}
        unassigned_double = {s: t for s, t in double.items()
                             if s not in assigned}
        # Need the pair requirement to still hold after de-dup
        if not unassigned_single or not unassigned_double:
            continue

        # Tags are identical for every member (the key encodes them)
        type_tags = next(iter(unassigned_double.values()))
        slots = [
            PrecisionSlot(position=pos, type_tag=tag)
            for pos, tag in zip(positions, type_tags)
        ]

        all_members = {**unassigned_single, **unassigned_double}
        member_dict: dict[str, str] = {}
        for sym in sorted(all_members):
            prec_key = ''.join(sym[p] for p in positions)
            member_dict[prec_key] = sym

        family = PrecisionFamily(
            pattern=pattern,
            members=member_dict,
            slots=slots,
        )
        result.families.append(family)

        for sym in all_members:
            assigned.add(sym)
            result._symbol_to_family[sym] = family

    result.independent = upper_symbols - assigned
    return result


# Convenience wrappers

def build_rename_map(symbols: set[str], target_mode) -> dict[str, str]:
    """Build old_name → new_name mapping."""
    classification = classify_symbols(symbols)
    return classification.build_rename_map(target_mode)


def target_prefix(target_mode, is_complex: bool) -> str:
    """Return the target prefix character for a given mode and type."""
    pmap = target_mode.prefix_map
    return pmap['C'] if is_complex else pmap['R']
