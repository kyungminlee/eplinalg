"""Stage upstream sources for a recipe and apply its patch list.

The prepared tree lives at ``<project_root>/build/staged-sources/<library>/``
and mirrors the contents of the recipe's ``source_dir``. Patches listed in
``recipe.patches`` are unified diffs (typically produced by ``git diff
--no-index``) applied with ``git apply`` from inside the staged tree.

Idempotency: a ``.prepared.stamp`` file inside the staged tree records when
preparation last completed and a fingerprint of the upstream source tree
(file names, sizes, mtimes). If the stamp is newer than every listed patch
file and the upstream fingerprint is unchanged, the stage is skipped.
Pass ``rebuild=True`` to force a clean re-stage.

This module is the input side of the pipeline reshape described in
``doc/archive/refactor-20260509.md``. Phase A wires it as a no-op-by-default CLI
command; subsequent phases route migration to read from the staged tree.
"""

from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path

from .config import RecipeConfig, load_recipe


STAMP_NAME = '.prepared.stamp'


def staged_root_for(project_root: Path, library: str) -> Path:
    """Where ``library``'s staged sources live."""
    return project_root / 'build' / 'staged-sources' / library


# ---------------------------------------------------------------------------
# Kind-normalization (post-patch, pre-migrate)
# ---------------------------------------------------------------------------

# Match floating-point literals with explicit E exponent (e.g. ``1.0E+0``)
# OR bare unsuffixed (``1.0``) OR with explicit kind4/8 trailers we don't
# rewrite (``1.0_8`` / ``1.0_wp``). The regex captures the mantissa and
# exponent groups so the substitution can rewrite to D-suffix form.
_E_LITERAL_RE = re.compile(r'(?<![\w.])(\d+\.\d*|\d*\.\d+)([Ee])([+-]?\d+)')
_BARE_LITERAL_RE = re.compile(
    r'(?<![\w.])(\d+\.\d*|\d*\.\d+)(?![DdEe\w]|_\d|_[A-Za-z])'
)
_INTEGER_PARAM_DECL_RE = re.compile(
    r'^\s*INTEGER\b.*\bPARAMETER\b.*::', re.IGNORECASE
)
_FORTRAN_OP_RE = re.compile(
    r'\.\s*(EQ|NE|LT|GT|LE|GE|AND|OR|NOT|TRUE|FALSE|EQV|NEQV)\s*\.',
    re.IGNORECASE,
)
_STRING_SPLIT_RE = re.compile(r"('(?:[^']|'')*'|\"(?:[^\"]|\"\")*\")")


def _normalize_kind8_literals_in_segment(seg: str) -> str:
    """Rewrite kind4-shaped literals to kind8 (D-suffix) form in
    ``seg``, which is a non-string-literal slice of one source line.
    Operator markers like ``.AND.`` are masked first so the dot-loose
    bare-literal regex doesn't consume their leading dot.
    """
    masked = _FORTRAN_OP_RE.sub(
        lambda m: '\x00' + m.group(1) + '\x00', seg,
    )
    # E-suffix → D-suffix, preserving sign.
    masked = _E_LITERAL_RE.sub(r'\1D\3', masked)
    # Bare → explicit D+0.
    masked = _BARE_LITERAL_RE.sub(r'\1D+0', masked)
    # Restore operator markers.
    return re.sub(r'\x00([A-Za-z]+)\x00', r'.\1.', masked)


def normalize_d_half_kinds(source: str, suffix: str) -> str:
    """Promote kind4-shaped literals to kind8 (D-suffix) form.

    Targeted at D/Z-half source files in the staged tree as a post-
    patch normalization step. The migrator's per-half rule (a) gate
    only promotes kind8 in D/Z halves — but upstream D/Z files routinely
    use bare ``1.0`` or `1.0E+0` for default-precision constants that
    Fortran auto-promotes via mixed-kind expression rules. After
    migration with rule (a), those tokens stay kind4 in the canonical
    D/Z output, while the corresponding S/C re-migration promotes them
    (because kind4 IS the source kind on the kind4 half), producing a
    spurious textual divergence.

    Normalizing kind4 literals to kind8 on the D/Z input *before*
    migration makes both halves' migrated output land on kind16/whatever
    target form uniformly. Comments and string literals are left alone.

    ``suffix`` is the source extension (``.f``, ``.F``, ``.f90``, etc.).
    Fixed-form (``.f`` / ``.F``) recognizes a comment marker (``C`` /
    ``*`` / ``!`` / ``c``) in column 1; free-form recognizes ``!``
    anywhere outside string literals.
    """
    is_free = suffix.lower() in ('.f90', '.f95')
    out: list[str] = []
    for line in source.splitlines(keepends=True):
        # Comment lines pass through unchanged.
        if not is_free and line[:1] in ('C', 'c', '*', '!'):
            out.append(line)
            continue
        if is_free and line.lstrip().startswith('!'):
            out.append(line)
            continue
        # Skip ``INTEGER, PARAMETER :: name = <fp_literal>`` — gfortran
        # tolerates the FP literal via integer coercion; rewriting to
        # D-suffix preserves that behavior but the migrator's
        # downstream literal-handler skips these too, so leaving them
        # untouched here keeps the round-trip clean.
        if _INTEGER_PARAM_DECL_RE.match(line):
            out.append(line)
            continue
        # Free-form: strip a trailing ``! …`` comment before normalizing
        # the code part, then reattach.
        if is_free:
            # Honor strings: only treat ``!`` outside strings as comment.
            parts = _STRING_SPLIT_RE.split(line)
            comment_started_at: int | None = None
            for idx in range(0, len(parts), 2):
                bang = parts[idx].find('!')
                if bang != -1:
                    parts[idx] = parts[idx][:bang]
                    comment_started_at = idx
                    break
            if comment_started_at is not None:
                # Discard everything after the comment marker on this line
                head = ''.join(parts[: comment_started_at + 1])
                # Find the original tail (including the ``!`` and trailing newline)
                head_len = len(head)
                tail = line[head_len:]
                code = head
            else:
                code = line
                tail = ''
        else:
            # Fixed-form: process whole line; inline ``!`` comments are
            # rare and the string-split below leaves them inside any
            # quoted region. Continuation char in column 6 is a single
            # non-blank char; the regex won't match it.
            code = line
            tail = ''
        parts = _STRING_SPLIT_RE.split(code)
        for idx in range(0, len(parts), 2):
            parts[idx] = _normalize_kind8_literals_in_segment(parts[idx])
        out.append(''.join(parts) + tail)
    return ''.join(out)


def _is_kind8_half_filename(name: str) -> bool:
    """Return True if the file's stem starts with d/z (single-letter)
    or pd/pz (ScaLAPACK two-letter). These are the kind8 half files
    whose literals should be normalized.
    """
    stem = Path(name).stem.lower()
    if not stem:
        return False
    if stem[0] == 'p' and len(stem) > 1 and stem[1] in 'dz':
        return True
    return stem[0] in 'dz'


def _normalize_staged_kinds(staged_root: Path) -> int:
    """Walk the staged tree and rewrite kind4 literals on D/Z-half
    files. Returns the number of files modified.
    """
    if not staged_root.is_dir():
        return 0
    modified = 0
    for path in staged_root.rglob('*'):
        if not path.is_file():
            continue
        suf = path.suffix.lower()
        if suf not in ('.f', '.f90', '.f95', '.for', '.h'):
            continue
        # Special-case .h: only Fortran-content headers (e.g. MUMPS
        # ``dmumps_struc.h``) — same as migrator policy for fortran
        # source dirs.
        if not _is_kind8_half_filename(path.name):
            continue
        try:
            text = path.read_text(errors='replace')
        except OSError:
            continue
        new_text = normalize_d_half_kinds(text, suf)
        if new_text != text:
            path.write_text(new_text)
            modified += 1
    return modified


def patch_dir_for(recipe_path: Path, library: str) -> Path:
    """Directory holding the recipe's patch files."""
    return recipe_path.parent / library / 'patches'


def _resolve_patches(recipe_path: Path,
                     config: RecipeConfig) -> list[Path]:
    pdir = patch_dir_for(recipe_path, config.library)
    resolved: list[Path] = []
    missing: list[str] = []
    for name in config.patches:
        p = pdir / name
        if not p.exists():
            missing.append(name)
        resolved.append(p)
    if missing:
        raise FileNotFoundError(
            f'{recipe_path.name}: patch file(s) not found in {pdir}: '
            f'{missing!r}'
        )
    return resolved


def _source_fingerprint(source_dir: Path) -> str:
    """Cheap content fingerprint of the upstream source tree.

    Hashes the sorted (relative-path, size, mtime_ns) triples of every
    file under ``source_dir``. The *name set* is part of the hash, so
    files added to or removed from upstream invalidate the cache even
    when their mtimes predate the stamp — ``git mv`` preserves mtime,
    so an mtime-only check misses files relocated into the tree.
    """
    h = hashlib.sha256()
    for p in sorted(source_dir.rglob('*')):
        if p.is_file():
            st = p.stat()
            h.update(
                f'{p.relative_to(source_dir)}\t'
                f'{st.st_size}\t{st.st_mtime_ns}\n'.encode()
            )
    return h.hexdigest()


def _is_fresh(stamp: Path, patches: list[Path], source_dir: Path) -> bool:
    if not stamp.exists():
        return False
    stamp_mtime = stamp.stat().st_mtime
    for p in patches:
        if p.stat().st_mtime > stamp_mtime:
            return False
    # Stamp body records the upstream fingerprint at stage time. An
    # empty body is a pre-fingerprint stamp — treat as stale once so
    # existing caches upgrade themselves.
    return stamp.read_text().strip() == _source_fingerprint(source_dir)


def run_prepare(recipe_path: Path,
                project_root: Path | None = None,
                *,
                rebuild: bool = False) -> Path:
    """Prepare the staged-sources tree for one recipe.

    Returns the staged root path.
    """
    if project_root is None:
        project_root = recipe_path.parent.parent.parent

    config = load_recipe(recipe_path, project_root)
    staged_root = staged_root_for(project_root, config.library)
    stamp = staged_root / STAMP_NAME

    patches = _resolve_patches(recipe_path, config)

    if not config.source_dir.exists():
        raise FileNotFoundError(
            f'{recipe_path.name}: source_dir does not exist: '
            f'{config.source_dir}'
        )

    if (not rebuild and staged_root.exists()
            and _is_fresh(stamp, patches, config.source_dir)):
        return staged_root

    if staged_root.exists():
        shutil.rmtree(staged_root)
    staged_root.parent.mkdir(parents=True, exist_ok=True)

    shutil.copytree(config.source_dir, staged_root)

    for patch in patches:
        # ``git apply`` is run from inside the staged tree so a patch
        # whose hunks reference top-level filenames (``a/dnrm2.f90``)
        # resolves directly. The patch path is made absolute first
        # because git apply re-resolves it against ``cwd=staged_root``,
        # not the caller's CWD.
        patch_abs = str(patch.resolve())
        # ``--check`` first: refusing on apply-conflict surfaces upstream
        # drift loud (vendor bumped a file the patch covers) instead of
        # silently producing a half-applied tree.
        subprocess.run(
            ['git', 'apply', '--whitespace=nowarn', '--check', patch_abs],
            cwd=staged_root, check=True,
        )
        subprocess.run(
            ['git', 'apply', '--whitespace=nowarn', patch_abs],
            cwd=staged_root, check=True,
        )

    # Post-patch kind normalization (rule (a) preconditioning). Rewrites
    # bare and E-suffix floating-point literals to D-suffix form on
    # every D/Z-half file in the staged tree, so that downstream
    # migration with rule (a) sees uniform kind8 literals on those
    # halves and the migrated output converges with the kind4 half's
    # re-migrated text. See ``normalize_d_half_kinds`` docstring for
    # the rationale and edge cases (string literals, INTEGER PARAMETER
    # idiom, comment markers).
    _normalize_staged_kinds(staged_root)

    stamp.write_text(_source_fingerprint(config.source_dir) + '\n')
    return staged_root


def _precision_sibling_set(filename: str,
                           upstream_root: Path | None = None) -> set[str] | None:
    """Return the co-family for ``filename`` (existing upstream files only).

    Generates candidate sibling names by precision-prefix swap, with the
    LAPACK orthogonal/unitary aliasing (``or`` ↔ ``un``) folded in:

      ``dorbdb3.f``  →  {``sorbdb3.f``, ``dorbdb3.f``, ``cunbdb3.f``, ``zunbdb3.f``}
      ``pdormrz.f``  →  {``psormrz.f``, ``pdormrz.f``, ``pcunmrz.f``, ``pzunmrz.f``}
      ``pdgeequ.f``  →  {``psgeequ.f``, ``pdgeequ.f``, ``pcgeequ.f``, ``pzgeequ.f``}

    When ``upstream_root`` is provided, candidates that don't correspond
    to a real upstream file are filtered out — so a candidate like
    ``dunbdb3.f`` (which doesn't exist; D-half is ``dorbdb3.f``) won't
    show up as a missing sibling.

    Returns ``None`` for files whose stem doesn't begin with a known
    precision prefix (e.g. ``utils.f``).
    """
    p = Path(filename)
    stem = p.stem.lower()
    suf = p.suffix

    # Detect prefix style. Try 2-letter (PS/PD/PC/PZ) first.
    for prefix_set, pfx_len in ((('ps', 'pd', 'pc', 'pz'), 2),
                                 (('s',  'd',  'c',  'z'),  1)):
        pfx = stem[:pfx_len]
        if pfx not in prefix_set or len(stem) <= pfx_len:
            continue
        tail = stem[pfx_len:]
        precision = pfx[-1]  # 's','d','c','z'
        is_real = precision in ('s', 'd')

        candidates: set[str] = set()
        for q in prefix_set:
            q_precision = q[-1]
            q_is_real = q_precision in ('s', 'd')
            # Apply or/un aliasing when crossing the real/complex boundary.
            new_tail = tail
            if is_real and not q_is_real and tail.startswith('or'):
                new_tail = 'un' + tail[2:]
            elif not is_real and q_is_real and tail.startswith('un'):
                new_tail = 'or' + tail[2:]
            candidates.add(f'{q}{new_tail}{suf}')

        if upstream_root is not None and upstream_root.is_dir():
            existing = {n for n in candidates
                        if (upstream_root / n).is_file()}
            return existing if existing else None
        return candidates
    return None


def _patch_touched_files(patch_path: Path) -> set[str]:
    """Return the set of filenames the patch's hunks add or modify."""
    touched: set[str] = set()
    for line in patch_path.read_text().splitlines():
        if line.startswith('+++ b/'):
            rel = line[len('+++ b/'):].strip().split('\t', 1)[0]
            touched.add(Path(rel).name)
    return touched


def verify_patches(recipe_path: Path,
                   project_root: Path | None = None) -> list[str]:
    """Symmetric-patch CI check.

    Aggregates every file touched by any patch in the recipe (excluding
    those whose patch is listed under ``asymmetric_patches:``). For
    each precision-prefixed file in the aggregate, requires all four
    co-family siblings to be touched by SOME patch — otherwise a fix
    in the D/Z half doesn't propagate to the S/C re-migration path.

    Sibling hunks may live in different patches; they often do
    (per-stem patches keep diffs reviewable).

    Returns a list of error strings; empty list means clean.
    """
    if project_root is None:
        project_root = recipe_path.parent.parent.parent
    config = load_recipe(recipe_path, project_root)
    pdir = patch_dir_for(recipe_path, config.library)
    if not pdir.is_dir():
        return []

    skip = set(config.asymmetric_patches) | set(config.one_sided_cleanup)
    aggregate: set[str] = set()
    file_to_patches: dict[str, list[str]] = {}
    errors: list[str] = []

    for patch_name in config.patches:
        if patch_name in skip:
            continue
        patch_path = pdir / patch_name
        if not patch_path.exists():
            errors.append(f'{patch_name}: not found in {pdir}')
            continue
        for f in _patch_touched_files(patch_path):
            aggregate.add(f)
            file_to_patches.setdefault(f, []).append(patch_name)

    upstream_root = config.source_dir
    for f in sorted(aggregate):
        sibs = _precision_sibling_set(f, upstream_root)
        if sibs is None:
            continue
        missing = sibs - aggregate
        if missing:
            errors.append(
                f'{f} (in {file_to_patches[f]}): precision-asymmetric — '
                f'co-family siblings {sorted(missing)} are not touched '
                f'by any patch. Add a patch that touches them, or list '
                f'the patch under asymmetric_patches: in the recipe.'
            )

    return errors


def prepare_recipe(recipe_path: Path,
                   project_root: Path | None = None,
                   *,
                   rebuild: bool = False) -> RecipeConfig:
    """Load a recipe with ``source_dir`` rewritten to its staged tree.

    Runs :func:`run_prepare` to ensure ``build/staged-sources/<library>/``
    exists and reflects the recipe's patch list, then loads the recipe
    and swaps ``config.source_dir`` to point at the staged tree. All
    downstream pipeline code (`scan_symbols`, `run_fortran_migration`,
    `run_c_migration`, …) reads through ``config.source_dir``, so this
    one swap routes the entire migration through the staged tree.
    """
    if project_root is None:
        project_root = recipe_path.parent.parent.parent
    staged_root = run_prepare(recipe_path, project_root, rebuild=rebuild)
    config = load_recipe(recipe_path, project_root)
    config.source_dir = staged_root
    return config
