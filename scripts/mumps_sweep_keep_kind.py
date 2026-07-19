"""Sweep every s/d and c/z pair under extern/MUMPS_5.9.0/src/ and
classify DOUBLE PRECISION declarations as keep-kind or promote.

Keep-kind rule (paired files only): a DP declaration in a d*/z* file is
keep-kind iff the same declaration (whitespace-normalized) also appears in
its s*/c* partner. Otherwise it is a working-precision DP to be promoted.

Shared files (un-prefixed, e.g. tools_common.F, mumps_load.F) are NOT in
the keep-kind manifest — they are copied through unchanged via copy_files
in the recipe. Every floating-point declaration in a shared file is
arithmetic-agnostic by construction (a shared file is linked into all four
arithmetic builds and therefore cannot depend on working precision), so no
per-line protection is needed: the file-level copy does the whole job.

This script also reports the list of shared files it found, so they can be
pasted into the recipe's copy_files section.

Usage:
    python mumps_sweep_keep_kind.py [--manifest <path>] [--no-write] [--verbose]

By default, the manifest is written to
    codegen/recipes/mumps/keep-kind.manifest
Re-run this script whenever the MUMPS source tree is updated.
"""

import argparse
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "extern" / "MUMPS_5.9.0" / "src"
INCLUDE = REPO / "extern" / "MUMPS_5.9.0" / "include"
DEFAULT_MANIFEST = REPO / "codegen" / "recipes" / "mumps" / "keep-kind.manifest"
DP_RE = re.compile(r"^\s*DOUBLE\s+PRECISION\b", re.IGNORECASE)
PAIRS = [("s", "d"), ("c", "z")]


def normalize(line: str) -> str:
    body = line.rstrip("\n")[:72]
    return re.sub(r"\s+", " ", body).strip().lower()


def dp_lines(path: Path) -> list[tuple[int, str, str]]:
    """Return [(lineno, raw, normalized)] for every DP declaration."""
    out = []
    for i, line in enumerate(path.read_text().splitlines(), start=1):
        if DP_RE.match(line):
            out.append((i, line, normalize(line)))
    return out


def classify_pair(lo: Path, hi: Path):
    """Classify DP decls in the (lo, hi) pair. A declaration is keep-kind
    iff the same normalized text appears in BOTH halves — in that case
    both halves must be protected from promotion to avoid cross-half
    divergence (the hi stays DP; an unprotected lo would promote and
    produce a different output).

    Returns
    -------
    keep_hi : list[(lineno, raw)]   lines in ``hi`` to protect
    keep_lo : list[(lineno, raw)]   lines in ``lo`` to protect
    promote_hi : list[(lineno, raw)] working-precision DPs in ``hi``
    """
    lo_by_norm: dict[str, list[tuple[int, str]]] = {}
    for ln, raw, norm in dp_lines(lo):
        lo_by_norm.setdefault(norm, []).append((ln, raw))

    keep_hi, promote_hi = [], []
    keep_lo_pairs: set[tuple[int, str]] = set()
    for ln, raw, norm in dp_lines(hi):
        if norm in lo_by_norm:
            keep_hi.append((ln, raw))
            for entry in lo_by_norm[norm]:
                keep_lo_pairs.add(entry)
        else:
            promote_hi.append((ln, raw))
    keep_lo = sorted(keep_lo_pairs)
    return keep_hi, keep_lo, promote_hi


def sweep_pairs(directory: Path, glob_suffix: str, totals, per_file,
                manifest_lines, paired_files):
    """Classify every lo/hi pair matching ``{hi_pref}{glob_suffix}`` in
    ``directory`` — rule is "DP in both sides → keep"."""
    for lo_pref, hi_pref in PAIRS:
        for hi_path in sorted(directory.glob(f"{hi_pref}{glob_suffix}")):
            lo_path = hi_path.with_name(lo_pref + hi_path.name[1:])
            if not lo_path.exists():
                continue  # not a real arithmetic variant; handled as shared
            paired_files.update({lo_path, hi_path})
            if not dp_lines(hi_path):
                continue
            keep_hi, keep_lo, promote_hi = classify_pair(lo_path, hi_path)
            totals["pairs"] += 1
            totals["keep"] += len(keep_hi) + len(keep_lo)
            totals["promote"] += len(promote_hi)
            per_file.append((hi_path.name, lo_path.name,
                             len(keep_hi) + len(keep_lo), len(promote_hi)))
            for ln, _ in keep_hi:
                manifest_lines.append(
                    f"{hi_path.relative_to(SRC.parents[2])}:{ln}")
            for ln, _ in keep_lo:
                manifest_lines.append(
                    f"{lo_path.relative_to(SRC.parents[2])}:{ln}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                    help=f"keep-kind manifest output "
                         f"(default: {DEFAULT_MANIFEST.relative_to(REPO)})")
    ap.add_argument("--no-write", action="store_true",
                    help="report only; don't write the manifest")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    totals = {"pairs": 0, "keep": 0, "promote": 0}
    per_file = []
    manifest_lines: list[str] = []
    paired_files: set[Path] = set()

    # Pass 1: paired s/d and c/z files.
    sweep_pairs(SRC, "*.F", totals, per_file, manifest_lines, paired_files)

    # Pass 1b: arithmetic-keyed headers under include/ (dmumps_struc.h
    # and siblings). Same rule — a DP decl in the d/z half is keep-kind
    # iff the identical decl appears in its s/c partner. Without this
    # pass, genuinely DP-stable struct fields (MEM_SUBTREE, COST_TRAV,
    # etc.) would be promoted by the migrator, breaking the interface
    # to DP-only shared routines in mumps_load.F etc.
    sweep_pairs(INCLUDE, "mumps_struc.h", totals, per_file, manifest_lines,
                paired_files)

    # Pass 2: discover shared files (for the copy_files report).
    shared_files = sorted(
        p.name for p in SRC.glob("*.F") if p not in paired_files
    )

    print(f"pairs classified:  {totals['pairs']}")
    print(f"keep-kind (paired): {totals['keep']}")
    print(f"promote (paired):   {totals['promote']}")
    print(f"shared files:       {len(shared_files)} "
          f"(should go in recipe `copy_files`)")

    if args.verbose:
        print("\nper-file keep-kind (hi, lo, keep, promote):")
        for row in per_file:
            print(f"  {row[0]:40s} {row[1]:40s} "
                  f"keep={row[2]:4d}  promote={row[3]:4d}")
        print("\nshared files for copy_files:")
        for name in shared_files:
            print(f"  - {Path(name).stem.upper()}")

    if not args.no_write:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text("\n".join(manifest_lines) + "\n")
        print(f"\nmanifest written: {args.manifest} "
              f"({len(manifest_lines)} entries)")


if __name__ == "__main__":
    main()
