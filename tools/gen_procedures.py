#!/usr/bin/env python3
"""Generate doc/output/procedures.md by walking the staged source trees.

Algorithm: anchor on the kind16 stage (whose prefixes — q/x/qx/iq/ix/xq —
don't overlap any other valid prefix). For each (family, stem) found in
kind16 via longest-prefix-first matching, construct the corresponding
name for every other target by direct string concatenation, and check
whether that file actually exists in that target's source tree. This
sidesteps the ambiguity of, say, treating Netlib `scopy.f` as either
's' + 'copy' or 'sc' + 'opy'.

Rows are dropped if the multifloats cell is missing — that filters out
universal utilities (`xerbla`, `xerbla_array`, `lsame`, …) and
upstream-only names (`dsdot`, `sdsdot`) that the migrator never re-emits.
"""
from pathlib import Path
import re

ROOT = Path("/home/kyungminlee/Code/fortran-migrator")
NETLIB = {
    "blas":      ROOT / "external/lapack-3.12.1/BLAS/SRC",
    "lapack":    ROOT / "external/lapack-3.12.1/SRC",
    "pblas":     ROOT / "external/scalapack-2.2.3/PBLAS/SRC",
    "scalapack": ROOT / "external/scalapack-2.2.3/SRC",
}

# Per-target precision-prefix maps. Families:
#   R  = real-prefix     (s/d/e/q/m)
#   C  = complex-prefix  (c/z/y/x/w)
#   RC = mixed real-from-complex (sc/dz/ey/qx/mw)  e.g. scnrm2 / dznrm2
#   CR = complex-vector with real scalar (cs/zd/ye/xq/wm)  e.g. csrot / zdscal
#   IR = integer-result, real-input (is/id/ie/iq/im)
#   IC = integer-result, complex-input (ic/iz/iy/ix/iw)
# Multifloats prefixes used to be t/v (and earlier dd/zz); the active
# scheme is m/w per ``targets/multifloats.yaml``.
PREFIXES = {
    "single":      {"R": "s",  "C": "c",  "RC": "sc", "CR": "cs", "IR": "is", "IC": "ic"},
    "double":      {"R": "d",  "C": "z",  "RC": "dz", "CR": "zd", "IR": "id", "IC": "iz"},
    "kind10":      {"R": "e",  "C": "y",  "RC": "ey", "CR": "ye", "IR": "ie", "IC": "iy"},
    "kind16":      {"R": "q",  "C": "x",  "RC": "qx", "CR": "xq", "IR": "iq", "IC": "ix"},
    "multifloats": {"R": "m",  "C": "w",  "RC": "mw", "CR": "wm", "IR": "im", "IC": "iw"},
}
TARGETS = list(PREFIXES.keys())
FAMILIES = ["R", "C", "RC", "CR", "IR", "IC"]
FAMILY_LABEL = {
    "R":  "Real-prefix family",
    "C":  "Complex-prefix family",
    "RC": "Mixed real-from-complex family",
    "CR": "Complex-vector with real-scalar family",
    "IR": "Integer-result, real-input",
    "IC": "Integer-result, complex-input",
}

# Longest-prefix-first ordering (per target) for unambiguous classification.
KIND16_BUCKETS = sorted(
    [(PREFIXES["kind16"][fam], fam) for fam in FAMILIES],
    key=lambda kv: (-len(kv[0]), kv[0]),
)


def list_basenames(dir_path: Path) -> set[str]:
    """Return the set of basenames (stripped of trailing '_' for PBLAS C
    wrappers and of `.f` / `.f90` / `.c` extensions) present in `dir_path`."""
    if not dir_path.is_dir():
        return set()
    out: set[str] = set()
    for f in dir_path.iterdir():
        m = re.match(r"^(.+?)(?:_)?\.(?:f|f90|c)$", f.name)
        if m:
            out.add(m.group(1))
    return out


def classify_kind16(filenames: set[str], p_prefix: str) -> set[tuple[str, str]]:
    """Map kind16 filenames to (family, stem) pairs."""
    out: set[tuple[str, str]] = set()
    for full in filenames:
        if not full.startswith(p_prefix):
            continue
        rest = full[len(p_prefix):]
        for pfx, fam in KIND16_BUCKETS:
            if rest.startswith(pfx):
                stem = rest[len(pfx):]
                if stem and re.match(r"^[a-z][a-z0-9_]*$", stem):
                    out.add((fam, stem))
                    break
    return out


def collect(library: str, p_prefix: str) -> dict[str, dict[str, dict[str, str | None]]]:
    """Return {family: {stem: {target: name_or_None}}}."""
    files = {tgt: list_basenames(
                NETLIB[library] if tgt in ("single", "double")
                else Path(f"/tmp/stg-{tgt}/{library}/src"))
             for tgt in TARGETS}
    rows: set[tuple[str, str]] = classify_kind16(files["kind16"], p_prefix)

    out: dict[str, dict[str, dict[str, str | None]]] = {fam: {} for fam in FAMILIES}
    for fam, stem in rows:
        cells: dict[str, str | None] = {}
        for tgt in TARGETS:
            name = f"{p_prefix}{PREFIXES[tgt][fam]}{stem}"
            cells[tgt] = name if name in files[tgt] else None
        # Drop universal utilities: keep only rows where the migrator
        # actually substituted across both kind16 and multifloats.
        if cells["multifloats"] is None:
            continue
        out[fam][stem] = cells
    return out


def write_table(out, library: str, by_fam: dict, p_prefix: str = "") -> None:
    out.write(f"\n## {library}\n\n")
    total = 0
    for fam in FAMILIES:
        stems = sorted(by_fam[fam])
        if not stems:
            continue
        out.write(f"\n### {library} — {FAMILY_LABEL[fam]} ({len(stems)} stems)\n\n")
        out.write("| stem | single | double | kind10 | kind16 | multifloats |\n")
        out.write("|------|--------|--------|--------|--------|-------------|\n")
        for stem in stems:
            cells = by_fam[fam][stem]
            row = [stem]
            for tgt in TARGETS:
                row.append(f"`{cells[tgt]}`" if cells[tgt] else "—")
            out.write("| " + " | ".join(row) + " |\n")
        total += len(stems)
    out.write(f"\n**{library} total: {total} unique (family, stem) rows.**\n")


def merge(a: dict, b: dict) -> dict:
    """Merge two {family: {stem: cells}} dicts (no overlap expected)."""
    out: dict[str, dict[str, dict[str, str | None]]] = {fam: {} for fam in FAMILIES}
    for fam in FAMILIES:
        out[fam].update(a.get(fam, {}))
        out[fam].update(b.get(fam, {}))
    return out


def main():
    out_path = ROOT / "doc/output/procedures.md"
    # Some libraries have multiple leading-letter conventions we want to
    # cover: ScaLAPACK has both `p<prec>...` distributed routines and
    # `b<prec>...` back-transformation helpers (bdlaapp, bdlaexc, …).
    libs: list[tuple[str, str, list[str]]] = [
        ("BLAS",      "blas",      [""]),
        ("LAPACK",    "lapack",    [""]),
        ("PBLAS",     "pblas",     ["p"]),
        ("ScaLAPACK", "scalapack", ["p", "b"]),
    ]
    with out_path.open("w") as f:
        f.write("# Routine cross-reference: Netlib → migrated targets\n\n")
        f.write("Per-routine name lookup across the upstream Netlib (single + double)\n")
        f.write("and the three migrator targets (kind10 / kind16 / multifloats).\n")
        f.write("Each row is a *stem* — the routine identity stripped of its\n")
        f.write("precision/family prefix; columns show the prefixed name as it\n")
        f.write("appears in source / object files. `—` means that target does not\n")
        f.write("expose this stem (the migrator did not emit it, or the upstream\n")
        f.write("has no counterpart in that precision).\n\n")
        f.write("Generated by `tools/gen_procedures.py` from staged source trees at\n")
        f.write("`/tmp/stg-{kind10,kind16,multifloats}/<lib>/src/` and the Netlib\n")
        f.write("references at `external/lapack-3.12.1/{BLAS,LAPACK}/SRC` and\n")
        f.write("`external/scalapack-2.2.3/{PBLAS,SRC}`.\n\n")
        f.write("**Prefix map**\n\n")
        f.write("| family | single | double | kind10 | kind16 | multifloats |\n")
        f.write("|--------|--------|--------|--------|--------|-------------|\n")
        labels = [("R", "real (R)"), ("C", "complex (C)"),
                  ("RC", "real-from-complex (RC)"),
                  ("CR", "complex-with-real-scalar (CR)"),
                  ("IR", "int-of-real (IR)"), ("IC", "int-of-complex (IC)")]
        for fam, label in labels:
            row = [label]
            for tgt in TARGETS:
                row.append(PREFIXES[tgt][fam])
            f.write("| " + " | ".join(row) + " |\n")
        f.write("\nFor PBLAS and ScaLAPACK every name additionally carries a\n")
        f.write("leading `p` (e.g. `pdgemm` / `pqgemm` / `pmgemm`). ScaLAPACK\n")
        f.write("also has a small `b`-prefix family for back-transformation\n")
        f.write("helpers (`bdlaapp`, `bdlaexc`, `bdtrexc`); their cells\n")
        f.write("therefore start with `b` rather than `p`.\n\n")
        f.write("Out of scope for this index: PBLAS-internal helper libraries\n")
        f.write("(`ptzblas`, `pbblas`) and ScaLAPACK's C-side wrappers\n")
        f.write("(`scalapack_c`). These ship in their own staged subtrees.\n")
        for name, sub, p_list in libs:
            by_fam: dict = {fam: {} for fam in FAMILIES}
            for p in p_list:
                by_fam = merge(by_fam, collect(sub, p))
            write_table(f, name, by_fam)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
