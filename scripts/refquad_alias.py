#!/usr/bin/env python3
"""Convert a ref_quad_*.f90 module to dispatch via a `*_quad`-renamed
external symbol.

Background: refblas_quad / reflapack_quad expose Fortran symbols
(``dgemm_``, ``dasum_``, …) that collide with the std blas / lapack
archive's same-named symbols when both are linked into the same test
executable (kind4/kind8 baselines). The post-build step in
tests/<lib>/refblas (or reflapack) renames refblas_quad's exports to
``*_quad_`` via objcopy. This script rewrites ref_quad_*.f90 so:

  - Every `interface … subroutine NAME(args) … end subroutine`
    inside the module's `interface` block is renamed to NAME_quad in
    the interface declaration. The body (intent declarations, dummy
    types) is unchanged.

  - A `contains` section is appended after the `end interface`. For
    each renamed routine, a thin module procedure is generated under
    the *original* name. It re-declares the same dummy args (so test
    files calling `dasum(...)` see the right signature via `use
    ref_quad_blas, only: dasum`) and forwards to the `_quad` external.

After this transform:
  - Test files keep `use ref_quad_blas, only: dasum` — `dasum` resolves
    to a module procedure (`__ref_quad_blas_MOD_dasum`) — no link
    collision with the std archive's `dasum_`.
  - The module procedure internally calls `dasum_quad`, which gfortran
    emits as `dasum_quad_` and the linker resolves from the renamed
    refblas_quad archive — quad-precision computation as before.

Usage:
    python scripts/refquad_alias.py <input.f90> <output.f90>
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


SIG_START_RE = re.compile(
    r"""^(?P<indent>\s+)
        (?P<type_prefix>
          (?:logical|integer|real|complex|character|
             double\s+precision|double\s+complex)
          (?:\s*\([^)]*\))?
          \s+)?
        (?P<attr_prefix>
          (?:(?:recursive|pure|impure|elemental)\s+)+)?
        (?P<kind>function|subroutine)
        \s+(?P<name>[a-zA-Z_][a-zA-Z_0-9]*)
        \s*\(""", re.VERBOSE)


def join_continuations(lines: list[str], start: int) -> tuple[str, int]:
    """Join an `&`-continued logical line starting at `lines[start]`.

    Returns (joined_line, next_index). The joined line strips trailing
    `&` markers but preserves everything else verbatim. Blank/comment
    continuation lines are accepted (stripped).
    """
    parts = []
    i = start
    while i < len(lines):
        raw = lines[i]
        stripped = raw.rstrip()
        if stripped.endswith('&'):
            parts.append(stripped[:-1])
            i += 1
            continue
        parts.append(stripped)
        i += 1
        break
    return ' '.join(p.strip() for p in parts), i


def extract_args(signature: str) -> list[str]:
    """Pull the argument list from a joined `subroutine NAME(...)` or
    `function NAME(...) result(...)` signature."""
    # First parenthesized group after the routine name.
    open_paren = signature.find('(')
    if open_paren < 0:
        return []
    depth = 1
    j = open_paren + 1
    while j < len(signature) and depth > 0:
        if signature[j] == '(':
            depth += 1
        elif signature[j] == ')':
            depth -= 1
        j += 1
    inner = signature[open_paren + 1:j - 1]
    if not inner.strip():
        return []
    return [a.strip() for a in inner.split(',') if a.strip()]


def find_result_var(signature: str) -> str | None:
    """For a function, extract the result variable name from
    `... result(r)`; default 'r' if no result clause."""
    m = re.search(r'result\s*\(\s*([a-zA-Z_][a-zA-Z_0-9]*)\s*\)',
                  signature, re.IGNORECASE)
    return m.group(1) if m else None


def split_module(text: str) -> tuple[str, str, str]:
    """Split into (preamble, interface_block, postamble) by locating
    the outermost `interface` … `end interface`. Tracks nesting depth
    so callback-typed routines (e.g. `dgees(jobvs, sort, select, ...)`,
    whose body contains a nested `interface` block declaring SELECT)
    don't close the outer block prematurely.
    """
    lines = text.splitlines(keepends=True)
    iface_open = None
    iface_close = None
    depth = 0
    for i, raw in enumerate(lines):
        s = raw.strip().lower()
        if s == 'interface' or s.startswith('interface '):
            if depth == 0 and iface_open is None:
                iface_open = i
            depth += 1
            continue
        if s == 'end interface' or s.startswith('end interface'):
            depth -= 1
            if depth == 0 and iface_open is not None:
                iface_close = i
                break
    if iface_open is None or iface_close is None:
        raise ValueError('no top-level interface block found')
    pre = ''.join(lines[:iface_open + 1])
    block = ''.join(lines[iface_open + 1:iface_close])
    post = ''.join(lines[iface_close:])
    return pre, block, post


def parse_routines(block: str) -> list[dict]:
    """Walk the interface block, returning a list of per-routine dicts:
        { 'kind': 'function'|'subroutine',
          'name': str,
          'sig_lines': [line, ...],     # original lines covering signature (may be multi-line)
          'body_lines': [line, ...],    # body declarations between sig and end
          'end_line': line,
          'sig_joined': str,
          'args': [str, ...],
          'result_var': str | None,     # for functions
          'block_text': str }           # full original text including blanks/comments
    """
    lines = block.splitlines(keepends=True)
    routines: list[dict] = []
    i = 0
    leading_buffer: list[str] = []  # blank/comment lines to attach to next routine
    while i < len(lines):
        m = SIG_START_RE.match(lines[i])
        if not m:
            leading_buffer.append(lines[i])
            i += 1
            continue
        # Capture multi-line signature.
        sig_start = i
        joined, next_i = join_continuations(lines, i)
        sig_lines = lines[sig_start:next_i]
        # Walk body lines until "end subroutine NAME" or "end function NAME".
        kind = m.group('kind').lower()
        name = m.group('name')
        end_pat = re.compile(rf'^\s*end\s+{kind}\b(\s+{re.escape(name)})?\s*$',
                             re.IGNORECASE)
        body_start = next_i
        j = body_start
        end_line = None
        while j < len(lines):
            if end_pat.match(lines[j].rstrip()):
                end_line = lines[j]
                break
            j += 1
        if end_line is None:
            raise ValueError(f'no matching end {kind} for {name} at line {sig_start + 1}')
        body_lines = lines[body_start:j]
        leading_text = ''.join(leading_buffer)
        block_text = leading_text + ''.join(sig_lines) + ''.join(body_lines) + end_line
        leading_buffer = []
        args = extract_args(joined)
        result_var = find_result_var(joined) if kind == 'function' else None
        type_prefix = (m.group('type_prefix') or '').strip()
        routines.append(dict(
            kind=kind, name=name,
            sig_lines=sig_lines, body_lines=body_lines,
            end_line=end_line, sig_joined=joined,
            args=args, result_var=result_var,
            type_prefix=type_prefix,
            block_text=block_text,
            leading=leading_text,
        ))
        i = j + 1
    # Any trailing comment/blank lines stay outside the routine list;
    # we drop them — interface blocks generally end clean.
    return routines


def rename_in_signature_lines(sig_lines: list[str], old: str, new: str) -> list[str]:
    """Replace the routine name token on the first signature line."""
    out = list(sig_lines)
    out[0] = re.sub(rf'(\b)({re.escape(old)})(\b)', rf'\1{new}\3', out[0], count=1)
    return out


def rename_in_end_line(end_line: str, old: str, new: str) -> str:
    """Rename `end (function|subroutine) NAME` to use NAME+suffix.
    Tolerates the missing-name form (`end subroutine`) by matching the
    name boundary if present; otherwise leaves the line alone."""
    pat = re.compile(rf'(\bend\s+(?:function|subroutine)\b)(\s+){re.escape(old)}\b',
                     re.IGNORECASE)
    return pat.sub(rf'\1\g<2>{new}', end_line, count=1)


def render_renamed_block(r: dict) -> str:
    """Reproduce the routine sub-block in the interface, but with name
    suffixed `_quad`."""
    new_sig = rename_in_signature_lines(r['sig_lines'], r['name'],
                                        r['name'] + '_quad')
    new_end = rename_in_end_line(r['end_line'], r['name'],
                                 r['name'] + '_quad')
    return ''.join(new_sig) + ''.join(r['body_lines']) + new_end


def _format_signature(indent: str, head: str, args: list[str],
                      tail: str = '', max_width: int = 100) -> str:
    """Format a function/subroutine signature, breaking the argument
    list across multiple `&`-continued lines if the single-line form
    would exceed gfortran's free-form line limit (132). Lines are
    capped to ``max_width`` for readability — well under the limit."""
    single = f"{indent}{head}({', '.join(args)}){tail}\n"
    if len(single.rstrip('\n')) <= max_width or not args:
        return single
    cont_indent = indent + ' ' * (len(head) + 1)
    out = [f"{indent}{head}({args[0]}"]
    cur = out[0]
    for a in args[1:]:
        candidate = cur + f", {a}"
        if len(candidate) > max_width:
            out[-1] = cur + ', &\n'
            cur = cont_indent + a
            out.append(cur)
        else:
            cur = candidate
            out[-1] = cur
    out[-1] = cur + f"){tail}\n"
    return ''.join(out)


def render_wrapper(r: dict) -> str:
    """Generate a module-level procedure wrapping NAME_quad under the
    original NAME. Re-declares the same dummy types so the wrapper has
    a complete signature visible to test callers."""
    lines: list[str] = []
    indent_outer = '    '
    indent_inner = '        '
    args = r['args']
    if r['kind'] == 'function':
        type_prefix = r.get('type_prefix') or ''
        if type_prefix:
            # Typed function (e.g. `logical function disnan`) — no
            # result clause; the function name itself is the result.
            head = f"{type_prefix} function {r['name']}"
            lines.append(_format_signature(indent_outer, head, args))
            rv = r['name']
        else:
            rv = r['result_var'] or 'r'
            lines.append(_format_signature(indent_outer, f"function {r['name']}",
                                           args, f" result({rv})"))
    else:
        lines.append(_format_signature(indent_outer, f"subroutine {r['name']}",
                                       args))
    # Re-declare the body. The wrapper is a module procedure (host
    # has `use prec_kinds, only: ep`), so `import :: ep` at the
    # wrapper-top level is invalid. But nested `interface` blocks
    # (e.g. callback signatures inside dgees / dgeesx / dgges / zgees /
    # zgeesx) DO need their own `import :: ep` because interface bodies
    # have isolated scope. Track depth and only drop top-level imports.
    nest = 0
    for raw in r['body_lines']:
        s = raw.strip()
        if not s:
            continue
        sl = s.lower()
        if sl == 'interface' or sl.startswith('interface '):
            nest += 1
            lines.append(f"{indent_inner}{raw.lstrip()}")
            continue
        if sl == 'end interface' or sl.startswith('end interface'):
            nest -= 1
            lines.append(f"{indent_inner}{raw.lstrip()}")
            continue
        if nest == 0 and sl.startswith('import'):
            continue
        lines.append(f"{indent_inner}{raw.lstrip()}")
    # Forward call.
    if r['kind'] == 'function':
        lines.append(_format_signature(indent_inner,
                                       f"{rv} = {r['name']}_quad",
                                       args))
        lines.append(f"{indent_outer}end function {r['name']}\n")
    else:
        lines.append(_format_signature(indent_inner,
                                       f"call {r['name']}_quad",
                                       args))
        lines.append(f"{indent_outer}end subroutine {r['name']}\n")
    lines.append('\n')
    return ''.join(lines)


def transform(text: str) -> str:
    pre, block, post = split_module(text)
    routines = parse_routines(block)

    # Rebuild interface block: keep the original block but with each
    # routine's sub-block replaced by the renamed version. We rebuild
    # by walking the original lines and substituting per-routine
    # block_text.
    rebuilt_block = block
    for r in routines:
        rebuilt_block = rebuilt_block.replace(r['block_text'],
                                              r['leading'] + render_renamed_block(r),
                                              1)
        # block_text starts with r['leading'] (the blank/comment lines
        # preceding the routine), so the replacement must re-prepend it
        # — otherwise those separator lines would be silently dropped
        # from the rebuilt interface block.

    # Build contains section.
    contains = ['\ncontains\n\n']
    for r in routines:
        contains.append(render_wrapper(r))
    contains_text = ''.join(contains)

    # `post` starts with the `end interface` line and goes through
    # `end module`. Insert contains_text just before `end module`.
    end_mod_re = re.compile(r'(^\s*end\s+module\b.*$)', re.IGNORECASE | re.MULTILINE)
    m = end_mod_re.search(post)
    if not m:
        raise ValueError('no `end module` found after interface block')
    new_post = post[:m.start()] + contains_text + post[m.start():]

    return pre + rebuilt_block + new_post


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    src = Path(sys.argv[1])
    dst = Path(sys.argv[2])
    text = src.read_text()
    if re.search(r'^\s*contains\s*$', text, re.MULTILINE):
        print(f'{src} already has a `contains` section — already '
              f'transformed; refusing to re-apply', file=sys.stderr)
        return 1
    out = transform(text)
    dst.write_text(out)
    print(f'wrote {dst} ({len(out.splitlines())} lines)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
