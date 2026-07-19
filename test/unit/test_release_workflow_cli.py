"""Guard the release workflow's migrator invocations against the real CLI.

The convergence-report step once called ``python -m migrator converge`` — a
subcommand that had been renamed to ``diverge`` — and passed an ``output_dir``
positional the new subcommand no longer accepts. The step's ``|| true`` masked
the argparse error, so every release quietly rendered a usage message instead
of a divergence report.

These tests extract every ``python -m migrator ...`` command line from
``.github/workflows/release.yml`` and parse it with the same
``argparse`` parser the CLI builds, so a stale subcommand name or an
argument the subcommand rejects fails here instead of silently in CI.
"""

import contextlib
import io
import re
import shlex
from pathlib import Path

import pytest

from migrator.__main__ import build_parser

_REPO_ROOT = Path(__file__).resolve().parents[2]
_RELEASE_YML = _REPO_ROOT / '.github' / 'workflows' / 'release.yml'

# Shell control operators / redirections that terminate a simple command.
_STOP_TOKENS = {'||', '&&', ';', '|', '&'}


def _placeholder_expand(text):
    """Replace shell / GitHub-Actions interpolations with a literal token so
    the command's *structure* (subcommand + flags + positionals) can be parsed
    without knowing the runtime values."""
    text = re.sub(r'\$\{\{[^}]*\}\}', 'X', text)      # ${{ matrix.target }}
    text = re.sub(r'\$\{[A-Za-z_]\w*\}', 'X', text)   # ${VAR}
    text = re.sub(r'\$[A-Za-z_]\w*', 'X', text)       # $VAR
    return text


def _extract_migrator_commands(yml_text):
    """Yield the argument list (subcommand + args) of every
    ``python -m migrator ...`` command in the workflow, honouring backslash
    line continuations and stopping at shell control / redirection tokens."""
    lines = yml_text.splitlines()
    commands = []
    for i, line in enumerate(lines):
        if 'python -m migrator' not in line:
            continue
        parts = []
        j = i
        while j < len(lines):
            cur = lines[j].rstrip()
            if cur.endswith('\\'):
                parts.append(cur[:-1])
                j += 1
                continue
            parts.append(cur)
            break
        cmd = ' '.join(parts)
        # Keep only what follows the ``migrator`` module token.
        rest = cmd.split('migrator', 1)[1]
        tokens = shlex.split(_placeholder_expand(rest))
        trimmed = []
        for tok in tokens:
            if tok in _STOP_TOKENS or '>' in tok or '<' in tok:
                break
            trimmed.append(tok)
        if trimmed:
            commands.append(trimmed)
    return commands


# Guard the read so a missing release.yml fails test_release_yml_is_present
# with its intended message instead of erroring the whole module at import.
_COMMANDS = (
    _extract_migrator_commands(_RELEASE_YML.read_text())
    if _RELEASE_YML.is_file()
    else []
)


def test_release_yml_is_present():
    assert _RELEASE_YML.is_file(), _RELEASE_YML


def test_workflow_has_migrator_commands():
    # If extraction silently found nothing the per-command tests below would
    # vacuously pass, so assert the workflow really does drive the migrator.
    subs = {c[0] for c in _COMMANDS}
    assert 'stage' in subs
    assert 'diverge' in subs


@pytest.mark.parametrize('argv', _COMMANDS, ids=lambda a: ' '.join(a))
def test_workflow_command_parses(argv):
    """Every workflow migrator call names a real subcommand whose parser
    accepts the exact arguments the workflow passes. An unknown subcommand
    also fails here: the subparsers are required, so argparse rejects it via
    SystemExit with its own "invalid choice" message on the captured stderr."""
    parser = build_parser()
    with contextlib.redirect_stderr(io.StringIO()) as err:
        try:
            parser.parse_args(argv)
        except SystemExit:
            pytest.fail(
                f'release.yml command `migrator {" ".join(argv)}` is rejected '
                f'by the CLI parser:\n{err.getvalue().strip()}'
            )
