"""Shared CLI helpers used across the migrator subcommands.

Small argument-to-object adapters (target-mode construction, parser
selection) that several ``cmd_*`` handlers and ``staging`` depend on.
Kept in their own module so both ``__main__`` and ``staging`` can import
them without a circular dependency. Extracted verbatim from ``__main__``.
"""
from .target_mode import load_target


def _get_target_mode(args):
    """Construct TargetMode based on CLI arguments."""
    target_str = getattr(args, 'target', None) or 'kind16'
    return load_target(target_str)


def _parser_args(args):
    """Extract parser/parser_cmd from CLI args."""
    parser = getattr(args, 'parser', None)
    parser_cmd = getattr(args, 'parser_cmd', None)
    return parser, parser_cmd
