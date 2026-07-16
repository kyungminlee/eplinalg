"""Make ``migrator`` importable from tests."""

import sys
from pathlib import Path

_CODEGEN = Path(__file__).resolve().parents[2] / 'codegen'
if str(_CODEGEN) not in sys.path:
    sys.path.insert(0, str(_CODEGEN))
