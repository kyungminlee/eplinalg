#!/usr/bin/env bash
# Build the eplinalg documentation: Sphinx/MyST (Markdown) -> HTML.
#
# One-time setup:
#     uv venv doc/.venv
#     uv pip install --python doc/.venv -r doc/requirements.txt
#
# Usage:  ./doc/build.sh   (run from the repo root or from doc/)
set -euo pipefail

DOC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$DOC_DIR/.." && pwd)"
VENV="$DOC_DIR/.venv"
SPHINX="$VENV/bin/sphinx-build"

if [[ ! -x "$SPHINX" ]]; then
    echo "error: $SPHINX not found. Run:" >&2
    echo "    uv venv doc/.venv && uv pip install --python doc/.venv -r doc/requirements.txt" >&2
    exit 1
fi

echo ">> configure (VERSION -> conf.py)"
EP_VERSION="$(tr -d '[:space:]' < "$ROOT_DIR/VERSION")"
sed "s/@PROJECT_VERSION@/$EP_VERSION/g" "$DOC_DIR/conf.py.in" > "$DOC_DIR/conf.py"

echo ">> sphinx-build (-> HTML)"
"$SPHINX" -b html "$DOC_DIR" "$DOC_DIR/_build/html" "$@"

echo ">> done: $DOC_DIR/_build/html/index.html"
