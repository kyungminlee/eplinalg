# eplinalg documentation

- [user/](user/index.md) — using the produced libraries: install, link, call. Includes the [API reference](user/api/index.md).
- [dev/](dev/index.md) — developing eplinalg: configure, build, codegen, test, debug, release.
- [upstream-bugs/](upstream-bugs/README.md) — bugs in the vendored Netlib sources that the migrator works around without editing `extern/`.
- [archive/](archive/) — historical surveys and timestamped reports (not rendered).

Test-suite map: [../test/integration/README.md](../test/integration/README.md).

## Building the HTML site

Built with **Sphinx + MyST** and the **Furo** theme; published to GitHub
Pages by `.github/workflows/docs.yml` on pushes to `main`.

One-time setup (via [uv](https://docs.astral.sh/uv/)):

```sh
uv venv doc/.venv
uv pip install --python doc/.venv -r doc/requirements.txt
```

Then:

```sh
./doc/build.sh
```

Output lands in `doc/_build/html/index.html`. `conf.py` is **generated**
from `conf.py.in` by `build.sh` (the version string comes from the root
`VERSION` file) and is git-ignored — edit the `.in` file, not the
generated one.
