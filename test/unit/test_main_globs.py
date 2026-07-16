"""Tests for ``_collect_source_files`` — the source-file discovery helper used
by ``cmd_stage`` to enumerate migrated outputs into the CMake manifest.

The helper must include capital-F preprocessed Fortran (``.F``) so that
LAPACK 2-stage routines (``iparam2stage.F``, ``dsytrd_sb2st.F``,
``zhetrd_hb2st.F``) reach the manifest. ``Path.glob`` is case-sensitive on
Linux, so ``*.f`` does **not** match ``*.F``; the helper must glob both.
"""

from migrator.staging import _collect_source_files


def test_collect_fortran_picks_up_capital_F(tmp_path):
    (tmp_path / 'foo.f').write_text('')
    (tmp_path / 'bar.F').write_text('')
    (tmp_path / 'baz.f90').write_text('')
    (tmp_path / 'qux.F90').write_text('')

    names = sorted(p.name for p in _collect_source_files(tmp_path, 'fortran'))
    assert names == ['bar.F', 'baz.f90', 'foo.f', 'qux.F90']


def test_collect_fortran_omits_c(tmp_path):
    (tmp_path / 'foo.f').write_text('')
    (tmp_path / 'bar.F').write_text('')
    (tmp_path / 'helper.c').write_text('')

    names = sorted(p.name for p in _collect_source_files(tmp_path, 'fortran'))
    assert names == ['bar.F', 'foo.f']


def test_collect_c_includes_fortran_helpers(tmp_path):
    """C recipes (PBLAS) carry ``copy_files`` Fortran helpers like
    ``pilaenv.f``; the helper must surface ``.c`` plus all Fortran extensions
    so the manifest captures both."""
    (tmp_path / 'main.c').write_text('')
    (tmp_path / 'helper.f').write_text('')
    (tmp_path / 'guarded.F').write_text('')
    (tmp_path / 'modern.f90').write_text('')

    names = sorted(p.name for p in _collect_source_files(tmp_path, 'c'))
    assert names == ['guarded.F', 'helper.f', 'main.c', 'modern.f90']


def test_collect_iparam2stage_F_present(tmp_path):
    """Direct regression for the LAPACK 2-stage manifest gap: iparam2stage.F
    must appear in the output set so linking does not fail with
    ``undefined reference to iparam2stage_``."""
    (tmp_path / 'iparam2stage.F').write_text('')
    (tmp_path / 'dsytrd_sb2st.F').write_text('')
    (tmp_path / 'zhetrd_hb2st.F').write_text('')
    (tmp_path / 'dgesv.f').write_text('')  # control: pre-existing routine

    names = sorted(p.name for p in _collect_source_files(tmp_path, 'fortran'))
    assert 'iparam2stage.F' in names
    assert 'dsytrd_sb2st.F' in names
    assert 'zhetrd_hb2st.F' in names
    assert 'dgesv.f' in names


def test_collect_dedupes_by_inode(tmp_path):
    """When a real file is reachable via two different glob patterns (e.g.
    on case-insensitive filesystems where ``foo.f`` and ``foo.F`` resolve to
    the same inode), the helper must yield it only once. Simulate via a
    hardlink so this test runs portably on Linux."""
    real = tmp_path / 'foo.f'
    real.write_text('')
    link = tmp_path / 'foo.F'
    try:
        link.hardlink_to(real)
    except OSError:
        # Some filesystems (e.g. tmpfs without hardlink support) refuse;
        # accept the test as a no-op rather than failing on env oddities.
        return

    files = _collect_source_files(tmp_path, 'fortran')
    assert len(files) == 1
