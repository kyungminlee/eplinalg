"""Tests for recipe loading and schema validation."""

import textwrap

import pytest

from migrator.config import load_recipe


def _write_recipe(tmp_path, body: str):
    recipe_dir = tmp_path / 'recipes'
    recipe_dir.mkdir()
    src_dir = tmp_path / 'src'
    src_dir.mkdir()
    recipe = recipe_dir / 'r.yaml'
    recipe.write_text(textwrap.dedent(body))
    return recipe


def test_minimal_recipe_loads(tmp_path):
    recipe = _write_recipe(tmp_path, """
        library: blas
        language: fortran
        source_dir: src
    """)
    config = load_recipe(recipe, project_root=tmp_path)
    assert config.library == 'blas'
    assert config.language == 'fortran'
    assert config.source_dir == tmp_path / 'src'


def test_missing_required_key_raises(tmp_path):
    recipe = _write_recipe(tmp_path, """
        library: blas
        # language and source_dir intentionally missing
    """)
    with pytest.raises(KeyError, match='language'):
        load_recipe(recipe, project_root=tmp_path)


def test_unknown_keys_warn_on_stderr(tmp_path, capsys):
    recipe = _write_recipe(tmp_path, """
        library: blas
        language: fortran
        source_dir: src
        copy-files: [FOO]   # typo: dash vs underscore
        spurious_key: 42
    """)
    config = load_recipe(recipe, project_root=tmp_path)
    captured = capsys.readouterr()
    assert 'unknown recipe key' in captured.err
    # The typo'd key should be mentioned by name, so the user can find
    # and fix the offending line.
    assert 'copy-files' in captured.err
    assert 'spurious_key' in captured.err
    # The mistyped key cannot have populated copy_files; that's the
    # exact silent-default-empty bug this warning guards against.
    assert config.copy_files == set()


def test_empty_symbols_block_loads(tmp_path):
    recipe = _write_recipe(tmp_path, """
        library: blas
        language: fortran
        source_dir: src
        symbols:
    """)
    # Empty ``symbols:`` block must load without raising; the field
    # (an nm-based archive scan) was retired along with
    # ``library_path``, but YAML loader tolerance is still a
    # regression-worthy property.
    load_recipe(recipe, project_root=tmp_path)


def test_empty_prefix_block_loads(tmp_path):
    recipe = _write_recipe(tmp_path, """
        library: blas
        language: fortran
        source_dir: src
        prefix:
    """)
    # Empty ``prefix:`` block must load without raising; the field
    # itself was retired in v0.3.3, but YAML loader tolerance is
    # still a regression-worthy property.
    load_recipe(recipe, project_root=tmp_path)


def test_non_mapping_top_level_raises(tmp_path):
    recipe = _write_recipe(tmp_path, """
        - not_a_mapping
    """)
    with pytest.raises(ValueError, match='must be a mapping'):
        load_recipe(recipe, project_root=tmp_path)
