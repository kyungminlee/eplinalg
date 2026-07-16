"""Tests for the post-migration ep_ symbol-privatization pass (task 44)."""

import textwrap

from migrator.config import load_recipe
from migrator.privatize import (
    build_c_map,
    build_fortran_map,
    privatize_c_text,
    privatize_fortran_text,
    privatize_tree,
)

NAMES = frozenset({
    # Fortran-callable (single trailing underscore)
    'blacs_gridinit_', 'descinit_', 'numroc_', 'pxerbla_',
    'sl_gridreshape_',
    # C-callable twins / C-only spellings
    'Cblacs_gridinit', 'BI_AuxBuff', 'PB_Ctop', 'getpbbuf',
    # decoration alias triple for the TOOLS reshape C function
    'reshape', 'reshape_', 'RESHAPE',
    # g77 double-underscore decoration alias (C-only)
    'sl_gridreshape__',
})


# ---------------------------------------------------------------- maps

def test_c_map_is_exact_and_case_sensitive():
    cmap = build_c_map(NAMES)
    assert cmap['Cblacs_gridinit'] == 'ep_Cblacs_gridinit'
    assert cmap['reshape'] == 'ep_reshape'
    assert cmap['RESHAPE'] == 'ep_RESHAPE'
    assert cmap['sl_gridreshape__'] == 'ep_sl_gridreshape__'
    assert 'cblacs_gridinit' not in cmap


def test_fortran_map_stems_and_exclusions():
    fmap = build_fortran_map(NAMES)
    # single-underscore entries enter under their upper-cased stem
    assert fmap['DESCINIT'] == 'EP_DESCINIT'
    assert fmap['SL_GRIDRESHAPE'] == 'EP_SL_GRIDRESHAPE'
    # RESHAPE enters too: the band solvers CALL RESHAPE (the REDIST C
    # routine), and no staged/upstream Fortran uses the intrinsic
    assert fmap['RESHAPE'] == 'EP_RESHAPE'
    # C-only names (no trailing underscore) never enter the Fortran map
    assert 'CBLACS_GRIDINIT' not in fmap
    assert 'GETPBBUF' not in fmap
    assert 'BI_AUXBUFF' not in fmap
    # double-underscore aliases are C decoration, not Fortran stems
    assert 'SL_GRIDRESHAPE_' not in fmap


# ------------------------------------------------------------- C text

def test_c_definitions_calls_and_globals_renamed():
    cmap = build_c_map(NAMES)
    src = textwrap.dedent("""\
        extern BLACBUFF BI_AuxBuff;
        F_VOID_FUNC blacs_gridinit_(Int *ConTxt, F_CHAR order, Int *nprow,
                                    Int *npcol)
        {
           void Cblacs_gridinit(Int *, char *, Int, Int);
           Cblacs_gridinit(ConTxt, F2C_CharTrans(order), Mpval(nprow),
                           Mpval(npcol));
        }
    """)
    out = privatize_c_text(src, cmap)
    assert 'ep_blacs_gridinit_(' in out
    assert 'void ep_Cblacs_gridinit(' in out
    assert 'ep_Cblacs_gridinit(ConTxt' in out
    assert 'ep_BI_AuxBuff;' in out
    # untouched: types, non-manifest identifiers
    assert 'BLACBUFF' in out
    assert 'F2C_CharTrans' in out


def test_c_rename_is_case_sensitive_and_whole_token():
    cmap = build_c_map(NAMES)
    src = 'void reshape(); void RESHAPE(); void Creshape(); int reshaped;\n'
    out = privatize_c_text(src, cmap)
    assert 'void ep_reshape();' in out
    assert 'void ep_RESHAPE();' in out
    # Creshape not in this test manifest; reshaped is a different token
    assert 'void Creshape();' in out
    assert 'int reshaped;' in out


def test_c_mangling_defines_renamed_but_includes_untouched():
    cmap = build_c_map(NAMES)
    src = textwrap.dedent("""\
        #include "reshape.h"
        # include <mpi.h>
        #define reshape_ RESHAPE
        #ifdef UpCase
        #define blacs_gridinit_ BLACS_GRIDINIT
        #endif
    """)
    out = privatize_c_text(src, cmap)
    assert '#include "reshape.h"' in out          # file names never change
    assert '# include <mpi.h>' in out
    assert '#define ep_reshape_ ep_RESHAPE' in out
    # LHS renamed; RHS (UpCase decoration, not in manifest) untouched
    assert '#define ep_blacs_gridinit_ BLACS_GRIDINIT' in out


def test_c_pass_is_idempotent():
    cmap = build_c_map(NAMES)
    src = 'Cblacs_gridinit(&ictxt, "R", 1, np);\n'
    once = privatize_c_text(src, cmap)
    assert privatize_c_text(once, cmap) == once


# ------------------------------------------------------- Fortran text

def test_fortran_call_sites_renamed_case_preserving():
    fmap = build_fortran_map(NAMES)
    src = textwrap.dedent("""\
              CALL DESCINIT( DESCA, M, N, MB, NB, 0, 0, ICTXT, LDA, INFO )
              np = numroc( n, nb, myrow, 0, nprow )
              CALL PXERBLA( ICTXT, 'PMGEMM', INFO )
    """)
    out = privatize_fortran_text(src, fmap)
    assert 'CALL EP_DESCINIT(' in out
    assert 'ep_numroc(' in out
    assert 'CALL EP_PXERBLA(' in out


def test_fortran_reshape_call_renamed_gridreshape_untouched():
    # CALL RESHAPE is provably the external REDIST routine (the intrinsic
    # is a function and cannot be CALLed); GRIDRESHAPE-family names must
    # not be clipped by the RESHAPE substring.
    fmap = build_fortran_map(NAMES)
    src = textwrap.dedent("""\
              EXTERNAL           PXERBLA, RESHAPE
              CALL RESHAPE( ICTXT, 0, ICTXTB, 1, IONE, NP )
              CTXT = EP_SL_GRIDRESHAPE( DESCA( 2 ), 0, 1, 1, NP, IONE )
    """)
    out = privatize_fortran_text(src, fmap)
    assert 'CALL EP_RESHAPE(' in out
    assert 'PXERBLA, EP_RESHAPE' in out
    assert 'EP_SL_GRIDRESHAPE(' in out
    assert 'EP_EP_' not in out


def test_fortran_function_definition_renamed():
    fmap = build_fortran_map(NAMES)
    src = textwrap.dedent("""\
              INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
              NUMROC = NEXTRA
              END
    """)
    out = privatize_fortran_text(src, fmap)
    assert 'INTEGER FUNCTION EP_NUMROC(' in out
    assert 'EP_NUMROC = NEXTRA' in out


# ------------------------------------------------------ tree dispatch

def test_privatize_tree_dispatches_by_extension(tmp_path):
    (tmp_path / 'grid.c').write_text(
        'void Cblacs_gridinit(Int *c, char *o, Int p, Int q);\n')
    (tmp_path / 'Bdef.h').write_text('extern BLACBUFF BI_AuxBuff;\n')
    (tmp_path / 'pmgetrf.f').write_text(
        '      CALL DESCINIT( DESCA, M, N, MB, NB, 0, 0, ICTXT, LDA, I )\n')
    (tmp_path / 'helper.F').write_text('      NP = NUMROC( N, NB, 0, 0, P )\n')
    (tmp_path / 'notes.txt').write_text('blacs_gridinit_ stays here\n')

    changed = privatize_tree(tmp_path, NAMES)
    assert changed == 4
    assert 'ep_Cblacs_gridinit' in (tmp_path / 'grid.c').read_text()
    assert 'ep_BI_AuxBuff' in (tmp_path / 'Bdef.h').read_text()
    assert 'EP_DESCINIT' in (tmp_path / 'pmgetrf.f').read_text()
    assert 'EP_NUMROC' in (tmp_path / 'helper.F').read_text()
    assert (tmp_path / 'notes.txt').read_text() == 'blacs_gridinit_ stays here\n'
    # second run is a no-op
    assert privatize_tree(tmp_path, NAMES) == 0


def test_privatize_tree_skips_pristine_split_originals(tmp_path):
    # The header-split step (task: Netlib-pristine public C headers)
    # restores split originals to upstream bytes and moves the
    # transformed content to a <pair>-prefixed sibling. The privatize
    # pass must leave the originals byte-identical and rename only the
    # sibling.
    pristine = 'void Cblacs_gridinit(Int *c, char *o, Int p, Int q);\n'
    (tmp_path / 'PBblacs.h').write_text(pristine)
    (tmp_path / 'mwPBblacs.h').write_text(pristine)
    (tmp_path / 'caller.c').write_text(
        '#include "mwPBblacs.h"\nCblacs_gridinit(&c, "R", p, q);\n')

    changed = privatize_tree(
        tmp_path, NAMES, skip={tmp_path / 'PBblacs.h'})
    assert changed == 2
    assert (tmp_path / 'PBblacs.h').read_text() == pristine
    assert 'ep_Cblacs_gridinit' in (tmp_path / 'mwPBblacs.h').read_text()
    caller = (tmp_path / 'caller.c').read_text()
    assert '#include "mwPBblacs.h"' in caller
    assert 'ep_Cblacs_gridinit(&c' in caller


def test_privatize_tree_two_files_share_one_map(tmp_path):
    # Regression guard for the id()-keyed rename-pattern cache class of
    # bug (see fortran/renames.py): consecutive files privatized with
    # the same logical map must all be rewritten, and a subsequent tree
    # with a different manifest must not reuse the previous pattern.
    (tmp_path / 'a.f').write_text('      CALL DESCINIT( D, M, N )\n')
    (tmp_path / 'b.f').write_text('      CALL DESCINIT( E, K, L )\n')
    assert privatize_tree(tmp_path, NAMES) == 2
    assert 'EP_DESCINIT' in (tmp_path / 'a.f').read_text()
    assert 'EP_DESCINIT' in (tmp_path / 'b.f').read_text()

    other = tmp_path / 'other'
    other.mkdir()
    (other / 'c.f').write_text('      CALL PXERBLA( ICTXT, SRNAME, INFO )\n')
    (other / 'd.f').write_text('      CALL DESCINIT( D, M, N )\n')
    assert privatize_tree(other, frozenset({'pxerbla_'})) == 1
    assert 'EP_PXERBLA' in (other / 'c.f').read_text()
    # descinit_ absent from the second manifest → untouched
    assert 'CALL DESCINIT' in (other / 'd.f').read_text()


# ------------------------------------------------------ recipe loader

def test_recipe_loads_privatize_manifest(tmp_path):
    recipe_dir = tmp_path / 'recipes'
    recipe_dir.mkdir()
    (tmp_path / 'src').mkdir()
    (recipe_dir / 'privatize.txt').write_text(textwrap.dedent("""\
        # comment line
        blacs_gridinit_
        Cblacs_gridinit

        numroc_
    """))
    recipe = recipe_dir / 'r.yaml'
    recipe.write_text(textwrap.dedent("""\
        library: blacs
        language: c
        source_dir: src
        privatize_symbols: privatize.txt
    """))
    config = load_recipe(recipe, project_root=tmp_path)
    assert config.privatize_symbols == frozenset(
        {'blacs_gridinit_', 'Cblacs_gridinit', 'numroc_'})


def test_recipe_without_privatize_key_defaults_empty(tmp_path):
    recipe_dir = tmp_path / 'recipes'
    recipe_dir.mkdir()
    (tmp_path / 'src').mkdir()
    recipe = recipe_dir / 'r.yaml'
    recipe.write_text(textwrap.dedent("""\
        library: blas
        language: fortran
        source_dir: src
    """))
    config = load_recipe(recipe, project_root=tmp_path)
    assert config.privatize_symbols == frozenset()
