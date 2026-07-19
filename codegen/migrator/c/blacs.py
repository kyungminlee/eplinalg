"""BLACS-specific C migration: hardcoded d/z directory cloning,
Bdef.h ABI patching, and the MPI_REAL16 availability check module.
"""

import shutil
from pathlib import Path

from ..target_mode import TargetMode
from ..templates import build_sub_vars
from .clone import (
    classify_blacs_stem,
    clone_c_file,
    derive_routine_renames,
    rename_c_file,
)


def migrate_blacs_c_directory(src_dir: Path, output_dir: Path,
                              target_mode: TargetMode,
                              skip_files: set[str] | None = None) -> dict:
    """BLACS-specific C migration (original behavior).

    Reads from ``src_dir`` (the staged source tree, already patched).

    ``skip_files`` (uppercased logical stems, no trailing Fortran
    underscore) are dropped from the migrated output — used to remove the
    genuine S/C entry points (``sgesd2d_``, ``BI_cvvsum``, …) that the
    plain ``blacs`` standard archive and MKL already provide, so the
    migrated archive stays strictly extended-precision. Skipped stems
    are also never used as clone sources (they are S/C, not the D/Z the
    clone loops read), so dropping them cannot disturb the e/y clones.
    """
    _skip = skip_files or set()
    template_vars = build_sub_vars(target_mode)
    rp = template_vars['RP']
    cp = template_vars['CP']

    output_dir.mkdir(parents=True, exist_ok=True)

    sources = sorted(src_dir.iterdir())

    # Copy headers always; drop precision-prefixed .c originals (the
    # std archive owns the S/D/C/Z entry points); copy precision-
    # independent dispatcher .c files (BI_BlacsErr, BI_GetBuff,
    # blacs_setup_, ...). The
    # dispatchers ride in the migrated archive because callers reach
    # them through the BLACS-style symbol map; the std archive picks
    # them up too via the upstream source dir, but the linker only
    # ever pulls one definition (whichever resolves an undefined
    # symbol first), so duplication is benign and the migrated copies
    # match the precision the migrated bodies expect.
    def _is_blacs_precision_prefix(stem: str) -> bool:
        # Mirrors the d/z clone discovery in the loops below (the s/c
        # variants classify_blacs_stem also knows are copied through —
        # only d/z are clone sources here).
        plan = classify_blacs_stem(stem)
        return plan is not None and plan[0] in ('d', 'z')

    for f in sources:
        ext = f.suffix.lower()
        # ``skip_files`` drops genuine S/C sources from the migrated
        # output entirely (headers are never precision-specific, so they
        # are exempt). Match on the underscore-stripped, uppercased stem
        # to align with the recipe's logical-name convention.
        stem_upper = (f.stem[:-1] if f.stem.endswith('_') else f.stem).upper()
        if ext == '.c' and stem_upper in _skip:
            continue
        if ext == '.h':
            shutil.copy2(f, output_dir / f.name)
        elif ext == '.c':
            if not _is_blacs_precision_prefix(f.stem):
                shutil.copy2(f, output_dir / f.name)

    cloned = []

    # Clone d-variant → real-extended
    for f in sources:
        if f.suffix.lower() != '.c':
            continue
        plan = classify_blacs_stem(f.stem)
        if plan is None or plan[0] != 'd':
            continue
        _src_prefix, _is_complex, subs = plan
        new_name = rename_c_file(f.name, 'd', rp)
        if new_name == f.name:
            continue
        renames = derive_routine_renames(f.stem, Path(new_name).stem)
        clone_c_file(f, output_dir / new_name,
                     subs, template_vars, renames)
        cloned.append(f'{f.name} → {new_name}')

    # Clone z-variant → complex-extended
    for f in sources:
        if f.suffix.lower() != '.c':
            continue
        plan = classify_blacs_stem(f.stem)
        if plan is None or plan[0] != 'z':
            continue
        _src_prefix, _is_complex, subs = plan
        new_name = rename_c_file(f.name, 'z', cp)
        if new_name == f.name:
            continue
        renames = derive_routine_renames(f.stem, Path(new_name).stem)
        clone_c_file(f, output_dir / new_name,
                     subs, template_vars, renames)
        cloned.append(f'{f.name} → {new_name}')

    # Patch Bdef.h with extended-precision type definitions and macros
    bdef_path = output_dir / 'Bdef.h'
    if bdef_path.exists():
        patch_bdef_header(bdef_path, target_mode, template_vars)


    # Generate MPI datatype availability check when the target requires
    # stock MPI datatypes that may not be universally available (e.g.
    # MPI_REAL16 / MPI_COMPLEX32 for KIND=16).
    if target_mode.c_needs_mpi_check:
        generate_mpi_real16_check(output_dir)

    return {
        'cloned': cloned,
        'divergences': [],
        'template_vars': template_vars,
        'split_headers': {},
    }


# The 11 BLACS user-facing routine suffixes (same for each type prefix).
_BLACS_ROUTINE_SUFFIXES = [
    'gesd2d', 'gerv2d', 'gebs2d', 'gebr2d',
    'trsd2d', 'trrv2d', 'trbs2d', 'trbr2d',
    'gsum2d', 'gamx2d', 'gamn2d',
]


def patch_bdef_header(bdef_path: Path, target_mode: TargetMode,
                      template_vars: dict[str, str]) -> None:
    """Add extended-precision type definitions and macros to Bdef.h."""
    rp = template_vars['RP']
    cp = template_vars['CP']
    real_type = template_vars['REAL_TYPE']      # e.g. QREAL / float64x2_t
    complex_type = template_vars['COMPLEX_TYPE']  # e.g. XCOMPLEX / complex128x2_t

    c_type = template_vars['C_REAL_TYPE']

    if target_mode.c_header_mode == 'include':
        # Module-based target: types come from an external header
        # (e.g. multifloats_bridge.h) which provides full operator
        # overloading via C++ templates.
        type_block = f"""
/*
 *  Companion types for migrated {rp}/{cp} routines.
 *  Provided by {target_mode.c_header}.
 */
#include "{target_mode.c_header}"

void BI_{rp}mvcopy(Int m, Int n, {real_type} *A, Int lda, {real_type} *buff);
void BI_{rp}vmcopy(Int m, Int n, {real_type} *A, Int lda, {real_type} *buff);
"""
    else:
        # typedef mode (KIND targets): emit primitive typedefs.
        type_block = f"""
/*
 *  Extended-precision types for migrated {{prefix}} routines.
 */
typedef {c_type} {real_type};
typedef struct {{{real_type} r, i;}} {complex_type};

void BI_{rp}mvcopy(Int m, Int n, {real_type} *A, Int lda, {real_type} *buff);
void BI_{rp}vmcopy(Int m, Int n, {real_type} *A, Int lda, {real_type} *buff);
""".replace('{prefix}', f'{rp}/{cp}')

    # --- Complex copy macros ---
    macro_block = f"""#define BI_{cp}mvcopy(m, n, A, lda, buff) \\
        BI_{rp}mvcopy(2*(m), (n), ({real_type} *) (A), 2*(lda), ({real_type} *) (buff))
#define BI_{cp}vmcopy(m, n, A, lda, buff) \\
        BI_{rp}vmcopy(2*(m), (n), ({real_type} *) (A), 2*(lda), ({real_type} *) (buff))
"""

    # --- Fortran name mangling defines ---
    def _mangling_block(prefix: str, transform) -> str:
        lines = []
        for suf in _BLACS_ROUTINE_SUFFIXES:
            src = f'{prefix}{suf}_'
            dst = transform(f'{prefix}{suf}')
            lines.append(f'#define {src:19s}{dst}')
        return '\n'.join(lines) + '\n'

    nochange_block = (_mangling_block(rp, lambda s: s) +
                      _mangling_block(cp, lambda s: s))
    upcase_block = (_mangling_block(rp, str.upper) +
                    _mangling_block(cp, str.upper))

    text = bdef_path.read_text(errors='replace')

    # Insert type definitions after the DCOMPLEX/SCOMPLEX typedefs
    text = text.replace(
        'typedef struct {float r, i;} SCOMPLEX;\n',
        'typedef struct {float r, i;} SCOMPLEX;\n' + type_block
    )

    # Insert copy macros after the existing BI_zvmcopy macro
    zvmcopy_line = '#define BI_zvmcopy(m, n, A, lda, buff) \\\n        BI_dvmcopy(2*(m), (n), (double *) (A), 2*(lda), (double *) (buff))'
    text = text.replace(zvmcopy_line, zvmcopy_line + '\n' + macro_block)

    # Insert name mangling defines in NOCHANGE section (after zgamn2d)
    nochange_marker = '#define zgamn2d_   zgamn2d\n'
    text = text.replace(nochange_marker, nochange_marker + nochange_block)

    # Insert name mangling defines in UPCASE section (after ZGAMN2D)
    upcase_marker = '#define zgamn2d_   ZGAMN2D\n'
    text = text.replace(upcase_marker, upcase_marker + upcase_block)

    bdef_path.write_text(text)


def generate_mpi_real16_check(output_dir: Path) -> None:
    """Generate a CMake module that verifies MPI_REAL16 support."""
    cmake = """\
# CheckMpiReal16.cmake — Verify that the MPI implementation provides
# MPI_REAL16 and MPI_COMPLEX32, which are required by the migrated
# quad-precision BLACS routines.
#
# Usage:
#   include(CheckMpiReal16)
#   check_mpi_real16()          # FATAL_ERROR if unsupported
#
# After a successful check the cache variable HAVE_MPI_REAL16 is set.

include(CheckCSourceCompiles)

function(check_mpi_real16)
    find_package(MPI REQUIRED COMPONENTS C)

    set(CMAKE_REQUIRED_INCLUDES  ${MPI_C_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS     ${MPI_C_COMPILE_FLAGS})

    check_c_source_compiles("
        #include <mpi.h>
        int main(void) {
            MPI_Datatype dt_real = MPI_REAL16;
            MPI_Datatype dt_cplx = MPI_COMPLEX32;
            (void)dt_real; (void)dt_cplx;
            return 0;
        }
    " HAVE_MPI_REAL16)

    if(NOT HAVE_MPI_REAL16)
        message(FATAL_ERROR
            "The MPI implementation does not provide MPI_REAL16 / MPI_COMPLEX32. "
            "Quad-precision BLACS requires an MPI library built with 128-bit real "
            "support (e.g. OpenMPI or MPICH compiled with __float128 enabled).")
    endif()

    set(HAVE_MPI_REAL16 ${HAVE_MPI_REAL16} PARENT_SCOPE)
endfunction()
"""
    (output_dir / 'CheckMpiReal16.cmake').write_text(cmake)
