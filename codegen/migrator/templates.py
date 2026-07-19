"""Template substitution engine shared across migration pipelines.

The ``{RP}/{CP}/{RPU}/{CPU}/{REAL_TYPE}/...`` placeholder vocabulary is
a contract spanning recipe YAML (``extra_renames``, ``c_type_aliases``,
``header_patches``), the C migrator's clone rules, and the Fortran
pipeline's extra-rename expansion — so the mechanism lives here rather
than as a private detail of any one consumer.
"""

from .prefix_classifier import target_prefix
from .target_mode import TargetMode


def build_sub_vars(target_mode: TargetMode) -> dict[str, str]:
    """Build template substitution variables for a given target mode.

    All values are read from the target's c_interop YAML section
    (populated in TargetMode by ``load_target()``).
    """
    assert target_mode.c_real_type is not None, (
        "target_mode missing c_real_type; ensure the target YAML has a "
        "c_interop section with real_type defined."
    )
    rp = target_prefix(target_mode, is_complex=False).lower()
    cp = target_prefix(target_mode, is_complex=True).lower()
    return {
        'REAL_TYPE': target_mode.c_real_type,
        'COMPLEX_TYPE': target_mode.c_complex_type,
        'C_REAL_TYPE': target_mode.c_c_real_type,
        'MPI_REAL': target_mode.c_mpi_real,
        'MPI_COMPLEX': target_mode.c_mpi_complex,
        'MPI_SUM_REAL': target_mode.c_mpi_sum_real,
        'MPI_SUM_COMPLEX': target_mode.c_mpi_sum_complex,
        'RP': rp,
        'CP': cp,
        'RPU': rp.upper(),
        'CPU': cp.upper(),
    }


def expand_template(s: str, template_vars: dict[str, str]) -> str:
    """Expand ``{KEY}`` placeholders in ``s`` using template_vars."""
    for key, val in template_vars.items():
        s = s.replace('{' + key + '}', val)
    return s
