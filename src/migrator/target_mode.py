from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml


@dataclass(frozen=True)
class TargetMode:
    """Specifies how the migrator should transform Fortran source code."""

    # Identity
    name: str

    # Type declaration output
    real_type: str
    complex_type: str

    # Literal replacement strategy
    literal_mode: str                  # 'kind_suffix' or 'constructor'
    kind_suffix: Optional[int]         # 10 or 16 for KIND, None for module-based
    real_constructor: Optional[str]    # None for KIND, 'real64x2' for MF
    complex_constructor: Optional[str] # None for KIND, 'cmplx64x2' for MF

    # Intrinsic call transformation
    intrinsic_mode: str                # 'add_kind' or 'wrap_constructor'

    # Prefix map for renaming
    prefix_map: dict[str, str]

    # Module / constants support
    module_name: Optional[str]         # None for KIND, 'multifloats' for MF
    known_constants: dict[str, str]    # {'ZERO': 'MF_ZERO', ...}
    la_constants_map: dict[str, str]   # la_constants rename map
    la_constants_suffix: str           # '_EY' (kind10), '_QX' (kind16), '_MW' (multifloats)

    # Module public names (for building USE...ONLY clauses)
    module_type_names: frozenset[str] = field(default_factory=frozenset)
    module_constant_names: frozenset[str] = field(default_factory=frozenset)
    module_generic_names: frozenset[str] = field(default_factory=frozenset)
    module_operator_generics: tuple[str, ...] = ()

    # C-interop fields.  Used by c_migrator to substitute types / MPI
    # handles / reduction ops in cloned BLACS and PBLAS C sources.
    # All targets (KIND and module-based) populate these via YAML.
    c_real_type: Optional[str] = None         # 'QREAL' / 'float64x2'
    c_complex_type: Optional[str] = None      # 'XCOMPLEX' / 'complex64x2'
    c_c_real_type: Optional[str] = None       # '__float128' / 'float64x2'
    c_mpi_real: Optional[str] = None          # 'MPI_REAL16' / 'MPI_FLOAT64X2'
    c_mpi_complex: Optional[str] = None       # 'MPI_COMPLEX32' / 'MPI_COMPLEX64X2'
    c_mpi_sum_real: Optional[str] = None      # 'MPI_SUM' / 'MPI_DD_SUM'
    c_mpi_sum_complex: Optional[str] = None   # 'MPI_SUM' / 'MPI_ZZ_SUM'
    c_mpi_max_real: Optional[str] = None      # 'MPI_MAX' / 'MPI_DD_AMX'
    c_mpi_max_complex: Optional[str] = None   # 'MPI_MAX' / 'MPI_ZZ_AMX'
    c_mpi_min_real: Optional[str] = None      # 'MPI_MIN' / 'MPI_DD_AMN'
    c_mpi_min_complex: Optional[str] = None   # 'MPI_MIN' / 'MPI_ZZ_AMN'
    c_mpi_module: Optional[str] = None        # 'multifloats_mpi_f' / 'quad_mpi_f'
    c_needs_mpi_check: bool = False           # True only for KIND=16
    c_header_mode: Optional[str] = None       # 'typedef' or 'include'
    c_header: Optional[str] = None            # 'multifloats_bridge.h'

    @property
    def is_kind_based(self) -> bool:
        return self.kind_suffix is not None

    @property
    def module_public_names(self) -> frozenset[str]:
        """All public names exported by the module (types + constants + generics)."""
        return self.module_type_names | self.module_constant_names | self.module_generic_names


def load_target(name_or_path: str, project_root: Path | None = None) -> TargetMode:
    """Load a TargetMode from a YAML file or built-in name.

    Resolution order:
      1. If name_or_path is a path to an existing .yaml file, load it directly.
      2. Look for targets/{name_or_path}.yaml relative to project_root.
      3. If name_or_path is a bare integer (e.g. "16"), try targets/kind{n}.yaml.
    """
    if project_root is None:
        # Default: repo root is two levels up from this file
        # (src/migrator/target_mode.py → repo root)
        project_root = Path(__file__).resolve().parent.parent.parent

    # 1. Direct path
    p = Path(name_or_path)
    if p.suffix in ('.yaml', '.yml') and p.exists():
        return _load_target_yaml(p)

    # 2. Named lookup
    targets_dir = project_root / 'targets'
    candidate = targets_dir / f'{name_or_path}.yaml'
    if candidate.exists():
        return _load_target_yaml(candidate)

    # 3. Bare integer → kindN
    if name_or_path.isdigit():
        candidate = targets_dir / f'kind{name_or_path}.yaml'
        if candidate.exists():
            return _load_target_yaml(candidate)

    raise FileNotFoundError(
        f"No target YAML found for '{name_or_path}'. "
        f"Searched: {targets_dir / f'{name_or_path}.yaml'}"
    )


def _load_target_yaml(path: Path) -> TargetMode:
    """Parse a target YAML file and construct a TargetMode."""
    with open(path) as f:
        d = yaml.safe_load(f)

    # Module section: can be null (KIND targets) or a dict (module-based)
    mod = d.get('module')
    if isinstance(mod, dict):
        module_name = mod['name']
        module_type_names = frozenset(mod.get('type_names') or [])
        module_constant_names = frozenset(mod.get('constant_names') or [])
        module_generic_names = frozenset(mod.get('generic_names') or [])
        module_operator_generics = tuple(mod.get('operator_generics') or [])
    else:
        module_name = mod  # None or a bare string (backward compat)
        module_type_names = frozenset()
        module_constant_names = frozenset()
        module_generic_names = frozenset()
        module_operator_generics = ()

    c = d.get('c_interop') or {}

    return TargetMode(
        name=d['name'],
        real_type=d['types']['real'],
        complex_type=d['types']['complex'],
        literal_mode=d['literals']['mode'],
        kind_suffix=d['literals'].get('kind_suffix'),
        real_constructor=d['literals'].get('real_constructor'),
        complex_constructor=d['literals'].get('complex_constructor'),
        intrinsic_mode=d['intrinsics']['mode'],
        prefix_map=d['prefixes'],
        module_name=module_name,
        known_constants=d.get('known_constants') or {},
        la_constants_map=d.get('la_constants_map') or {},
        la_constants_suffix=d.get('la_constants_suffix', ''),
        module_type_names=module_type_names,
        module_constant_names=module_constant_names,
        module_generic_names=module_generic_names,
        module_operator_generics=module_operator_generics,
        c_real_type=c.get('real_type'),
        c_complex_type=c.get('complex_type'),
        c_c_real_type=c.get('c_real_type'),
        c_mpi_real=c.get('mpi_real'),
        c_mpi_complex=c.get('mpi_complex'),
        c_mpi_sum_real=c.get('mpi_sum_real'),
        c_mpi_sum_complex=c.get('mpi_sum_complex'),
        c_mpi_max_real=c.get('mpi_max_real'),
        c_mpi_max_complex=c.get('mpi_max_complex'),
        c_mpi_min_real=c.get('mpi_min_real'),
        c_mpi_min_complex=c.get('mpi_min_complex'),
        c_mpi_module=c.get('mpi_module'),
        c_needs_mpi_check=c.get('needs_mpi_check', False),
        c_header_mode=c.get('header_mode'),
        c_header=c.get('header'),
    )
