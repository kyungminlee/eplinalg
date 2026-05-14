# Vendored soft-fp headers from GCC's libgcc

These are unmodified copies of the soft-fp template headers from
GCC's `libgcc/soft-fp/` directory (and the x86_64-specific
`libgcc/config/i386/{,64/}sfp-machine.h`). They implement IEEE 754
binary128 (`__float128` / TFmode) arithmetic from a set of macros
parameterized by precision width.

In libgcc, these headers are combined with thin entry-point `.c`
files (`multf3.c`, `addtf3.c`, …) to produce the exported symbols
`__multf3`, `__addtf3`, etc. that gcc emits calls to for every
`__float128` operator. We use the same headers a different way: we
include them once into our overlay and wrap the macros in
`static inline __attribute__((always_inline))` functions, so the
arithmetic body inlines into the GEMM hot loop instead of going
through libgcc's external entry point.

## License

LGPL-2.1+ with the **libgcc linking exception** (see the copyright
header at the top of every file). The linking exception explicitly
permits including these files in any program — closed-source, MIT,
whatever — without invoking the LGPL's source-disclosure
requirement on the rest of the codebase. Modifications to these
files would still be LGPL'd; we keep them unmodified to preserve
that boundary.

## Files

| file | purpose |
|---|---|
| `soft-fp.h`    | top-level entry point; defines `_FP_DECL_EX`, the round-mode machinery, exception flag handling |
| `op-1.h`       | 1-limb fraction ops (used for single-precision) |
| `op-2.h`       | 2-limb fraction ops (used for `__float128` on 64-bit hosts) |
| `op-common.h`  | precision-agnostic `_FP_MUL`, `_FP_ADD`, `_FP_DIV` macros |
| `quad.h`       | TFmode parameters: exponent bits, fraction bits, etc. |
| `sfp-machine.h` + `config/i386/64/sfp-machine.h` | x86_64 machine-specific hooks (64x64→128 multiply via `mulq` inline asm) |

Source: GCC 12.x (matches Ubuntu 24.04 toolchain ABI). If we move to
a different toolchain the headers can be re-vendored from
`gcc/libgcc/soft-fp/` and `gcc/libgcc/config/i386/`.

## Inline shim

`src/parallel_blas/kind16/qmath_inline.h` defines the user-facing
`qmul`, `qadd`, `qsub`, `qfma` wrappers — these are the
inline-everywhere replacements for `__multf3` / `__addtf3` /
`__subtf3` calls.
