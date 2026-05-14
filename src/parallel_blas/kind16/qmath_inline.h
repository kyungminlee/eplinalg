/*
 * qmath_inline.h — always-inline soft-float __float128 arithmetic.
 *
 * Wraps the soft-fp template macros (vendored under
 * external/libgcc-softfp/) into static inline functions, so the
 * arithmetic body inlines into the GEMM hot loop instead of going
 * through libgcc's `__multf3` / `__addtf3` external entry points.
 *
 * Same numerical behavior as the libgcc entry points (same source
 * macros, same exception/rounding model). Tradeoff: significant
 * code growth in the consumer (each qmul / qadd expands to ~700 /
 * ~1200 x86_64 insns). Worth it iff the inlining lets gcc share
 * unpack/normalize work across adjacent ops, and/or the
 * PLT-trampoline cost was a meaningful share of the per-op time.
 *
 * The shim itself is GPL-incompatible-only-in-the-style-of-libgcc:
 * the underlying macros are LGPL-2.1+ with libgcc's linking
 * exception, so including this header in any program is fine.
 */
#ifndef QMATH_INLINE_H
#define QMATH_INLINE_H

#include "soft-fp.h"
#include "quad.h"

static inline __attribute__((always_inline)) __float128
qmul(__float128 a, __float128 b)
{
    FP_DECL_EX;
    FP_DECL_Q(A);
    FP_DECL_Q(B);
    FP_DECL_Q(R);
    __float128 r;

    FP_INIT_ROUNDMODE;
    FP_UNPACK_Q(A, a);
    FP_UNPACK_Q(B, b);
    FP_MUL_Q(R, A, B);
    FP_PACK_Q(r, R);
    FP_HANDLE_EXCEPTIONS;
    return r;
}

static inline __attribute__((always_inline)) __float128
qadd(__float128 a, __float128 b)
{
    FP_DECL_EX;
    FP_DECL_Q(A);
    FP_DECL_Q(B);
    FP_DECL_Q(R);
    __float128 r;

    FP_INIT_ROUNDMODE;
    FP_UNPACK_SEMIRAW_Q(A, a);
    FP_UNPACK_SEMIRAW_Q(B, b);
    FP_ADD_Q(R, A, B);
    FP_PACK_SEMIRAW_Q(r, R);
    FP_HANDLE_EXCEPTIONS;
    return r;
}

static inline __attribute__((always_inline)) __float128
qsub(__float128 a, __float128 b)
{
    FP_DECL_EX;
    FP_DECL_Q(A);
    FP_DECL_Q(B);
    FP_DECL_Q(R);
    __float128 r;

    FP_INIT_ROUNDMODE;
    FP_UNPACK_SEMIRAW_Q(A, a);
    FP_UNPACK_SEMIRAW_Q(B, b);
    FP_SUB_Q(R, A, B);
    FP_PACK_SEMIRAW_Q(r, R);
    FP_HANDLE_EXCEPTIONS;
    return r;
}

/* FMA via soft-fp's _FP_FMA template would land here, but that
 * macro needs scratch storage the wrapper-style entry point pattern
 * doesn't provide cleanly. Keep the two-rounding form
 * `qadd(qmul(a,b), c)` as the inner-kernel fma surrogate. */

#endif /* QMATH_INLINE_H */
