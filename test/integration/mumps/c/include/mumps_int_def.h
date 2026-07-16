/*
 * Stand-in for upstream's generated `mumps_int_def.h`.
 *
 * Upstream's build system synthesizes this file from
 * `extern/MUMPS_5.8.2/src/mumps_int_def{32,64}_h.in` based on whether
 * `-DINTSIZE64` was passed. The migrator and the fm-mumps cmake glue
 * do NOT pass `-DINTSIZE64`, so 32-bit MUMPS_INT is the only path the
 * bridge needs to support — pin it here so the override
 * `mumps_c_types.h` resolves `#include "mumps_int_def.h"` cleanly.
 */

#ifndef MUMPS_INT_H
#define MUMPS_INT_H
#define MUMPS_INTSIZE32
#endif
