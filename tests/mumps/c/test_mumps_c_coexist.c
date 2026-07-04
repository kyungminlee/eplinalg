/* Coexistence link check for the two C bridges this release ships.
 *
 * A build emits FOUR ?mumps_c entry points:
 *   - the migrated real bridge   TARGET_REAL_MUMPS_C    (emumps_c / qmumps_c /
 *     mmumps_c, depending on the target width — always the upstream-`d` lineage)
 *   - the migrated complex bridge TARGET_COMPLEX_MUMPS_C (ymumps_c / xmumps_c /
 *     wmumps_c — the upstream-`z` lineage)
 *   - the genuine double real bridge    dmumps_c
 *   - the genuine double complex bridge  zmumps_c
 * The migrated pair are bundled into ${LIB_PREFIX}mumps_c; the genuine pair are
 * separate libdmumps_c / libzmumps_c. Both share the arithmetic-independent C
 * runtime factored into libmumps_c_runtime (mumps_common.c, mumps_addr.c, IO,
 * thread, numa, ... — none of which carry an arithmetic or width macro).
 *
 * This single translation unit references all four so the linker must resolve
 * (and, under --whole-archive, pull) every bridge object. The check is purely a
 * LINK check: with the runtime shared, no two archives may define the same
 * runtime symbol. The build wires this source TWICE — once per coherent family:
 *
 *   genuine  : libdmumps_c + libzmumps_c + libmumps_c_runtime
 *   migrated : ${LIB_PREFIX}mumps_c      + libmumps_c_runtime
 *
 * Each family must link with NO "multiple definition" diagnostic. Before the
 * runtime was factored out, both bridges bundled their own copy and the genuine
 * pair (the standard upstream "use dmumps_c and zmumps_c in one program" case)
 * could not co-link. The two families are deliberately NOT linked together: the
 * migrated real bridge is compiled at MUMPS_ARITH_d and so emits the very same
 * dmumps_* Fortran-callback helpers as the genuine d bridge (and likewise z),
 * and their Fortran commons collide — exactly upstream's "one libmumps_common
 * per program" constraint. Cross-family exclusion is expected, not a defect.
 *
 * Forward-declared with void* (not the real *MUMPS_STRUC_C) on purpose: only the
 * symbol references matter here, no struct layout is needed, and this keeps the
 * probe independent of which width the migrated headers describe.
 */

void TARGET_REAL_MUMPS_C(void *);
void TARGET_COMPLEX_MUMPS_C(void *);
void dmumps_c(void *);
void zmumps_c(void *);

/* Address-taken so nothing is optimized away and every bridge object is
 * genuinely referenced. Never invoked — this is a link-only probe. */
static void (*const all_bridges[])(void *) = {
    TARGET_REAL_MUMPS_C,
    TARGET_COMPLEX_MUMPS_C,
    dmumps_c,
    zmumps_c,
};

int main(void)
{
    /* Return 0 (success) whenever the entry points have real addresses, which
     * they always do once linked. The point is that we got here at all — i.e.
     * the link succeeded without a multiple-definition error. */
    unsigned long acc = 0;
    for (unsigned i = 0; i < sizeof(all_bridges) / sizeof(all_bridges[0]); ++i)
        acc |= (unsigned long)(void *)all_bridges[i];
    return acc == 0;
}
