/* Task 44: the extended targets (kind10 / kind16 / multifloats) link the
 * ep_-privatized BLACS/ScaLAPACK-common engine, whose family-independent
 * entry points carry an ep_ prefix so they can never collide with an MKL
 * linked into the same consumer. Hand-written test sources that call
 * those routines by their pristine names include this header (and are
 * named .F90 so the preprocessor runs); the build defines
 * EPLINALG_TEST_EP_BLACS when the staged stack is privatized (see
 * tests/CMakeLists.txt). Baseline (kind4/kind8) stagings leave the macro
 * undefined and keep the pristine Netlib names.
 *
 * All call sites are lowercase; cpp is case-sensitive, so keep any new
 * uses lowercase too.
 */
#ifdef EPLINALG_TEST_EP_BLACS
#define blacs_pinfo    ep_blacs_pinfo
#define blacs_get      ep_blacs_get
#define blacs_set      ep_blacs_set
#define blacs_barrier  ep_blacs_barrier
#define blacs_pcoord   ep_blacs_pcoord
#define blacs_pnum     ep_blacs_pnum
#define blacs_gridinit ep_blacs_gridinit
#define blacs_gridinfo ep_blacs_gridinfo
#define blacs_gridexit ep_blacs_gridexit
#define blacs_exit     ep_blacs_exit
#define igesd2d        ep_igesd2d
#define igerv2d        ep_igerv2d
#define igebs2d        ep_igebs2d
#define igebr2d        ep_igebr2d
#define igsum2d        ep_igsum2d
#define igamx2d        ep_igamx2d
#define igamn2d        ep_igamn2d
#endif
