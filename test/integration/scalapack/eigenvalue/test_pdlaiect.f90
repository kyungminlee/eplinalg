program test_pdlaiect
    ! PDLAIECT — direct unit test for the C-side IEEE sign-bit
    ! helpers (p?lasnbt, p?laiectb, p?laiectl, p?lachkieee). These
    ! are the Sturm-sequence count routines that pdstebz / pdsyevx
    ! fast-path on; upstream pdlaiect.c does sign extraction by
    ! pointer-casting double* to unsigned int* and shifting bit 31.
    !
    ! The migrator's auto-clones for kind10 (EREAL = long double),
    ! kind16 (QREAL = __float128), and multifloats (float64x2)
    ! inherit that bit-twiddle, which is meaningless on those types.
    ! The recipe overrides (codegen/recipes/scalapack/{e,q,m}fc_overrides/
    ! p{e,q,m}laiect.c) replace the body with native ``<`` comparison.
    !
    ! Two checks per routine:
    !   1. pdlasnbt   returns IEFLAG = 1  (was 0 from the auto-clone;
    !                                       proves override is wired in)
    !   2. pdlachkieee returns ISIEEE = 1 (was 0 / UB from the auto-clone)
    !   3. pdlaiectb on a 3x3 identity tridiagonal returns the
    !      mathematically correct count for two shift values
    !      (0.5 → 0; 1.5 → 3)
    !   4. pdlaiectl returns the same answer (both routines
    !      forward to the same body in the override)
    use prec_kinds,        only: ep
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank
    use target_scalapack,  only: target_name, target_pdlasnbt,         &
                                 target_pdlaiectb, target_pdlaiectl,   &
                                 target_pdlachkieee
    implicit none

    integer  :: ieflag, isieee, count_b, count_l
    real(ep) :: d(5), d2(3), d1(1)
    real(ep) :: rmax, rmin
    real(ep) :: err

    call grid_init()
    call report_init('pdlaiect', target_name, my_rank)

    if (my_rank == 0) then
        ! 1. pdlasnbt must report fast-path is available.
        call target_pdlasnbt(ieflag)
        err = real(abs(ieflag - 1), ep)
        call report_case('pdlasnbt:ieflag=1', err, 0.0_ep)

        ! 2. pdlachkieee must report IEEE compliance.
        !    rmax / rmin are unused by the override; pass plausible
        !    values to stay clear of FPU traps in case any future
        !    override re-introduces actual checks.
        rmax = huge(0.0_ep) * 0.5_ep
        rmin = tiny(0.0_ep)
        call target_pdlachkieee(isieee, rmax, rmin)
        err = real(abs(isieee - 1), ep)
        call report_case('pdlachkieee:isieee=1', err, 0.0_ep)

        ! 3. pdlaiectb on tridiag I_3.
        !    Interleaved layout: d(1)=D1, d(2)=E1^2, d(3)=D2,
        !                       d(4)=E2^2, d(5)=D3. Length 2N-1 = 5.
        d = [ 1.0_ep, 0.0_ep, 1.0_ep, 0.0_ep, 1.0_ep ]

        !    Eigenvalues of T = I_3 are {1, 1, 1}.
        !      sigma=0.5 → 0 eigenvalues ≤ sigma → count = 0
        !      sigma=1.5 → 3 eigenvalues ≤ sigma → count = 3
        call target_pdlaiectb(0.5_ep, 3, d, count_b)
        err = real(abs(count_b - 0), ep)
        call report_case('pdlaiectb:sigma=0.5:count=0', err, 0.0_ep)

        call target_pdlaiectb(1.5_ep, 3, d, count_b)
        err = real(abs(count_b - 3), ep)
        call report_case('pdlaiectb:sigma=1.5:count=3', err, 0.0_ep)

        ! 4. pdlaiectl agrees with pdlaiectb (overrides forward both
        !    big-endian and little-endian variants to the same body).
        call target_pdlaiectl(0.5_ep, 3, d, count_l)
        err = real(abs(count_l - 0), ep)
        call report_case('pdlaiectl:sigma=0.5:count=0', err, 0.0_ep)

        call target_pdlaiectl(1.5_ep, 3, d, count_l)
        err = real(abs(count_l - 3), ep)
        call report_case('pdlaiectl:sigma=1.5:count=3', err, 0.0_ep)

        ! 5. Distinct eigenvalues with non-zero off-diagonal — exercises
        !    the Sturm recurrence ``tmp - pe2/prev_tmp`` properly. I_3
        !    (case 3) has E_i^2 = 0, so the recurrence's off-diagonal term
        !    vanishes; a regression to the buggy bit-twiddle could return
        !    the same answer for I_3 by coincidence. Distinct eigenvalues
        !    force the loop body to produce diverse tmp values and detect
        !    that case.
        !
        !    T = tridiag([2,2], diag=[2,2,2], [2,2]) sized n=2 with D=2,
        !    off-diagonal=1 — eigenvalues are {1, 3} exactly.
        d2 = [ 2.0_ep, 1.0_ep, 2.0_ep ]
        call target_pdlaiectb(0.5_ep, 2, d2, count_b)
        err = real(abs(count_b - 0), ep)
        call report_case('pdlaiectb:distinct:sigma=0.5:count=0', err, 0.0_ep)

        call target_pdlaiectb(2.5_ep, 2, d2, count_b)
        err = real(abs(count_b - 1), ep)
        call report_case('pdlaiectb:distinct:sigma=2.5:count=1', err, 0.0_ep)

        call target_pdlaiectb(3.5_ep, 2, d2, count_b)
        err = real(abs(count_b - 2), ep)
        call report_case('pdlaiectb:distinct:sigma=3.5:count=2', err, 0.0_ep)

        ! 6. n=1 edge case — Sturm loop body never executes; only the
        !    pre-loop initialization runs. Guards against refactor bugs
        !    in loop-bound handling.
        d1 = [ 5.0_ep ]
        call target_pdlaiectb(0.0_ep, 1, d1, count_b)
        err = real(abs(count_b - 0), ep)
        call report_case('pdlaiectb:n=1:sigma<lambda:count=0', err, 0.0_ep)

        call target_pdlaiectb(10.0_ep, 1, d1, count_b)
        err = real(abs(count_b - 1), ep)
        call report_case('pdlaiectb:n=1:sigma>lambda:count=1', err, 0.0_ep)
    end if

    call report_finalize()
    call grid_exit()
end program test_pdlaiect
