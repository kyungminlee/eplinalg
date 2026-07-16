program test_drotmg
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: rel_err_scalar
    use target_blas,   only: target_name, target_eps, target_drotmg
    use ref_quad_blas, only: drotmg
    implicit none

    ! Each row is (d1, d2, x1, y1) chosen to drive a specific flag /
    ! rescale branch in DROTMG:
    !   1–4: ordinary inputs spanning flag=0 / flag=1 territory
    !   5  : d1 < 0  → forces flag = -1 (full rescale)
    !   6  : d2 = 0  → degenerate axis, exercises early return
    !   7  : y1 = 0  → flag = -2 (identity rotation)
    !   8  : extreme |q1| → triggers GAMSQ rescale path inside drotmg
    !   9  : extreme 1/|q1| → triggers 1/GAMSQ rescale path
    integer, parameter :: ncases = 9
    real(ep), parameter :: gam = 4096.0_ep
    real(ep), parameter :: cases(4, ncases) = reshape([ &
        1.5_ep, 3.2_ep, 0.5_ep, 0.3_ep, &
        3.0_ep, 2.4_ep, 1.0_ep, 0.6_ep, &
        4.5_ep, 1.6_ep, 1.5_ep, 0.9_ep, &
        6.0_ep, 0.8_ep, 2.0_ep, 1.2_ep, &
        -1.0_ep, 1.0_ep, 0.5_ep, 0.3_ep, &
        2.0_ep, 0.0_ep, 0.5_ep, 0.3_ep, &
        2.0_ep, 1.0_ep, 0.5_ep, 0.0_ep, &
        1.0_ep, 1.0_ep, 1.0_ep / gam, gam, &
        1.0_ep, 1.0_ep, gam, 1.0_ep / gam], &
        [4, ncases])
    integer :: i
    real(ep) :: d1_0, d2_0, x1_0, y1_0
    real(ep) :: d1_ref, d2_ref, x1_ref, param_ref(5)
    real(ep) :: d1_got, d2_got, x1_got, param_got(5)
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('drotmg', target_name)
    do i = 1, ncases
        d1_0 = cases(1, i); d2_0 = cases(2, i)
        x1_0 = cases(3, i); y1_0 = cases(4, i)
        d1_ref = d1_0; d2_ref = d2_0; x1_ref = x1_0
        d1_got = d1_0; d2_got = d2_0; x1_got = x1_0
        call drotmg(d1_ref, d2_ref, x1_ref, y1_0, param_ref)
        call target_drotmg(d1_got, d2_got, x1_got, y1_0, param_got)
        ! Compare d1, d2, x1, and the FLAG (param(1)). Other param
        ! slots are flag-conditional (BLAS spec: p(2)/p(5) only used
        ! when flag=1; p(3)/p(4) only when flag=0; etc.) — comparing
        ! the unused slots is a false positive.
        err = max(rel_err_scalar(d1_got, d1_ref), &
                  rel_err_scalar(d2_got, d2_ref), &
                  rel_err_scalar(x1_got, x1_ref), &
                  rel_err_scalar(param_got(1), param_ref(1)))
        call compare_meaningful_params(param_ref, param_got, err)
        tol = 32.0_ep * target_eps
        write(label, '(a,i0)') 'case=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()

contains

    ! Compare only the param slots that are meaningful for the
    ! returned flag value. Updates `err` with the max relative error
    ! across the meaningful slots.
    subroutine compare_meaningful_params(p_ref, p_got, err)
        real(ep), intent(in)    :: p_ref(5), p_got(5)
        real(ep), intent(inout) :: err
        integer :: flag

        flag = nint(p_ref(1))
        select case (flag)
        case (-1)
            ! All four off-diagonal/diagonal entries are meaningful
            err = max(err, rel_err_scalar(p_got(2), p_ref(2)), &
                          rel_err_scalar(p_got(3), p_ref(3)), &
                          rel_err_scalar(p_got(4), p_ref(4)), &
                          rel_err_scalar(p_got(5), p_ref(5)))
        case (0)
            ! H = [1, h12; h21, 1] — only p(3)=h12, p(4)=h21 used
            err = max(err, rel_err_scalar(p_got(3), p_ref(3)), &
                          rel_err_scalar(p_got(4), p_ref(4)))
        case (1)
            ! H = [h11, 1; -1, h22] — only p(2)=h11, p(5)=h22 used
            err = max(err, rel_err_scalar(p_got(2), p_ref(2)), &
                          rel_err_scalar(p_got(5), p_ref(5)))
        case (-2)
            ! H = identity — no off-diagonal params used
        end select
    end subroutine compare_meaningful_params

end program test_drotmg
