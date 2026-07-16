program test_dpocon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_spd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dpotrf, target_dpocon
    use ref_quad_lapack, only: dpotrf, dpocon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info
    integer, allocatable :: iwork(:)
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('dpocon', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_spd_matrix_quad(n, A0, seed = 112001 + 47 * i)
        anorm = maxval(sum(abs(A0), dim=1))
        allocate(A_ref(n,n), A_got(n,n), work(3*n), iwork(n))
        A_ref = A0; A_got = A0
        call dpotrf('U', n, A_ref, n, info)
        call target_dpotrf('U', n, A_got, n, info)
        call dpocon('U', n, A_ref, n, anorm, rcond_ref, work, iwork, info)
        call target_dpocon('U', n, A_got, n, anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 100.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, work, iwork)
    end do
    call report_finalize()
end program test_dpocon
