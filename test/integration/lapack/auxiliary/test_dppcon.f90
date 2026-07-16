program test_dppcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_spd_matrix_quad, pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, target_dpptrf, target_dppcon
    use ref_quad_lapack, only: dpptrf, dppcon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info, np
    integer, allocatable :: iwork(:)
    real(ep), allocatable :: A0(:,:), AP_ref(:), AP_got(:), work(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('dppcon', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_spd_matrix_quad(n, A0, seed = 122001 + 47 * i)
        anorm = maxval(sum(abs(A0), dim=1))
        allocate(AP_ref(np), AP_got(np), work(3*n), iwork(n))
        call pack_sym_packed_quad('U', n, A0, AP_ref); AP_got = AP_ref
        call dpptrf('U', n, AP_ref, info)
        call target_dpptrf('U', n, AP_got, info)
        call dppcon('U', n, AP_ref, anorm, rcond_ref, work, iwork, info)
        call target_dppcon('U', n, AP_got, anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(AP_ref, AP_got, work, iwork)
    end do
    call report_finalize()
end program test_dppcon
