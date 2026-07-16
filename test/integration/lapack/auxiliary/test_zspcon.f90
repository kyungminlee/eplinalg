program test_zspcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_complex_symmetric_quad, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zsptrf, target_zspcon
    use ref_quad_lapack, only: zsptrf, zspcon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info, np
    integer, allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep), allocatable :: A0(:,:), AP_ref(:), AP_got(:), work(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('zspcon', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_complex_symmetric_quad(n, A0, seed = 127001 + 47 * i)
        anorm = maxval(sum(abs(A0), dim=1))
        allocate(AP_ref(np), AP_got(np), ipiv_ref(n), ipiv_got(n), work(2*n))
        call pack_herm_packed_quad('U', n, A0, AP_ref); AP_got = AP_ref
        call zsptrf('U', n, AP_ref, ipiv_ref, info)
        call target_zsptrf('U', n, AP_got, ipiv_got, info)
        call zspcon('U', n, AP_ref, ipiv_ref, anorm, rcond_ref, work, info)
        call target_zspcon('U', n, AP_got, ipiv_got, anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(AP_ref, AP_got, ipiv_ref, ipiv_got, work)
    end do
    call report_finalize()
end program test_zspcon
