program test_dptcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_dpttrf, target_dptcon
    use ref_quad_lapack, only: dpttrf, dptcon
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer :: i, n, info, j
    real(ep), allocatable :: d0(:), e0(:), d_ref(:), e_ref(:), d_got(:), e_got(:), work(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('dptcon', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n,   d0, seed = 124001 + 47 * i)
        call gen_vector_quad(n-1, e0, seed = 124011 + 47 * i)
        do j = 1, n; d0(j) = abs(d0(j)) + real(4, ep); end do
        anorm = 0.0_ep
        do j = 1, n
            anorm = max(anorm, &
                merge(abs(e0(j-1)), 0.0_ep, j > 1) + abs(d0(j)) + &
                merge(abs(e0(j)),   0.0_ep, j < n))
        end do
        allocate(d_ref(n), e_ref(n-1), d_got(n), e_got(n-1), work(n))
        d_ref = d0; e_ref = e0; d_got = d0; e_got = e0
        call dpttrf(n, d_ref, e_ref, info)
        call target_dpttrf(n, d_got, e_got, info)
        call dptcon(n, d_ref, e_ref, anorm, rcond_ref, work, info)
        call target_dptcon(n, d_got, e_got, anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(d0, e0, d_ref, e_ref, d_got, e_got, work)
    end do
    call report_finalize()
end program test_dptcon
