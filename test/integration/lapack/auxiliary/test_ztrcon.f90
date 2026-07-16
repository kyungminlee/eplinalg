program test_ztrcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztrcon
    use ref_quad_lapack, only: ztrcon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info, j
    complex(ep), allocatable :: A(:,:), work(:)
    real(ep), allocatable :: rwork(:)
    real(ep) :: rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('ztrcon', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 129001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + cmplx(real(2*n, ep), 0.0_ep, ep); end do
        allocate(work(2*n), rwork(n))
        call ztrcon('1', 'U', 'N', n, A, n, rcond_ref, work, rwork, info)
        call target_ztrcon('1', 'U', 'N', n, A, n, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, work, rwork)
    end do
    call report_finalize()
end program test_ztrcon
