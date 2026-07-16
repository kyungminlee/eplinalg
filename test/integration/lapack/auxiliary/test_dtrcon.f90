program test_dtrcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtrcon
    use ref_quad_lapack, only: dtrcon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info, j
    integer, allocatable :: iwork(:)
    real(ep), allocatable :: A(:,:), work(:)
    real(ep) :: rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('dtrcon', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 128001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + real(2*n, ep); end do
        allocate(work(3*n), iwork(n))
        call dtrcon('1', 'U', 'N', n, A, n, rcond_ref, work, iwork, info)
        call target_dtrcon('1', 'U', 'N', n, A, n, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, work, iwork)
    end do
    call report_finalize()
end program test_dtrcon
