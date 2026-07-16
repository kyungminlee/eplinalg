program test_dtpcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_matrix_quad, pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, target_dtpcon
    use ref_quad_lapack, only: dtpcon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info, np, j
    integer, allocatable :: iwork(:)
    real(ep), allocatable :: A(:,:), AP(:), work(:)
    real(ep) :: rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('dtpcon', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_matrix_quad(n, n, A, seed = 130001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + real(2*n, ep); end do
        allocate(AP(np), work(3*n), iwork(n))
        call pack_sym_packed_quad('U', n, A, AP)
        call dtpcon('1', 'U', 'N', n, AP, rcond_ref, work, iwork, info)
        call target_dtpcon('1', 'U', 'N', n, AP, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AP, work, iwork)
    end do
    call report_finalize()
end program test_dtpcon
