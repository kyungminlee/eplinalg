program test_dgecon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgetrf, target_dgecon
    use ref_quad_lapack, only: dgetrf, dgecon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info
    integer, allocatable :: ipiv_ref(:), ipiv_got(:), iwork(:)
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('dgecon', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A0, seed = 110001 + 47 * i)
        anorm = maxval(sum(abs(A0), dim=1))
        allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
        allocate(work(4*n), iwork(n))
        A_ref = A0; A_got = A0
        call dgetrf(n, n, A_ref, n, ipiv_ref, info)
        call target_dgetrf(n, n, A_got, n, ipiv_got, info)
        call dgecon('1', n, A_ref, n, anorm, rcond_ref, work, iwork, info)
        call target_dgecon('1', n, A_got, n, anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 100.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, ipiv_ref, ipiv_got, work, iwork)
    end do
    call report_finalize()
end program test_dgecon
