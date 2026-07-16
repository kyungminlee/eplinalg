program test_zhecon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zhetrf, target_zhecon
    use ref_quad_lapack, only: zhetrf, zhecon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info, lwork
    integer, allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:), cwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('zhecon', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 119001 + 47 * i)
        anorm = maxval(sum(abs(A0), dim=1))
        allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
        A_ref = A0; A_got = A0
        call zhetrf('U', n, A_ref, n, ipiv_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zhetrf('U', n, A_ref, n, ipiv_ref, work, lwork, info)
        deallocate(work)
        call target_zhetrf('U', n, A_got, n, ipiv_got, info)
        allocate(cwork(2*n))
        call zhecon('U', n, A_ref, n, ipiv_ref, anorm, rcond_ref, cwork, info)
        call target_zhecon('U', n, A_got, n, ipiv_got, anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, ipiv_ref, ipiv_got, cwork)
    end do
    call report_finalize()
end program test_zhecon
