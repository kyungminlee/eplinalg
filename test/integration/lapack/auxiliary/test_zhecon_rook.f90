program test_zhecon_rook
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zhetrf_rook, target_zhecon_rook
    use ref_quad_lapack, only: zhetrf_rook, zhecon_rook
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info, lwork, j
    integer, allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:), cwork(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    complex(ep) :: wopt(1)
    character(len=48) :: label

    call report_init('zhecon_rook', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 144101 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + cmplx(real(n, ep), 0.0_ep, ep); end do
        anorm = maxval(sum(abs(A0), dim=1))
        allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
        A_ref = A0; A_got = A0
        call zhetrf_rook('U', n, A_ref, n, ipiv_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zhetrf_rook('U', n, A_ref, n, ipiv_ref, work, lwork, info)
        deallocate(work)
        call target_zhetrf_rook('U', n, A_got, n, ipiv_got, info)
        allocate(cwork(2*n))
        call zhecon_rook('U', n, A_ref, n, ipiv_ref, anorm, rcond_ref, cwork, info)
        call target_zhecon_rook('U', n, A_got, n, ipiv_got, anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 100.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, ipiv_ref, ipiv_got, cwork)
    end do
    call report_finalize()
end program test_zhecon_rook
