! zgesvdq: SVD with pivoted QR preconditioning (complex).
program test_zgesvdq
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgesvdq
    use ref_quad_lapack, only: zgesvdq
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [8,  16]
    integer :: i, m, n, info, numrank_r, numrank_g, lcwork, lrwork, liwork
    complex(ep), allocatable :: A0(:,:), A_r(:,:), A_g(:,:), U(:,:), V(:,:), cwork(:)
    real(ep),    allocatable :: S_r(:), S_g(:), rwork(:)
    integer, allocatable :: iwork(:)
    ! zgesvdq workspace query writes CWORK(1) and CWORK(2); allocate at least 2.
    complex(ep) :: wopt(2)
    real(ep) :: ropt(1), err, tol
    character(len=48) :: label

    call report_init('zgesvdq', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A0, seed = 25051 + 71 * i)
        allocate(A_r(m, n), A_g(m, n), S_r(min(m, n)), S_g(min(m, n)), U(m, n), V(n, n))
        liwork = m + n
        allocate(iwork(liwork))
        A_r = A0; A_g = A0
        call zgesvdq('A', 'N', 'N', 'N', 'N', m, n, A_r, m, S_r, U, m, V, n, &
                     numrank_r, iwork, -1, wopt, -1, ropt, -1, info)
        lcwork = max(1, int(real(wopt(1), ep)))
        lrwork = max(1, int(ropt(1)))
        liwork = max(liwork, iwork(1))
        deallocate(iwork); allocate(iwork(liwork), cwork(lcwork), rwork(lrwork))
        call zgesvdq('A', 'N', 'N', 'N', 'N', m, n, A_r, m, S_r, U, m, V, n, &
                     numrank_r, iwork, liwork, cwork, lcwork, rwork, lrwork, info)
        deallocate(cwork, rwork)
        call target_zgesvdq('A', 'N', 'N', 'N', 'N', m, n, A_g, m, S_g, U, m, V, n, &
                            numrank_g, info)
        err = max_rel_err_vec(S_g, S_r)
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ', n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_r, A_g, S_r, S_g, U, V, iwork)
    end do
    call report_finalize()
end program test_zgesvdq
