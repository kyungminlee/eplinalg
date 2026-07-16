program test_zgesvd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgesvd
    use ref_quad_lapack, only: zgesvd
    implicit none

    integer, parameter :: ms(*) = [12, 32, 48]
    integer, parameter :: ns(*) = [8,  24, 40]
    integer :: i, m, n, kmn, info, lwork, lrwork, j
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: U_ref(:,:), U_got(:,:), VT_ref(:,:), VT_got(:,:)
    complex(ep), allocatable :: work_ref(:), US(:,:), USVt(:,:)
    real(ep),    allocatable :: s_ref(:), s_got(:), rwork(:)
    complex(ep) :: wopt(1)
    real(ep)    :: err_s, err_r, tol, anorm, rnorm
    character(len=48) :: label

    call report_init('zgesvd', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); kmn = min(m, n)
        call gen_matrix_complex(m, n, A0, seed = 12101 + 79 * i)

        allocate(A_ref(m, n), A_got(m, n), s_ref(kmn), s_got(kmn))
        allocate(U_ref(m, m), U_got(m, m), VT_ref(n, n), VT_got(n, n))
        lrwork = max(1, 5 * kmn)
        allocate(rwork(lrwork))
        A_ref = A0; A_got = A0

        call zgesvd('A', 'A', m, n, A_ref, m, s_ref, U_ref, m, VT_ref, n, &
                    wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work_ref(lwork))
        call zgesvd('A', 'A', m, n, A_ref, m, s_ref, U_ref, m, VT_ref, n, &
                    work_ref, lwork, rwork, info)

        call target_zgesvd('A', 'A', m, n, A_got, m, s_got, U_got, m, &
                           VT_got, n, info)

        err_s = max_rel_err_vec(s_got, s_ref)
        allocate(US(m, kmn), USVt(m, n))
        do j = 1, kmn
            US(:, j) = U_got(:, j) * s_got(j)
        end do
        USVt = matmul(US, VT_got(1:kmn, 1:n))
        anorm = maxval(abs(A0))
        rnorm = maxval(abs(A0 - USVt))
        err_r = rnorm / max(anorm, tiny(1.0_ep))
        tol   = 16.0_ep * real(max(m, n), ep)**3 * target_eps

        write(label, '(a,i0,a,i0,a)') 'm=', m, ',n=', n, ',out=S'
        call report_case(trim(label), err_s, tol)
        write(label, '(a,i0,a,i0,a)') 'm=', m, ',n=', n, ',out=residual'
        call report_case(trim(label), err_r, tol)

        deallocate(A_ref, A_got, s_ref, s_got, U_ref, U_got, VT_ref, VT_got, &
                   work_ref, rwork, US, USVt)
    end do
    call report_finalize()
end program test_zgesvd
