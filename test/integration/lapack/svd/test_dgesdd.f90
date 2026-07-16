program test_dgesdd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgesdd
    use ref_quad_lapack, only: dgesdd
    implicit none

    integer, parameter :: ms(*) = [12, 32, 48]
    integer, parameter :: ns(*) = [8,  24, 40]
    integer :: i, m, n, kmn, info, lwork, j
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: s_ref(:), s_got(:), U_ref(:,:), U_got(:,:), VT_ref(:,:), VT_got(:,:), work(:)
    real(ep), allocatable :: US(:,:), USVt(:,:)
    integer,  allocatable :: iwork(:)
    real(ep) :: wopt(1), err_s, err_r, tol, anorm, rnorm
    character(len=48) :: label

    call report_init('dgesdd', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); kmn = min(m, n)
        call gen_matrix_quad(m, n, A0, seed = 21001 + 47 * i)
        allocate(A_ref(m,n), A_got(m,n), s_ref(kmn), s_got(kmn))
        allocate(U_ref(m,m), U_got(m,m), VT_ref(n,n), VT_got(n,n), iwork(8*kmn))
        A_ref = A0; A_got = A0
        call dgesdd('A', m, n, A_ref, m, s_ref, U_ref, m, VT_ref, n, &
                    wopt, -1, iwork, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgesdd('A', m, n, A_ref, m, s_ref, U_ref, m, VT_ref, n, &
                    work, lwork, iwork, info)
        deallocate(work)
        call target_dgesdd('A', m, n, A_got, m, s_got, U_got, m, VT_got, n, info)
        err_s = max_rel_err_vec(s_got, s_ref)
        allocate(US(m,kmn), USVt(m,n))
        do j = 1, kmn
            US(:, j) = U_got(:, j) * s_got(j)
        end do
        USVt = matmul(US, VT_got(1:kmn, 1:n))
        anorm = maxval(abs(A0))
        rnorm = maxval(abs(A0 - USVt))
        err_r = rnorm / max(anorm, tiny(1.0_ep))
        tol   = 16.0_ep * real(max(m,n), ep)**3 * target_eps
        write(label, '(a,i0,a,i0,a)') 'm=', m, ',n=', n, ',out=S'
        call report_case(trim(label), err_s, tol)
        write(label, '(a,i0,a,i0,a)') 'm=', m, ',n=', n, ',out=residual'
        call report_case(trim(label), err_r, tol)
        deallocate(A_ref, A_got, s_ref, s_got, U_ref, U_got, VT_ref, VT_got, iwork, US, USVt)
    end do
    call report_finalize()
end program test_dgesdd
