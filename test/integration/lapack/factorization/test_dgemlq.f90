program test_dgemlq
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgemlq
    use ref_quad_lapack, only: dgelq, dgemlq
    implicit none

    integer, parameter :: ms(*) = [12, 20, 28]
    integer, parameter :: ns(*) = [16, 32, 48]
    integer, parameter :: ncols = 5
    integer :: i, m, n, info, lwork, tsize
    real(ep), allocatable :: A(:,:), T(:), work(:)
    real(ep), allocatable :: C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: tquery(5), wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgemlq', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A, seed = 370081 + 47 * i)
        call dgelq(m, n, A, m, tquery, -1, wopt, -1, info)
        tsize = max(1, int(tquery(1)))
        lwork = max(1, int(wopt(1)))
        allocate(T(tsize), work(lwork))
        call dgelq(m, n, A, m, T, tsize, work, lwork, info)
        deallocate(work)
        call gen_matrix_quad(n, ncols, C0, seed = 370091 + 47 * i)
        allocate(C_ref(n, ncols), C_got(n, ncols))
        C_ref = C0; C_got = C0
        call dgemlq('L', 'N', n, ncols, m, A, m, T, tsize, C_ref, n, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgemlq('L', 'N', n, ncols, m, A, m, T, tsize, C_ref, n, work, lwork, info)
        call target_dgemlq('L', 'N', n, ncols, m, A, m, T, tsize, C_got, n, info)
        err = max_rel_err_mat(C_got, C_ref)
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, T, C0, C_ref, C_got, work)
    end do
    call report_finalize()
end program test_dgemlq
