program test_dormtr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_symmetric_matrix_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dormtr
    use ref_quad_lapack, only: dsytrd, dormtr
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: ncols = 5
    integer :: i, n, info, lwork
    real(ep), allocatable :: A(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep), allocatable :: D(:), E(:), tau(:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dormtr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A, seed = 230021 + 47 * i)
        call gen_matrix_quad(n, ncols, C0, seed = 230031 + 47 * i)
        allocate(D(n), E(n-1), tau(n-1))
        call dsytrd('U', n, A, n, D, E, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dsytrd('U', n, A, n, D, E, tau, work, lwork, info)
        deallocate(work)
        allocate(C_ref(n, ncols), C_got(n, ncols))
        C_ref = C0; C_got = C0
        call dormtr('L', 'U', 'N', n, ncols, A, n, tau, C_ref, n, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dormtr('L', 'U', 'N', n, ncols, A, n, tau, C_ref, n, work, lwork, info)
        call target_dormtr('L', 'U', 'N', n, ncols, A, n, tau, C_got, n, info)
        err = max_rel_err_mat(C_got, C_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, C0, D, E, tau, work, C_ref, C_got)
    end do
    call report_finalize()
end program test_dormtr
