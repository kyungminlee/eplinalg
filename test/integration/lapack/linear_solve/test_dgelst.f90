! dgelst: T-factor variant of dgels (LAPACK 3.10+).
program test_dgelst
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgelst
    use ref_quad_lapack, only: dgelst
    implicit none

    integer, parameter :: ms(*) = [16, 24]
    integer, parameter :: ns(*) = [8, 16]
    integer, parameter :: nrhs  = 2
    integer :: i, m, n, ldb, info, lwork
    real(ep), allocatable :: A0(:,:), B0(:,:), A_r(:,:), B_r(:,:), A_g(:,:), B_g(:,:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgelst', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); ldb = max(m, n)
        call gen_matrix_quad(m, n,    A0, seed = 153001 + 47 * i)
        call gen_matrix_quad(ldb, nrhs, B0, seed = 153011 + 47 * i)
        allocate(A_r(m,n), B_r(ldb,nrhs), A_g(m,n), B_g(ldb,nrhs))
        A_r = A0; B_r = B0; A_g = A0; B_g = B0
        call dgelst('N', m, n, nrhs, A_r, m, B_r, ldb, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgelst('N', m, n, nrhs, A_r, m, B_r, ldb, work, lwork, info)
        deallocate(work)
        call target_dgelst('N', m, n, nrhs, A_g, m, B_g, ldb, info)
        err = max_rel_err_mat(B_g(1:n,:), B_r(1:n,:))
        tol = 16.0_ep * real(max(m,n), ep)**3 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_r, B_r, A_g, B_g)
    end do
    call report_finalize()
end program test_dgelst
