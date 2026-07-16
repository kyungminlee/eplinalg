! dgetsls solves overdetermined LS via TS-QR. m > n, TRANS='N'.
program test_dgetsls
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgetsls
    use ref_quad_lapack, only: dgetsls
    implicit none

    integer, parameter :: ms(*)   = [16, 32, 48]
    integer, parameter :: ns(*)   = [12, 20, 28]
    integer, parameter :: nrhs    = 3
    integer :: i, m, n, info, lwork
    real(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgetsls', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 390001 + 47 * i)
        call gen_matrix_quad(m, nrhs, B0, seed = 390011 + 47 * i)
        allocate(A_ref(m, n), B_ref(m, nrhs), A_got(m, n), B_got(m, nrhs))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        call dgetsls('N', m, n, nrhs, A_ref, m, B_ref, m, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgetsls('N', m, n, nrhs, A_ref, m, B_ref, m, work, lwork, info)
        call target_dgetsls('N', m, n, nrhs, A_got, m, B_got, m, info)
        ! Compare the LS solution X (top n rows of B), which is unique.
        err = max_rel_err_mat(B_got(1:n, 1:nrhs), B_ref(1:n, 1:nrhs))
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, A_got, B_got, work)
    end do
    call report_finalize()
end program test_dgetsls
