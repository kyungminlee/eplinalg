! dgemqr applies Q from dgeqr output. Run reference dgeqr first to
! produce A and T, then apply Q to a random C from the left.
program test_dgemqr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgemqr
    use ref_quad_lapack, only: dgeqr, dgemqr
    implicit none

    integer, parameter :: ms(*) = [16, 32, 48]
    integer, parameter :: ns(*) = [12, 20, 28]
    integer, parameter :: ncols = 5
    integer :: i, m, n, info, lwork, tsize
    real(ep), allocatable :: A(:,:), T(:), work(:)
    real(ep), allocatable :: C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: tquery(5), wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgemqr', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A, seed = 370041 + 47 * i)
        call dgeqr(m, n, A, m, tquery, -1, wopt, -1, info)
        tsize = max(1, int(tquery(1)))
        lwork = max(1, int(wopt(1)))
        allocate(T(tsize), work(lwork))
        call dgeqr(m, n, A, m, T, tsize, work, lwork, info)
        deallocate(work)
        call gen_matrix_quad(m, ncols, C0, seed = 370051 + 47 * i)
        allocate(C_ref(m, ncols), C_got(m, ncols))
        C_ref = C0; C_got = C0
        call dgemqr('L', 'N', m, ncols, n, A, m, T, tsize, C_ref, m, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgemqr('L', 'N', m, ncols, n, A, m, T, tsize, C_ref, m, work, lwork, info)
        call target_dgemqr('L', 'N', m, ncols, n, A, m, T, tsize, C_got, m, info)
        err = max_rel_err_mat(C_got, C_ref)
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, T, C0, C_ref, C_got, work)
    end do
    call report_finalize()
end program test_dgemqr
