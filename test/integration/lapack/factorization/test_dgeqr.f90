! dgeqr: workspace-driven blocked QR. Compares the full overwritten A
! and T(1:3) (TSIZE/MB/NB metadata) plus T(6:TSIZE) (DLATSQR/DGEQRT
! block data). T(4) and T(5) are unused per the spec and are skipped.
program test_dgeqr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgeqr
    use ref_quad_lapack, only: dgeqr
    implicit none

    integer, parameter :: ms(*) = [16, 32, 48]
    integer, parameter :: ns(*) = [12, 20, 28]
    integer :: i, m, n, info, lwork, tsize, j
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:)
    real(ep), allocatable :: T_ref(:), T_got(:)
    real(ep) :: tquery(5), wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgeqr', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 370001 + 47 * i)
        allocate(A_ref(m, n), A_got(m, n))
        A_ref = A0; A_got = A0
        call dgeqr(m, n, A_ref, m, tquery, -1, wopt, -1, info)
        tsize = max(1, int(tquery(1)))
        lwork = max(1, int(wopt(1)))
        allocate(T_ref(tsize), T_got(tsize), work(lwork))
        call dgeqr(m, n, A_ref, m, T_ref, tsize, work, lwork, info)
        call target_dgeqr(m, n, A_got, m, T_got, tsize, info)
        err = max_rel_err_mat(A_got, A_ref)
        do j = 1, 3
            err = max(err, abs(T_got(j) - T_ref(j)) / max(abs(T_ref(j)), 1.0e-50_ep))
        end do
        do j = 6, tsize
            err = max(err, abs(T_got(j) - T_ref(j)) / max(abs(T_ref(j)), 1.0e-50_ep))
        end do
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, T_ref, T_got, work)
    end do
    call report_finalize()
end program test_dgeqr
