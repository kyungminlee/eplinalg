! dtzrzf: RZ factorization of a trapezoidal matrix. Compare only the
! triangular factor R = A(1:m,1:m) since the reflector storage in
! A(:,m+1:n) and tau is not canonical between blocking heuristics.
program test_dtzrzf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtzrzf
    use ref_quad_lapack, only: dtzrzf
    implicit none

    integer, parameter :: ms(*) = [6, 16]
    integer, parameter :: ns(*) = [12, 32]
    integer :: i, m, n, info, lwork, j
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: tau_ref(:), tau_got(:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dtzrzf', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 20101 + 83 * i)
        ! Pre-shape into trapezoidal: zero below diagonal of leading m x m
        do j = 1, m-1
            A0(j+1:m, j) = 0.0_ep
        end do
        allocate(A_ref(m,n), A_got(m,n), tau_ref(m), tau_got(m))
        A_ref = A0; A_got = A0
        call dtzrzf(m, n, A_ref, m, tau_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dtzrzf(m, n, A_ref, m, tau_ref, work, lwork, info)
        deallocate(work)
        call target_dtzrzf(m, n, A_got, m, tau_got, info)
        ! Compare only upper triangle of the leading m x m block.
        do j = 1, m
            A_ref(j+1:m, j) = 0.0_ep
            A_got(j+1:m, j) = 0.0_ep
        end do
        err = max_rel_err_mat(A_got(1:m,1:m), A_ref(1:m,1:m))
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_dtzrzf
