! dlacn2: 1-norm estimator for op(A) using reverse-communication
! interface. Compare estimates given the same A and same iteration.
program test_dlacn2
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dlacn2
    use ref_quad_lapack, only: dlacn2
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, kase_r, kase_g, isave_r(3), isave_g(3), iter
    real(ep), allocatable :: A(:,:), V_r(:), V_g(:), X_r(:), X_g(:), Y(:)
    integer, allocatable :: isgn_r(:), isgn_g(:)
    real(ep) :: est_r, est_g, err, tol
    character(len=48) :: label

    call report_init('dlacn2', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 19601 + 47 * i)
        allocate(V_r(n), V_g(n), X_r(n), X_g(n), Y(n))
        allocate(isgn_r(n), isgn_g(n))
        kase_r = 0; kase_g = 0
        isave_r = 0; isave_g = 0
        est_r = 0.0_ep; est_g = 0.0_ep
        ! Reverse comm loop — call until kase=0. Interleave ref + target.
        do iter = 1, 10
            call dlacn2(n, V_r, X_r, isgn_r, est_r, kase_r, isave_r)
            call target_dlacn2(n, V_g, X_g, isgn_g, est_g, kase_g, isave_g)
            if (kase_r == 0) exit
            if (kase_g == 0) exit
            ! kase=1: compute Y = A * X
            ! kase=2: compute Y = A^T * X
            if (kase_r == 1) then
                Y = matmul(A, X_r); X_r = Y
            else
                Y = matmul(transpose(A), X_r); X_r = Y
            end if
            if (kase_g == 1) then
                Y = matmul(A, X_g); X_g = Y
            else
                Y = matmul(transpose(A), X_g); X_g = Y
            end if
        end do
        err = rel_err_scalar(est_g, est_r)
        tol = 64.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, V_r, V_g, X_r, X_g, Y, isgn_r, isgn_g)
    end do
    call report_finalize()
end program test_dlacn2
