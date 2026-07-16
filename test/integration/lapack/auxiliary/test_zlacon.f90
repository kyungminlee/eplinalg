! zlacon: legacy 1-norm estimator (complex). Same SAVE-state pattern as
! dlacon — drive ref + target serially, compare final est.
program test_zlacon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zlacon
    use ref_quad_lapack, only: zlacon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, kase, iter
    complex(ep), allocatable :: A(:,:), V(:), X(:), Y(:)
    real(ep) :: est_r, est_g, err, tol
    character(len=48) :: label

    call report_init('zlacon', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 19751 + 47 * i)
        allocate(V(n), X(n), Y(n))
        kase = 0; est_r = 0.0_ep
        do iter = 1, 20
            call zlacon(n, V, X, est_r, kase)
            if (kase == 0) exit
            if (kase == 1) then
                Y = matmul(A, X); X = Y
            else
                Y = matmul(transpose(conjg(A)), X); X = Y
            end if
        end do
        kase = 0; est_g = 0.0_ep
        do iter = 1, 20
            call target_zlacon(n, V, X, est_g, kase)
            if (kase == 0) exit
            if (kase == 1) then
                Y = matmul(A, X); X = Y
            else
                Y = matmul(transpose(conjg(A)), X); X = Y
            end if
        end do
        err = rel_err_scalar(est_g, est_r)
        tol = 64.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, V, X, Y)
    end do
    call report_finalize()
end program test_zlacon
