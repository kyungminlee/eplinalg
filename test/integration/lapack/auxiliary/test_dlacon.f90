! dlacon: legacy 1-norm estimator for op(A) using reverse-communication.
! Compare estimates given the same A. Note: dlacon stores its iteration
! state in internal SAVE variables (not via an ISAVE arg), so the ref and
! target paths each carry their own SAVE state — they cannot be
! interleaved within the same process. We therefore drive ref to
! completion, then target, and compare the final est value.
program test_dlacon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dlacon
    use ref_quad_lapack, only: dlacon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, kase, iter
    real(ep), allocatable :: A(:,:), V(:), X(:), Y(:)
    integer, allocatable :: isgn(:)
    real(ep) :: est_r, est_g, err, tol
    character(len=48) :: label

    call report_init('dlacon', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 19651 + 47 * i)
        allocate(V(n), X(n), Y(n), isgn(n))
        ! Reference path.
        kase = 0; est_r = 0.0_ep
        do iter = 1, 20
            call dlacon(n, V, X, isgn, est_r, kase)
            if (kase == 0) exit
            if (kase == 1) then
                Y = matmul(A, X); X = Y
            else
                Y = matmul(transpose(A), X); X = Y
            end if
        end do
        ! Target path (independent SAVE state).
        kase = 0; est_g = 0.0_ep
        do iter = 1, 20
            call target_dlacon(n, V, X, isgn, est_g, kase)
            if (kase == 0) exit
            if (kase == 1) then
                Y = matmul(A, X); X = Y
            else
                Y = matmul(transpose(A), X); X = Y
            end if
        end do
        err = rel_err_scalar(est_g, est_r)
        tol = 64.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, V, X, Y, isgn)
    end do
    call report_finalize()
end program test_dlacon
