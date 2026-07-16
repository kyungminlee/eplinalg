! Compares Hermitian eigenvalues W; validates eigenvectors via the
! invariant ||A*Z - Z*diag(W)|| / ||A|| (sign-flip immune).
program test_zheev
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zheev
    use ref_quad_lapack, only: zheev
    implicit none

    integer, parameter :: ns(*) = [8, 24, 48]
    integer :: i, n, j, info, lwork, lrwork
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: work_ref(:), AZ(:,:), ZD(:,:)
    real(ep),    allocatable :: w_ref(:), w_got(:), rwork_ref(:)
    complex(ep) :: wopt(1)
    real(ep) :: err_w, err_r, tol, anorm, rnorm
    character(len=48) :: label

    call report_init('zheev', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 11001 + 73 * i)

        allocate(A_ref(n, n), A_got(n, n), w_ref(n), w_got(n))
        lrwork = max(1, 3*n - 2)
        allocate(rwork_ref(lrwork))
        A_ref = A0;  A_got = A0

        call zheev('V', 'U', n, A_ref, n, w_ref, wopt, -1, rwork_ref, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work_ref(lwork))
        call zheev('V', 'U', n, A_ref, n, w_ref, work_ref, lwork, rwork_ref, info)

        call target_zheev('V', 'U', n, A_got, n, w_got, info)

        err_w = max_rel_err_vec(w_got, w_ref)
        allocate(AZ(n, n), ZD(n, n))
        AZ = matmul(A0, A_got)
        do j = 1, n
            ZD(:, j) = A_got(:, j) * cmplx(w_got(j), 0.0_ep, ep)
        end do
        anorm = maxval(abs(A0))
        rnorm = maxval(abs(AZ - ZD))
        err_r = rnorm / max(anorm, tiny(1.0_ep))
        tol   = 16.0_ep * real(n, ep)**3 * target_eps

        write(label, '(a,i0,a)') 'n=', n, ',out=W'
        call report_case(trim(label), err_w, tol)
        write(label, '(a,i0,a)') 'n=', n, ',out=residual'
        call report_case(trim(label), err_r, tol)

        deallocate(A_ref, A_got, w_ref, w_got, work_ref, rwork_ref, AZ, ZD)
    end do
    call report_finalize()
end program test_zheev
