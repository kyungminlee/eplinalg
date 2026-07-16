! Compares eigenvalue vector W (unique, ascending-sorted).
! Eigenvectors Z have a sign convention that's precision-dependent
! (a near-zero sign comparison inside LAPACK can flip), so we
! validate them via the invariant ||A*Z - Z*diag(W)|| / ||A|| which
! is immune to column-sign flips.
program test_dsyev
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsyev
    use ref_quad_lapack, only: dsyev
    implicit none

    integer, parameter :: ns(*) = [8, 32, 64]
    ! Cycle JOBZ and UPLO so each shape exercises an independent path:
    ! ('V','U') drives the back-transformation through DORMTR and the
    ! upper tridiagonal reduction; ('V','L') the lower variant; ('N','U')
    ! skips the back-transformation entirely (eigenvalues only).
    character(len=1), parameter :: jobzs(*) = ['V', 'V', 'N']
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n, j, info, lwork
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: w_ref(:), w_got(:), work_ref(:)
    real(ep), allocatable :: AZ(:,:), ZD(:,:)
    real(ep) :: wopt(1), err_w, err_r, tol, anorm, rnorm
    character(len=48) :: label

    call report_init('dsyev', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 10001 + 71 * i)

        allocate(A_ref(n, n), A_got(n, n), w_ref(n), w_got(n))
        A_ref = A0;  A_got = A0

        call dsyev(jobzs(i), uplos(i), n, A_ref, n, w_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work_ref(lwork))
        call dsyev(jobzs(i), uplos(i), n, A_ref, n, w_ref, &
                   work_ref, lwork, info)

        call target_dsyev(jobzs(i), uplos(i), n, A_got, n, w_got, info)

        err_w = max_rel_err_vec(w_got, w_ref)
        tol   = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,a,a,a,a,i0,a)') &
            'jobz=', jobzs(i), ',uplo=', uplos(i), ',n=', n, ',out=W'
        call report_case(trim(label), err_w, tol)

        ! Invariant residual ||A0 * Z_got - Z_got * diag(w_got)|| only
        ! makes sense when eigenvectors were requested.
        if (jobzs(i) == 'V') then
            allocate(AZ(n, n), ZD(n, n))
            AZ = matmul(A0, A_got)
            do j = 1, n
                ZD(:, j) = A_got(:, j) * w_got(j)
            end do
            anorm = maxval(abs(A0))
            rnorm = maxval(abs(AZ - ZD))
            err_r = rnorm / max(anorm, tiny(1.0_ep))
            write(label, '(a,a,a,a,a,i0,a)') &
                'jobz=', jobzs(i), ',uplo=', uplos(i), ',n=', n, &
                ',out=residual'
            call report_case(trim(label), err_r, tol)
            deallocate(AZ, ZD)
        end if

        deallocate(A_ref, A_got, w_ref, w_got, work_ref)
    end do
    call report_finalize()
end program test_dsyev
