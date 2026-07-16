! dgedmd: Dynamic Mode Decomposition (real). Each invocation extracts
! a Ritz spectrum from snapshot pair (X, Y = A*X). Layer (a) compares
! the canonicalized complex spectrum (sorted lex by Re,Im); layer (b)
! checks the routine's own residual output stays within tol_res*||A||_F.
! Layer (c) (full reconstruction) is intentionally omitted in v1 — it
! requires a quad pinv(Z) solve and is gated by the same JOBZ='V'/JOBF='N'
! combo already covered indirectly by layer (b).
program test_dgedmd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec, max_rel_err_vec_z, sort_eig_lex_z
    use test_data,       only: gen_orthogonal_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgedmd
    use ref_quad_lapack, only: dgedmd
    implicit none

    integer, parameter :: nsizes = 3
    integer, parameter :: ms(nsizes) = [8, 12, 24]
    integer, parameter :: ns(nsizes) = [5, 8, 16]
    character(len=1), parameter :: jobs_set(2) = ['N', 'S']
    character(len=1), parameter :: jobz_set(2) = ['N', 'V']
    character(len=1), parameter :: jobf_set(2) = ['N', 'E']

    integer :: isz, ijobs, ijobz, ijobr, ijobf
    integer :: m, n, info_ref, info_got, k_ref, k_got
    integer :: lwork_ref, liwork_ref, j
    integer, parameter :: whtsvd = 1, nrnk = -1
    real(ep) :: tol_dmd, tol_eig, tol_res, normA, err_a, err_b
    real(ep) :: wopt(2)
    integer :: iwopt(2)
    real(ep), allocatable :: A_quad(:,:), X0(:,:), Y0(:,:)
    real(ep), allocatable :: U(:,:), lam_re(:), lam_im(:)
    real(ep), allocatable :: X_ref(:,:), Y_ref(:,:), X_got(:,:), Y_got(:,:)
    real(ep), allocatable :: reig_ref(:), imeig_ref(:), reig_got(:), imeig_got(:)
    real(ep), allocatable :: Z_ref(:,:), Z_got(:,:), B_ref(:,:), B_got(:,:)
    real(ep), allocatable :: W_ref(:,:), W_got(:,:), S_ref(:,:), S_got(:,:)
    real(ep), allocatable :: res_ref(:), res_got(:), work(:)
    integer,  allocatable :: iwork(:)
    complex(ep), allocatable :: spec_ref(:), spec_got(:)
    character(len=1) :: jobs, jobz, jobr, jobf
    character(len=80) :: label

    call report_init('dgedmd', target_name)

    do isz = 1, nsizes
        m = ms(isz); n = ns(isz)
        ! A = U * diag(lam) * U^T (orthogonal U so U^-1 = U^T).
        ! Spectrum mixes one near-zero, one negative, and one
        ! complex-conjugate pair — encoded via the 2x2 block trick.
        call gen_orthogonal_quad(n, U, seed = 91 + 17*isz)
        allocate(lam_re(n), lam_im(n))
        do j = 1, n
            lam_re(j) = 0.5_ep + 0.1_ep * real(j, ep)
            lam_im(j) = 0.0_ep
        end do
        if (n >= 4) then
            lam_re(1) = 0.6_ep;   lam_im(1) =  0.2_ep
            lam_re(2) = 0.6_ep;   lam_im(2) = -0.2_ep
            lam_re(n) = 1.0e-3_ep
        end if
        allocate(A_quad(n, n))
        ! A = U * diag(lam) * U^T (real Schur for the pair: a 2x2 block
        ! [[Re, Im],[-Im, Re]] is a rotation with eigenvalues Re ± i*Im).
        block
            real(ep), allocatable :: T(:,:)
            integer :: ii, jj
            allocate(T(n, n)); T = 0.0_ep
            ii = 1
            do while (ii <= n)
                if (ii < n .and. lam_im(ii) /= 0.0_ep) then
                    T(ii,   ii  ) = lam_re(ii)
                    T(ii,   ii+1) = lam_im(ii)
                    T(ii+1, ii  ) = -lam_im(ii)
                    T(ii+1, ii+1) = lam_re(ii)
                    ii = ii + 2
                else
                    T(ii, ii) = lam_re(ii)
                    ii = ii + 1
                end if
            end do
            A_quad = matmul(matmul(U, T), transpose(U))
            deallocate(T)
        end block
        normA = sqrt(sum(A_quad**2))

        call gen_matrix_quad(m, n, X0, seed = 401 + 23*isz)
        ! Y = X * A^T  (so for each row i: Y(i,:) = A * X(i,:)^T).
        ! Actually m >= n; here treat X as (m x n) snapshot matrix and
        ! Y = X * A^T satisfies Y(i, j) = sum_k X(i, k) * A(j, k) — but
        ! the DMD convention is X[:,k] is a snapshot vector at time k.
        ! Use: m = state dim, n = number of snapshots. A is n-by-n.
        ! X0 is m-by-n; let A_state operate as Y = X * A^T so that columns
        ! of Y = columns of X advanced by A^T (n-by-n). Spectrum recovered
        ! is sigma(A^T) = sigma(A) (real spectrum).
        allocate(Y0(m, n))
        Y0 = matmul(X0, transpose(A_quad))

        do ijobs = 1, 2
            jobs = jobs_set(ijobs)
            do ijobz = 1, 2
                jobz = jobz_set(ijobz)
                do ijobr = 1, 2
                    if (ijobr == 1) then
                        jobr = 'N'
                    else
                        ! Doc: JOBR='R' requires JOBZ='V'.
                        if (jobz == 'V') then
                            jobr = 'R'
                        else
                            cycle
                        end if
                    end if
                    do ijobf = 1, 2
                        jobf = jobf_set(ijobf)
                        ! Allocate per-combo working arrays
                        allocate(X_ref(m, n), Y_ref(m, n), X_got(m, n), Y_got(m, n))
                        allocate(reig_ref(n), imeig_ref(n), reig_got(n), imeig_got(n))
                        allocate(Z_ref(m, n), Z_got(m, n))
                        allocate(B_ref(m, n), B_got(m, n))
                        allocate(W_ref(n, n), W_got(n, n))
                        allocate(S_ref(n, n), S_got(n, n))
                        allocate(res_ref(n), res_got(n))
                        X_ref = X0; Y_ref = Y0
                        X_got = X0; Y_got = Y0
                        tol_dmd = 10.0_ep * target_eps
                        ! Reference workspace query
                        call dgedmd(jobs, jobz, jobr, jobf, whtsvd, m, n,    &
                                    X_ref, m, Y_ref, m, nrnk, tol_dmd,        &
                                    k_ref, reig_ref, imeig_ref, Z_ref, m,     &
                                    res_ref, B_ref, m, W_ref, n, S_ref, n,    &
                                    wopt, -1, iwopt, -1, info_ref)
                        if (info_ref == 0) then
                            lwork_ref  = max(1, int(wopt(1)))
                            liwork_ref = max(1, iwopt(1))
                            allocate(work(lwork_ref), iwork(liwork_ref))
                            call dgedmd(jobs, jobz, jobr, jobf, whtsvd, m, n,    &
                                        X_ref, m, Y_ref, m, nrnk, tol_dmd,       &
                                        k_ref, reig_ref, imeig_ref, Z_ref, m,    &
                                        res_ref, B_ref, m, W_ref, n, S_ref, n,   &
                                        work, lwork_ref, iwork, liwork_ref,      &
                                        info_ref)
                            deallocate(work, iwork)
                        end if
                        call target_dgedmd(jobs, jobz, jobr, jobf, whtsvd, m, n, &
                                           X_got, m, Y_got, m, nrnk, tol_dmd,    &
                                           k_got, reig_got, imeig_got, Z_got, m, &
                                           res_got, B_got, m, W_got, n,          &
                                           S_got, n, info_got)
                        write(label, '(a,i0,a,i0,a,a,a,a,a,a,a,a)') &
                            'm=', m, ',n=', n,                                &
                            ',JOBS=', jobs, ',JOBZ=', jobz,                    &
                            ',JOBR=', jobr, ',JOBF=', jobf

                        if (info_ref /= 0 .or. info_got /= 0 .or. k_ref /= k_got) then
                            ! Surface as a layer-(a) failure with err=1
                            tol_eig = 100.0_ep * real(n, ep)**2 * target_eps
                            call report_case(trim(label) // ',layer=a,info-mismatch', &
                                              1.0_ep, tol_eig)
                        else if (k_got > 0) then
                            ! Layer (a) — sorted spectrum
                            allocate(spec_ref(k_got), spec_got(k_got))
                            do j = 1, k_got
                                spec_ref(j) = cmplx(reig_ref(j), imeig_ref(j), ep)
                                spec_got(j) = cmplx(reig_got(j), imeig_got(j), ep)
                            end do
                            call sort_eig_lex_z(spec_ref)
                            call sort_eig_lex_z(spec_got)
                            err_a = max_rel_err_vec_z(spec_got, spec_ref)
                            tol_eig = 100.0_ep * real(n, ep)**2 * target_eps
                            call report_case(trim(label) // ',layer=a', err_a, tol_eig)
                            deallocate(spec_ref, spec_got)
                            ! Layer (b) — residuals (only when JOBR='R').
                            ! Absolute delta scaled by normA — robust to
                            ! near-zero residuals (where pure-relative err
                            ! collapses) and well-conditioned data alike.
                            if (jobr == 'R') then
                                tol_res = 100.0_ep * real(n, ep) * target_eps
                                err_b = maxval(abs(res_got(1:k_got) -      &
                                                    res_ref(1:k_got))) /  &
                                        max(normA, target_eps)
                                call report_case(trim(label) // ',layer=b', err_b, tol_res)
                            end if
                        end if
                        deallocate(X_ref, Y_ref, X_got, Y_got)
                        deallocate(reig_ref, imeig_ref, reig_got, imeig_got)
                        deallocate(Z_ref, Z_got, B_ref, B_got)
                        deallocate(W_ref, W_got, S_ref, S_got)
                        deallocate(res_ref, res_got)
                    end do
                end do
            end do
        end do

        deallocate(A_quad, U, lam_re, lam_im, X0, Y0)
    end do

    call report_finalize()
end program test_dgedmd
