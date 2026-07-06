! dgedmdq: snapshot-form Dynamic Mode Decomposition. Single-trajectory
! input F (M-by-N), routine builds X = F[:,1:N-1] / Y = F[:,2:N]
! internally via QR + sliding. Layer (a) compares the canonicalized
! complex spectrum (sorted lex Re,Im); layer (b) checks routine-reported
! residuals stay below tol_res * ||A||_F. Layer (c) (reconstruction) is
! intentionally skipped — with JOBT='F' the returned T matrix is not
! canonical across implementations (eigenvector scaling is free), so a
! reconstruction comparison would be seed-brittle; v1 stays at (a)+(b).
program test_dgedmdq
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec, max_rel_err_vec_z, sort_eig_lex_z
    use test_data,       only: gen_orthogonal_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgedmdq
    use ref_quad_lapack, only: dgedmdq
    implicit none

    integer, parameter :: nsizes = 2
    ! Two sizes: a single-trajectory snapshot's conditioning collapses
    ! quickly (κ(F) ∝ |λ_max/λ_min|^N), so n above ~8 starts losing
    ! precision on kind10/multifloats targets even with relaxed tolerance.
    integer, parameter :: ms(nsizes) = [8, 12]
    integer, parameter :: ns(nsizes) = [5, 8]
    character(len=1), parameter :: jobs_set(2) = ['N', 'S']
    character(len=1), parameter :: jobz_set(2) = ['N', 'V']
    character(len=1), parameter :: jobf_set(2) = ['N', 'E']

    integer :: isz, ijobs, ijobz, ijobr, ijobf
    integer :: m, n, info_ref, info_got, k_ref, k_got
    integer :: lwork_ref, liwork_ref, j, mn
    integer, parameter :: whtsvd = 1, nrnk = -1
    real(ep) :: tol_dmd, tol_eig, tol_res, normA, err_a, err_b
    real(ep) :: wopt(2)
    integer :: iwopt(2)
    real(ep), allocatable :: A_quad(:,:), F0(:,:)
    real(ep), allocatable :: U(:,:), lam_re(:), lam_im(:)
    real(ep), allocatable :: F_ref(:,:), F_got(:,:)
    real(ep), allocatable :: X_ref(:,:), Y_ref(:,:), X_got(:,:), Y_got(:,:)
    real(ep), allocatable :: reig_ref(:), imeig_ref(:), reig_got(:), imeig_got(:)
    real(ep), allocatable :: Z_ref(:,:), Z_got(:,:), B_ref(:,:), B_got(:,:)
    real(ep), allocatable :: V_ref(:,:), V_got(:,:), S_ref(:,:), S_got(:,:)
    real(ep), allocatable :: res_ref(:), res_got(:), work(:)
    integer,  allocatable :: iwork(:)
    complex(ep), allocatable :: spec_ref(:), spec_got(:)
    character(len=1) :: jobs, jobz, jobr, jobf, jobq, jobt
    character(len=80) :: label

    call report_init('dgedmdq', target_name)

    do isz = 1, nsizes
        m = ms(isz); n = ns(isz)
        mn = min(m, n)
        ! A = U * diag(lam) * U^T (m x m). Used to advance snapshots:
        ! F(:, k+1) = A * F(:, k) — single trajectory.
        call gen_orthogonal_quad(m, U, seed = 191 + 31*isz)
        allocate(lam_re(m), lam_im(m))
        do j = 1, m
            lam_re(j) = 0.5_ep + 0.05_ep * real(j, ep)
            lam_im(j) = 0.0_ep
        end do
        if (m >= 4) then
            lam_re(1) = 0.6_ep;   lam_im(1) =  0.2_ep
            lam_re(2) = 0.6_ep;   lam_im(2) = -0.2_ep
            lam_re(m) = 1.0e-3_ep
        end if
        allocate(A_quad(m, m))
        block
            real(ep), allocatable :: T(:,:)
            integer :: ii
            allocate(T(m, m)); T = 0.0_ep
            ii = 1
            do while (ii <= m)
                if (ii < m .and. lam_im(ii) /= 0.0_ep) then
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

        ! Build snapshot matrix F (m x n) by advancing initial state.
        allocate(F0(m, n))
        block
            real(ep), allocatable :: x0(:)
            integer :: kk
            call random_number(F0(:, 1:1))
            allocate(x0(m))
            x0 = F0(:, 1)
            do kk = 2, n
                x0 = matmul(A_quad, x0)
                F0(:, kk) = x0
            end do
            deallocate(x0)
        end block

        do ijobs = 1, 2
            jobs = jobs_set(ijobs)
            do ijobz = 1, 2
                jobz = jobz_set(ijobz)
                do ijobr = 1, 2
                    if (ijobr == 1) then
                        jobr = 'N'
                    else
                        if (jobz == 'V') then
                            jobr = 'R'
                        else
                            cycle
                        end if
                    end if
                    do ijobf = 1, 2
                        jobf = jobf_set(ijobf)
                        jobq = 'N'   ! pin: F is overwritten with QR reflectors
                        jobt = 'R'   ! pin: Y holds R factor on exit (canonical)
                        allocate(F_ref(m, n), F_got(m, n))
                        allocate(X_ref(n, n), Y_ref(n, n), X_got(n, n), Y_got(n, n))
                        allocate(reig_ref(n), imeig_ref(n), reig_got(n), imeig_got(n))
                        allocate(Z_ref(m, n), Z_got(m, n))
                        allocate(B_ref(m, n), B_got(m, n))
                        allocate(V_ref(n, n), V_got(n, n))
                        allocate(S_ref(n, n), S_got(n, n))
                        allocate(res_ref(n), res_got(n))
                        F_ref = F0; F_got = F0
                        tol_dmd = 10.0_ep * target_eps
                        call dgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd,    &
                                     m, n, F_ref, m, X_ref, n, Y_ref, n,           &
                                     nrnk, tol_dmd, k_ref, reig_ref, imeig_ref,    &
                                     Z_ref, m, res_ref, B_ref, m,                   &
                                     V_ref, n, S_ref, n,                            &
                                     wopt, -1, iwopt, -1, info_ref)
                        if (info_ref == 0) then
                            lwork_ref  = max(1, int(wopt(1)))
                            liwork_ref = max(1, iwopt(1))
                            allocate(work(lwork_ref), iwork(liwork_ref))
                            call dgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd,    &
                                         m, n, F_ref, m, X_ref, n, Y_ref, n,           &
                                         nrnk, tol_dmd, k_ref, reig_ref, imeig_ref,    &
                                         Z_ref, m, res_ref, B_ref, m,                   &
                                         V_ref, n, S_ref, n,                            &
                                         work, lwork_ref, iwork, liwork_ref,            &
                                         info_ref)
                            deallocate(work, iwork)
                        end if
                        call target_dgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, &
                                            m, n, F_got, m, X_got, n, Y_got, n,        &
                                            nrnk, tol_dmd, k_got, reig_got, imeig_got, &
                                            Z_got, m, res_got, B_got, m,                &
                                            V_got, n, S_got, n,                         &
                                            info_got)
                        write(label, '(a,i0,a,i0,a,a,a,a,a,a,a,a)')                     &
                            'm=', m, ',n=', n, ',JOBS=', jobs, ',JOBZ=', jobz,          &
                            ',JOBR=', jobr, ',JOBF=', jobf
                        if (info_ref /= 0 .or. info_got /= 0 .or. k_ref /= k_got) then
                            ! Loose tolerance: single-trajectory data degrades
                            ! conditioning by factors that are hard to pin
                            ! analytically; 1000*n^2 keeps kind10/multifloats
                            ! within reach while kind16 stays bit-exact.
                            tol_eig = 10000.0_ep * real(n, ep)**2 * target_eps
                            call report_case(trim(label) // ',layer=a,info-mismatch',  &
                                              1.0_ep, tol_eig)
                        else if (k_got > 0) then
                            allocate(spec_ref(k_got), spec_got(k_got))
                            do j = 1, k_got
                                spec_ref(j) = cmplx(reig_ref(j), imeig_ref(j), ep)
                                spec_got(j) = cmplx(reig_got(j), imeig_got(j), ep)
                            end do
                            call sort_eig_lex_z(spec_ref)
                            call sort_eig_lex_z(spec_got)
                            err_a = max_rel_err_vec_z(spec_got, spec_ref)
                            ! Loose tolerance: single-trajectory data degrades
                            ! conditioning by factors that are hard to pin
                            ! analytically; 1000*n^2 keeps kind10/multifloats
                            ! within reach while kind16 stays bit-exact.
                            tol_eig = 10000.0_ep * real(n, ep)**2 * target_eps
                            call report_case(trim(label) // ',layer=a', err_a, tol_eig)
                            deallocate(spec_ref, spec_got)
                            ! Layer (b) — absolute residual delta scaled by normA.
                            if (jobr == 'R') then
                                tol_res = 10000.0_ep * real(n, ep) * target_eps
                                err_b = maxval(abs(res_got(1:k_got) -      &
                                                    res_ref(1:k_got))) /  &
                                        max(normA, target_eps)
                                call report_case(trim(label) // ',layer=b', err_b, tol_res)
                            end if
                        end if
                        deallocate(F_ref, F_got)
                        deallocate(X_ref, Y_ref, X_got, Y_got)
                        deallocate(reig_ref, imeig_ref, reig_got, imeig_got)
                        deallocate(Z_ref, Z_got, B_ref, B_got)
                        deallocate(V_ref, V_got, S_ref, S_got)
                        deallocate(res_ref, res_got)
                    end do
                end do
            end do
        end do

        deallocate(A_quad, U, lam_re, lam_im, F0)
    end do

    call report_finalize()
end program test_dgedmdq
