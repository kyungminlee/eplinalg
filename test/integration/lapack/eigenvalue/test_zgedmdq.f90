! zgedmdq: complex snapshot-form Dynamic Mode Decomposition. Mirrors
! test_dgedmdq with EIGS(*) complex spectrum (no REIG/IMEIG split) and
! a real WORK + complex ZWORK + integer IWORK trio. Layer (c) skipped
! for the same T-matrix canonicality reason as test_dgedmdq.
program test_zgedmdq
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec, max_rel_err_vec_z, sort_eig_lex_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgedmdq
    use ref_quad_lapack, only: zgedmdq
    implicit none

    integer, parameter :: nsizes = 2
    integer, parameter :: ms(nsizes) = [8, 12]
    integer, parameter :: ns(nsizes) = [5, 8]
    character(len=1), parameter :: jobs_set(2) = ['N', 'S']
    character(len=1), parameter :: jobz_set(2) = ['N', 'V']
    character(len=1), parameter :: jobf_set(2) = ['N', 'E']

    integer :: isz, ijobs, ijobz, ijobr, ijobf
    integer :: m, n, info_ref, info_got, k_ref, k_got
    integer :: lzwork_ref, lwork_ref, liwork_ref, j, mn
    integer, parameter :: whtsvd = 1, nrnk = -1
    real(ep) :: tol_dmd, tol_eig, tol_res, normA, err_a, err_b
    complex(ep) :: zwopt(2)
    real(ep) :: wopt(2)
    integer :: iwopt(2)
    complex(ep), allocatable :: A_quad(:,:), F0(:,:)
    complex(ep), allocatable :: F_ref(:,:), F_got(:,:)
    complex(ep), allocatable :: X_ref(:,:), Y_ref(:,:), X_got(:,:), Y_got(:,:)
    complex(ep), allocatable :: eigs_ref(:), eigs_got(:)
    complex(ep), allocatable :: Z_ref(:,:), Z_got(:,:), B_ref(:,:), B_got(:,:)
    complex(ep), allocatable :: V_ref(:,:), V_got(:,:), S_ref(:,:), S_got(:,:)
    complex(ep), allocatable :: zwork(:)
    real(ep),    allocatable :: res_ref(:), res_got(:), work(:)
    integer,     allocatable :: iwork(:)
    complex(ep), allocatable :: spec_ref(:), spec_got(:)
    character(len=1) :: jobs, jobz, jobr, jobf, jobq, jobt
    character(len=80) :: label

    call report_init('zgedmdq', target_name)

    do isz = 1, nsizes
        m = ms(isz); n = ns(isz)
        mn = min(m, n)
        ! Random complex A (m-by-m), normA for layer-(b) scale.
        call gen_matrix_complex(m, m, A_quad, seed = 23 + 53*isz)
        normA = sqrt(real(sum(abs(A_quad)**2), ep))

        ! Snapshot trajectory: F(:,1) random, F(:,k+1) = A * F(:,k).
        allocate(F0(m, n))
        block
            real(ep), allocatable :: rR(:), rI(:)
            complex(ep), allocatable :: x0(:)
            integer :: kk
            allocate(rR(m), rI(m), x0(m))
            call random_number(rR); call random_number(rI)
            x0 = cmplx(rR, rI, ep)
            F0(:, 1) = x0
            do kk = 2, n
                x0 = matmul(A_quad, x0)
                F0(:, kk) = x0
            end do
            deallocate(rR, rI, x0)
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
                        jobq = 'N'
                        jobt = 'R'
                        allocate(F_ref(m, n), F_got(m, n))
                        allocate(X_ref(n, n), Y_ref(n, n), X_got(n, n), Y_got(n, n))
                        allocate(eigs_ref(n), eigs_got(n))
                        allocate(Z_ref(m, n), Z_got(m, n))
                        allocate(B_ref(m, n), B_got(m, n))
                        allocate(V_ref(n, n), V_got(n, n))
                        allocate(S_ref(n, n), S_got(n, n))
                        allocate(res_ref(n), res_got(n))
                        F_ref = F0; F_got = F0
                        tol_dmd = 10.0_ep * target_eps
                        call zgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd,    &
                                     m, n, F_ref, m, X_ref, n, Y_ref, n,           &
                                     nrnk, tol_dmd, k_ref, eigs_ref,                &
                                     Z_ref, m, res_ref, B_ref, m,                   &
                                     V_ref, n, S_ref, n,                            &
                                     zwopt, -1, wopt, -1, iwopt, -1, info_ref)
                        if (info_ref == 0) then
                            lzwork_ref = max(1, int(real(zwopt(1), ep)))
                            lwork_ref  = max(1, int(wopt(1)))
                            liwork_ref = max(1, iwopt(1))
                            allocate(zwork(lzwork_ref), work(lwork_ref), &
                                     iwork(liwork_ref))
                            call zgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, &
                                         m, n, F_ref, m, X_ref, n, Y_ref, n,         &
                                         nrnk, tol_dmd, k_ref, eigs_ref,             &
                                         Z_ref, m, res_ref, B_ref, m,                &
                                         V_ref, n, S_ref, n,                         &
                                         zwork, lzwork_ref, work, lwork_ref,         &
                                         iwork, liwork_ref, info_ref)
                            deallocate(zwork, work, iwork)
                        end if
                        call target_zgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, &
                                            m, n, F_got, m, X_got, n, Y_got, n,        &
                                            nrnk, tol_dmd, k_got, eigs_got,             &
                                            Z_got, m, res_got, B_got, m,                &
                                            V_got, n, S_got, n,                         &
                                            info_got)
                        write(label, '(a,i0,a,i0,a,a,a,a,a,a,a,a)')                     &
                            'm=', m, ',n=', n, ',JOBS=', jobs, ',JOBZ=', jobz,          &
                            ',JOBR=', jobr, ',JOBF=', jobf
                        if (info_ref /= 0 .or. info_got /= 0 .or. k_ref /= k_got) then
                            tol_eig = 10000.0_ep * real(n, ep)**2 * target_eps
                            call report_case(trim(label) // ',layer=a,info-mismatch',  &
                                              1.0_ep, tol_eig)
                        else if (k_got > 0) then
                            allocate(spec_ref(k_got), spec_got(k_got))
                            do j = 1, k_got
                                spec_ref(j) = eigs_ref(j)
                                spec_got(j) = eigs_got(j)
                            end do
                            call sort_eig_lex_z(spec_ref)
                            call sort_eig_lex_z(spec_got)
                            err_a = max_rel_err_vec_z(spec_got, spec_ref)
                            tol_eig = 10000.0_ep * real(n, ep)**2 * target_eps
                            call report_case(trim(label) // ',layer=a', err_a, tol_eig)
                            deallocate(spec_ref, spec_got)
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
                        deallocate(eigs_ref, eigs_got)
                        deallocate(Z_ref, Z_got, B_ref, B_got)
                        deallocate(V_ref, V_got, S_ref, S_got)
                        deallocate(res_ref, res_got)
                    end do
                end do
            end do
        end do

        deallocate(A_quad, F0)
    end do

    call report_finalize()
end program test_zgedmdq
