! zgedmd: complex Dynamic Mode Decomposition. Eigenvalues come back
! as a complex array EIGS(*) — no REIG/IMEIG assemble step. Same
! comparison strategy as test_dgedmd: layer (a) sorted spectrum, layer
! (b) residual-norm sanity (only when JOBR='R', which requires JOBZ='V').
program test_zgedmd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec, max_rel_err_vec_z, sort_eig_lex_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgedmd
    use ref_quad_lapack, only: zgedmd
    implicit none

    integer, parameter :: nsizes = 3
    integer, parameter :: ms(nsizes) = [8, 12, 24]
    integer, parameter :: ns(nsizes) = [5, 8, 16]
    character(len=1), parameter :: jobs_set(2) = ['N', 'S']
    character(len=1), parameter :: jobz_set(2) = ['N', 'V']
    character(len=1), parameter :: jobf_set(2) = ['N', 'E']

    integer :: isz, ijobs, ijobz, ijobr, ijobf
    integer :: m, n, info_ref, info_got, k_ref, k_got
    integer :: lzwork_ref, lrwork_ref, liwork_ref, j
    integer, parameter :: whtsvd = 1, nrnk = -1
    real(ep) :: tol_dmd, tol_eig, tol_res, normA, err_a, err_b
    complex(ep) :: zwopt(2)
    real(ep) :: rwopt(2)
    integer :: iwopt(2)
    complex(ep), allocatable :: A_quad(:,:), X0(:,:), Y0(:,:)
    complex(ep), allocatable :: X_ref(:,:), Y_ref(:,:), X_got(:,:), Y_got(:,:)
    complex(ep), allocatable :: eigs_ref(:), eigs_got(:)
    complex(ep), allocatable :: Z_ref(:,:), Z_got(:,:), B_ref(:,:), B_got(:,:)
    complex(ep), allocatable :: WW_ref(:,:), WW_got(:,:), S_ref(:,:), S_got(:,:)
    complex(ep), allocatable :: zwork(:)
    real(ep),    allocatable :: res_ref(:), res_got(:), rwork(:)
    integer,     allocatable :: iwork(:)
    complex(ep), allocatable :: spec_ref(:), spec_got(:)
    character(len=1) :: jobs, jobz, jobr, jobf
    character(len=80) :: label

    call report_init('zgedmd', target_name)

    do isz = 1, nsizes
        m = ms(isz); n = ns(isz)
        ! Random complex A; X random complex; Y = X * A^T (so the spectrum
        ! recovered from (X, Y) is sigma(A^T) = sigma(A)). normA used as
        ! the layer-(b) scale.
        call gen_matrix_complex(n, n, A_quad, seed = 33 + 41*isz)
        normA = sqrt(real(sum(abs(A_quad)**2), ep))

        call gen_matrix_complex(m, n, X0, seed = 711 + 19*isz)
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
                        if (jobz == 'V') then
                            jobr = 'R'
                        else
                            cycle
                        end if
                    end if
                    do ijobf = 1, 2
                        jobf = jobf_set(ijobf)
                        allocate(X_ref(m, n), Y_ref(m, n), X_got(m, n), Y_got(m, n))
                        allocate(eigs_ref(n), eigs_got(n))
                        allocate(Z_ref(m, n), Z_got(m, n))
                        allocate(B_ref(m, n), B_got(m, n))
                        allocate(WW_ref(n, n), WW_got(n, n))
                        allocate(S_ref(n, n), S_got(n, n))
                        allocate(res_ref(n), res_got(n))
                        X_ref = X0; Y_ref = Y0
                        X_got = X0; Y_got = Y0
                        tol_dmd = 10.0_ep * target_eps
                        call zgedmd(jobs, jobz, jobr, jobf, whtsvd, m, n,    &
                                    X_ref, m, Y_ref, m, nrnk, tol_dmd,        &
                                    k_ref, eigs_ref, Z_ref, m, res_ref,       &
                                    B_ref, m, WW_ref, n, S_ref, n,            &
                                    zwopt, -1, rwopt, -1, iwopt, -1, info_ref)
                        if (info_ref == 0) then
                            lzwork_ref = max(1, int(real(zwopt(1), ep)))
                            lrwork_ref = max(1, int(rwopt(1)))
                            liwork_ref = max(1, iwopt(1))
                            allocate(zwork(lzwork_ref), rwork(lrwork_ref), &
                                     iwork(liwork_ref))
                            call zgedmd(jobs, jobz, jobr, jobf, whtsvd, m, n,    &
                                        X_ref, m, Y_ref, m, nrnk, tol_dmd,       &
                                        k_ref, eigs_ref, Z_ref, m, res_ref,      &
                                        B_ref, m, WW_ref, n, S_ref, n,           &
                                        zwork, lzwork_ref, rwork, lrwork_ref,    &
                                        iwork, liwork_ref, info_ref)
                            deallocate(zwork, rwork, iwork)
                        end if
                        call target_zgedmd(jobs, jobz, jobr, jobf, whtsvd, m, n, &
                                           X_got, m, Y_got, m, nrnk, tol_dmd,    &
                                           k_got, eigs_got, Z_got, m, res_got,   &
                                           B_got, m, WW_got, n, S_got, n,        &
                                           info_got)
                        write(label, '(a,i0,a,i0,a,a,a,a,a,a,a,a)')              &
                            'm=', m, ',n=', n, ',JOBS=', jobs, ',JOBZ=', jobz,   &
                            ',JOBR=', jobr, ',JOBF=', jobf
                        if (info_ref /= 0 .or. info_got /= 0 .or. k_ref /= k_got) then
                            tol_eig = 100.0_ep * real(n, ep)**2 * target_eps
                            call report_case(trim(label) // ',layer=a,info-mismatch', &
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
                            tol_eig = 100.0_ep * real(n, ep)**2 * target_eps
                            call report_case(trim(label) // ',layer=a', err_a, tol_eig)
                            deallocate(spec_ref, spec_got)
                            if (jobr == 'R') then
                                tol_res = 100.0_ep * real(n, ep) * target_eps
                                err_b = maxval(abs(res_got(1:k_got) -      &
                                                    res_ref(1:k_got))) /  &
                                        max(normA, target_eps)
                                call report_case(trim(label) // ',layer=b', err_b, tol_res)
                            end if
                        end if
                        deallocate(X_ref, Y_ref, X_got, Y_got)
                        deallocate(eigs_ref, eigs_got)
                        deallocate(Z_ref, Z_got, B_ref, B_got)
                        deallocate(WW_ref, WW_got, S_ref, S_got)
                        deallocate(res_ref, res_got)
                    end do
                end do
            end do
        end do

        deallocate(A_quad, X0, Y0)
    end do

    call report_finalize()
end program test_zgedmd
