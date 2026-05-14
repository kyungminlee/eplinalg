! fuzz_egemm — consistency fuzz test for the kind10 egemm overlay.
!
! Each iteration:
!   1. Picks random M, N, K (log-uniform in [1, 256]).
!   2. Picks random LDA, LDB, LDC ≥ required, adding 0..4 extra rows
!      to stress stride handling.
!   3. Picks random TRANSA, TRANSB ∈ {'N','T','C'}.
!   4. Picks random ALPHA, BETA (mix of 0, ±1, uniform in [-1,1]).
!   5. Generates random matrices A, B, C0 in [-1, 1] at REAL(KIND=10).
!   6. Calls overlay egemm_ (resolved via composite eblas → C kernel)
!      and migrated egemm_migrated_ (resolved via the renamed archive,
!      Fortran kernel) on identical inputs.
!   7. Asserts max_rel_err ≤ max(10, 4*K) * eps10 — Wilkinson-style
!      bound, since reduction order differs between blocked-parallel
!      and unblocked-serial summation. Fails print the seed so the
!      case is reproducible.
!
! Knobs (env):
!   BLAS_FUZZ_SEED   — integer; default time-derived (printed at start)
!   BLAS_FUZZ_CASES  — number of cases; default 200

program fuzz_egemm
    use fuzz_util
    implicit none

    interface
        subroutine egemm(transa, transb, m, n, k, alpha, a, lda, &
                         b, ldb, beta, c, ldc)
            import :: rk10
            character, intent(in) :: transa, transb
            integer,   intent(in) :: m, n, k, lda, ldb, ldc
            real(rk10), intent(in) :: alpha, beta
            real(rk10), intent(in) :: a(lda, *), b(ldb, *)
            real(rk10), intent(inout) :: c(ldc, *)
        end subroutine egemm
        subroutine egemm_migrated(transa, transb, m, n, k, alpha, a, lda, &
                                  b, ldb, beta, c, ldc)
            import :: rk10
            character, intent(in) :: transa, transb
            integer,   intent(in) :: m, n, k, lda, ldb, ldc
            real(rk10), intent(in) :: alpha, beta
            real(rk10), intent(in) :: a(lda, *), b(ldb, *)
            real(rk10), intent(inout) :: c(ldc, *)
        end subroutine egemm_migrated
    end interface

    integer :: seed, ncases, i
    integer :: m, n, k, lda, ldb, ldc, rows_a, cols_a, rows_b, cols_b
    integer :: pad_a, pad_b, pad_c, fails
    character :: ta, tb
    real(rk10) :: alpha, beta, err, tol
    real(rk10), allocatable :: A(:), B(:), C0(:), C_ov(:), C_mig(:)
    real(8) :: u

    seed = read_seed()
    ncases = read_cases()
    print '(a,i0,a,i0)', 'fuzz_egemm: seed=', seed, '  cases=', ncases
    call seed_rng(seed)

    fails = 0
    do i = 1, ncases
        m = rand_int_log(1, 256)
        n = rand_int_log(1, 256)
        k = rand_int_log(1, 256)
        ! Tolerance scales with K: a length-K dot product (the inner
        ! reduction of GEMM) can legitimately differ between
        ! summation orders by up to ~2*K*eps (Wilkinson). Use a
        ! 4*K*eps envelope with a 10-ulp absolute floor for small K.
        tol = max(10.0_rk10, 4.0_rk10 * real(k, rk10)) * eps10
        ta = rand_trans()
        tb = rand_trans()
        call rand_alpha_beta(alpha, beta)

        if (ta == 'N') then
            rows_a = m; cols_a = k
        else
            rows_a = k; cols_a = m
        end if
        if (tb == 'N') then
            rows_b = k; cols_b = n
        else
            rows_b = n; cols_b = k
        end if
        call random_number(u); pad_a = int(u * 5.0_8)
        call random_number(u); pad_b = int(u * 5.0_8)
        call random_number(u); pad_c = int(u * 5.0_8)
        lda = rows_a + pad_a
        ldb = rows_b + pad_b
        ldc = m       + pad_c
        if (lda < 1) lda = 1
        if (ldb < 1) ldb = 1
        if (ldc < 1) ldc = 1

        allocate(A(lda * cols_a))
        allocate(B(ldb * cols_b))
        allocate(C0(ldc * n))
        allocate(C_ov(ldc * n))
        allocate(C_mig(ldc * n))
        call fill_matrix_10(A,  rows_a, cols_a, lda)
        call fill_matrix_10(B,  rows_b, cols_b, ldb)
        call fill_matrix_10(C0, m,      n,      ldc)
        C_ov  = C0
        C_mig = C0

        call egemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C_ov,  ldc)
        call egemm_migrated(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C_mig, ldc)
        err = max_rel_err_10(C_ov, C_mig, m, n, ldc)

        if (err > tol .or. err /= err) then  ! NaN-safe
            fails = fails + 1
            print '(a,i0,a,2a1,a,3i6,a,2es12.4,a,es12.4)', &
                'FAIL case ', i, '  trans=', ta, tb, &
                '  mnk=', m, n, k, &
                '  a,b=', alpha, beta, '  err=', err
        end if

        deallocate(A, B, C0, C_ov, C_mig)
    end do

    print '(a,i0,a,i0,a,i0)', &
        'fuzz_egemm: ', ncases - fails, '/', ncases, ' passed; fails=', fails
    if (fails > 0) error stop 1
end program fuzz_egemm
