! fuzz_ygemm — consistency fuzz for the kind10 complex GEMM overlay.
!
! Same shape as fuzz_egemm; uses complex variants of the fuzz helpers.
! Both impls receive identical inputs; tolerance is the Wilkinson
! bound for length-K complex dot products (4*K*eps for the absolute
! diff between two orderings, with a 10-ulp floor for tiny K).

program fuzz_ygemm
    use fuzz_util
    implicit none

    interface
        subroutine ygemm(transa, transb, m, n, k, alpha, a, lda, &
                         b, ldb, beta, c, ldc)
            import :: rk10
            character, intent(in) :: transa, transb
            integer,   intent(in) :: m, n, k, lda, ldb, ldc
            complex(rk10), intent(in) :: alpha, beta
            complex(rk10), intent(in) :: a(lda, *), b(ldb, *)
            complex(rk10), intent(inout) :: c(ldc, *)
        end subroutine ygemm
        subroutine ygemm_migrated(transa, transb, m, n, k, alpha, a, lda, &
                                  b, ldb, beta, c, ldc)
            import :: rk10
            character, intent(in) :: transa, transb
            integer,   intent(in) :: m, n, k, lda, ldb, ldc
            complex(rk10), intent(in) :: alpha, beta
            complex(rk10), intent(in) :: a(lda, *), b(ldb, *)
            complex(rk10), intent(inout) :: c(ldc, *)
        end subroutine ygemm_migrated
    end interface

    integer :: seed, ncases, i
    integer :: m, n, k, lda, ldb, ldc, rows_a, cols_a, rows_b, cols_b
    integer :: fails, sentinel_fails, bad_ov, bad_mig
    character :: ta, tb
    complex(rk10) :: alpha, beta
    real(rk10) :: err, tol
    complex(rk10), allocatable :: A(:), B(:), C0(:), C_ov(:), C_mig(:)

    seed = read_seed()
    ncases = read_cases()
    print '(a,i0,a,i0)', 'fuzz_ygemm: seed=', seed, '  cases=', ncases
    call seed_rng(seed)

    fails = 0
    sentinel_fails = 0
    do i = 1, ncases
        m = rand_int_log(1, 256)
        n = rand_int_log(1, 256)
        k = rand_int_log(1, 256)
        ta = rand_trans_complex()
        tb = rand_trans_complex()
        call rand_alpha_beta_c10(alpha, beta)
        tol = max(10.0_rk10, 4.0_rk10 * real(k, rk10)) * eps10

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
        lda = rows_a + rand_pad()
        ldb = rows_b + rand_pad()
        ldc = m      + rand_pad()
        if (lda < 1) lda = 1
        if (ldb < 1) ldb = 1
        if (ldc < 1) ldc = 1

        allocate(A(lda * cols_a))
        allocate(B(ldb * cols_b))
        allocate(C0(ldc * n))
        allocate(C_ov(ldc * n))
        allocate(C_mig(ldc * n))
        call fill_matrix_c10(A,  rows_a, cols_a, lda)
        call fill_matrix_c10(B,  rows_b, cols_b, ldb)
        call fill_matrix_c10(C0, m,      n,      ldc)
        C_ov  = C0
        C_mig = C0

        call ygemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C_ov,  ldc)
        call ygemm_migrated(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C_mig, ldc)
        err = max_rel_err_c10(C_ov, C_mig, m, n, ldc)
        bad_ov  = check_sentinels_c10(C_ov,  m, n, ldc)
        bad_mig = check_sentinels_c10(C_mig, m, n, ldc)

        if (err > tol .or. err /= err) then
            fails = fails + 1
            print '(a,i0,a,2a1,a,3i6,a,3i6,a,es12.4,a,es12.4)', &
                'FAIL case ', i, '  trans=', ta, tb, &
                '  mnk=', m, n, k, &
                '  lda/ldb/ldc=', lda, ldb, ldc, &
                '  err=', err, '  tol=', tol
        end if
        if (bad_ov /= 0) then
            sentinel_fails = sentinel_fails + 1
            print '(a,i0,a,2a1,a,3i6,a,3i6,a,i0)', &
                'SENTINEL overlay case ', i, '  trans=', ta, tb, &
                '  mnk=', m, n, k, &
                '  lda/ldb/ldc=', lda, ldb, ldc, &
                '  bad col=', bad_ov
        end if
        if (bad_mig /= 0) then
            sentinel_fails = sentinel_fails + 1
            print '(a,i0,a,2a1,a,3i6,a,3i6,a,i0)', &
                'SENTINEL migrated case ', i, '  trans=', ta, tb, &
                '  mnk=', m, n, k, &
                '  lda/ldb/ldc=', lda, ldb, ldc, &
                '  bad col=', bad_mig
        end if

        deallocate(A, B, C0, C_ov, C_mig)
    end do

    print '(a,i0,a,i0,a,i0,a,i0)', &
        'fuzz_ygemm: ', ncases - fails, '/', ncases, &
        ' passed; tol_fails=', fails, ' sentinel_fails=', sentinel_fails
    if (fails > 0 .or. sentinel_fails > 0) error stop 1
end program fuzz_ygemm
