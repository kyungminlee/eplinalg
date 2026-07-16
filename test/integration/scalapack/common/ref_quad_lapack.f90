! Explicit interfaces to reflapack_quad — the vendored Netlib LAPACK
! compiled with gfortran's -freal-8-real-16 so REAL(KIND=8) and
! DOUBLE PRECISION entities are promoted to REAL(KIND=16) without
! routine renaming. The D/Z routines below operate at quad precision
! and serve as the universal reference for the differential precision
! tests.
!
! Calling these via implicit typing would be wrong: a caller without
! the explicit interface would default DGESV etc. to DOUBLE PRECISION
! (KIND=8), mismatching the promoted KIND=16 symbol. Always `use
! ref_quad_lapack`.

module ref_quad_lapack
    use prec_kinds, only: ep
    implicit none

    interface
        ! ── Linear solve / factorization — real ──────────────────────
        subroutine dgesv_quad(n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, lda, ldb
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine dgesv_quad

        subroutine dgetrf_quad(m, n, A, lda, ipiv, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda
            real(ep), intent(inout) :: A(lda,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine dgetrf_quad

        subroutine dgetrs_quad(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dgetrs_quad

        subroutine dpotrf_quad(uplo, n, A, lda, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine dpotrf_quad

        subroutine dpotrs_quad(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dpotrs_quad

        subroutine dgtsv_quad(n, nrhs, dl, d, du, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, ldb
            real(ep), intent(inout) :: dl(*), d(*), du(*), B(ldb,*)
            integer,  intent(out)   :: info
        end subroutine dgtsv_quad

        subroutine zgtsv_quad(n, nrhs, dl, d, du, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(inout) :: dl(*), d(*), du(*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zgtsv_quad

        subroutine dptsv_quad(n, nrhs, d, e, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, ldb
            real(ep), intent(inout) :: d(*), e(*), B(ldb,*)
            integer,  intent(out)   :: info
        end subroutine dptsv_quad

        subroutine zptsv_quad(n, nrhs, d, e, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, ldb
            real(ep),    intent(inout) :: d(*)
            complex(ep), intent(inout) :: e(*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zptsv_quad

        subroutine dpttrf_quad(n, d, e, info)
            import :: ep
            integer,  intent(in)    :: n
            real(ep), intent(inout) :: d(*), e(*)
            integer,  intent(out)   :: info
        end subroutine dpttrf_quad

        subroutine dpttrs_quad(n, nrhs, d, e, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, ldb
            real(ep), intent(in)    :: d(*), e(*)
            real(ep), intent(inout) :: B(ldb,*)
            integer,  intent(out)   :: info
        end subroutine dpttrs_quad

        subroutine zpttrf_quad(n, d, e, info)
            import :: ep
            integer,     intent(in)    :: n
            real(ep),    intent(inout) :: d(*)
            complex(ep), intent(inout) :: e(*)
            integer,     intent(out)   :: info
        end subroutine zpttrf_quad

        subroutine zpttrs_quad(uplo, n, nrhs, d, e, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            real(ep),    intent(in)    :: d(*)
            complex(ep), intent(in)    :: e(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpttrs_quad

        subroutine dgecon_quad(norm, n, A, lda, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)    :: norm
            integer,   intent(in)    :: n, lda
            real(ep),  intent(in)    :: A(lda,*), anorm
            real(ep),  intent(out)   :: rcond
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgecon_quad

        subroutine dgebal_quad(job, n, A, lda, ilo, ihi, scale, info)
            import :: ep
            character, intent(in)    :: job
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: ilo, ihi, info
            real(ep),  intent(out)   :: scale(*)
        end subroutine dgebal_quad

        subroutine dpoequ_quad(n, A, lda, S, scond, amax, info)
            import :: ep
            integer,  intent(in)    :: n, lda
            real(ep), intent(in)    :: A(lda,*)
            real(ep), intent(out)   :: S(*), scond, amax
            integer,  intent(out)   :: info
        end subroutine dpoequ_quad

        subroutine zgeequ_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda
            complex(ep), intent(in)    :: A(lda,*)
            real(ep),    intent(out)   :: R(*), C(*), rowcnd, colcnd, amax
            integer,     intent(out)   :: info
        end subroutine zgeequ_quad

        subroutine zpoequ_quad(n, A, lda, S, scond, amax, info)
            import :: ep
            integer,     intent(in)    :: n, lda
            complex(ep), intent(in)    :: A(lda,*)
            real(ep),    intent(out)   :: S(*), scond, amax
            integer,     intent(out)   :: info
        end subroutine zpoequ_quad

        subroutine dgeequ_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda
            real(ep), intent(in)    :: A(lda,*)
            real(ep), intent(out)   :: R(*), C(*), rowcnd, colcnd, amax
            integer,  intent(out)   :: info
        end subroutine dgeequ_quad

        subroutine dpocon_quad(uplo, n, A, lda, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(in)    :: A(lda,*), anorm
            real(ep),  intent(out)   :: rcond
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dpocon_quad

        subroutine dstebz_quad(range, order, n, vl, vu, il, iu, abstol, &
                          d, e, m, nsplit, w, iblock, isplit, &
                          work, iwork, info)
            import :: ep
            character, intent(in)    :: range, order
            integer,   intent(in)    :: n, il, iu
            real(ep),  intent(in)    :: vl, vu, abstol, d(*), e(*)
            integer,   intent(out)   :: m, nsplit
            real(ep),  intent(out)   :: w(*)
            integer,   intent(out)   :: iblock(*), isplit(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dstebz_quad

        ! ── Linear solve — complex ───────────────────────────────────
        subroutine zgesv_quad(n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgesv_quad

        ! ── QR factorization — real ──────────────────────────────────
        subroutine dgeqrf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqrf_quad

        subroutine dorgqr_quad(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, k, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorgqr_quad

        subroutine dgeqr2_quad(m, n, A, lda, tau, work, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqr2_quad

        subroutine dormqr_quad(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormqr_quad

        subroutine dgelqf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgelqf_quad

        subroutine dorglq_quad(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, k, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorglq_quad

        subroutine dormlq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormlq_quad

        subroutine dgeqlf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqlf_quad

        subroutine dorgql_quad(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, k, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorgql_quad

        subroutine dormql_quad(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormql_quad

        subroutine dgerqf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgerqf_quad

        subroutine dorgrq_quad(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, k, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorgrq_quad

        subroutine dormrq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormrq_quad

        subroutine dgels_quad(trans, m, n, nrhs, A, lda, B, ldb, &
                         work, lwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgels_quad

        subroutine dposv_quad(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dposv_quad

        subroutine dpotri_quad(uplo, n, A, lda, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine dpotri_quad

        subroutine dtrtri_quad(uplo, diag, n, A, lda, info)
            import :: ep
            character, intent(in)    :: uplo, diag
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine dtrtri_quad

        subroutine dgetri_quad(n, A, lda, ipiv, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: n, lda, lwork
            integer,  intent(in)    :: ipiv(*)
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dgetri_quad

        subroutine dgebrd_quad(m, n, A, lda, d, e, tauq, taup, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: d(*), e(*), tauq(*), taup(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgebrd_quad

        subroutine dormbr_quad(vect, side, trans, m, n, k, A, lda, tau, &
                          C, ldc, work, lwork, info)
            import :: ep
            character, intent(in)    :: vect, side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormbr_quad

        subroutine dtzrzf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dtzrzf_quad

        subroutine dormrz_quad(side, trans, m, n, k, l, A, lda, tau, &
                          C, ldc, work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, l, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormrz_quad

        subroutine dggqrf_quad(n, m, p, A, lda, taua, B, ldb, taub, &
                          work, lwork, info)
            import :: ep
            integer,  intent(in)    :: n, m, p, lda, ldb, lwork
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dggqrf_quad

        subroutine dggrqf_quad(m, p, n, A, lda, taua, B, ldb, taub, &
                          work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, p, n, lda, ldb, lwork
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dggrqf_quad

        subroutine dtrrfs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, &
                          X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo, trans, diag
            integer,   intent(in)  :: n, nrhs, lda, ldb, ldx
            real(ep),  intent(in)  :: A(lda,*), B(ldb,*), X(ldx,*)
            real(ep),  intent(out) :: ferr(*), berr(*), work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtrrfs_quad

        subroutine dtrcon_quad(norm, uplo, diag, n, A, lda, rcond, &
                          work, iwork, info)
            import :: ep
            character, intent(in)  :: norm, uplo, diag
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtrcon_quad

        subroutine ztrcon_quad(norm, uplo, diag, n, A, lda, rcond, &
                          work, rwork, info)
            import :: ep
            character,   intent(in)  :: norm, uplo, diag
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            complex(ep), intent(out) :: work(*)
            real(ep),    intent(out) :: rcond, rwork(*)
            integer,     intent(out) :: info
        end subroutine ztrcon_quad

        subroutine ztrevc_quad(side, howmny, select, n, T, ldt, &
                          VL, ldvl, VR, ldvr, mm, m, work, rwork, info)
            import :: ep
            character,   intent(in)    :: side, howmny
            logical,     intent(in)    :: select(*)
            integer,     intent(in)    :: n, ldt, ldvl, ldvr, mm
            complex(ep), intent(in)    :: T(ldt,*)
            complex(ep), intent(inout) :: VL(ldvl,*), VR(ldvr,*)
            integer,     intent(out)   :: m
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine ztrevc_quad

        subroutine dtrsen_quad(job, compq, select, n, T, ldt, Q, ldq, &
                          WR, WI, m, S, SEP, work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: job, compq
            logical,   intent(in)    :: select(*)
            integer,   intent(in)    :: n, ldt, ldq, lwork, liwork
            real(ep),  intent(inout) :: T(ldt,*), Q(ldq,*)
            real(ep),  intent(out)   :: WR(*), WI(*), S, SEP, work(*)
            integer,   intent(out)   :: m, iwork(*), info
        end subroutine dtrsen_quad

        subroutine dhseqr_quad(job, compz, n, ilo, ihi, H, ldh, WR, WI, &
                          Z, ldz, work, lwork, info)
            import :: ep
            character, intent(in)    :: job, compz
            integer,   intent(in)    :: n, ilo, ihi, ldh, ldz, lwork
            real(ep),  intent(inout) :: H(ldh,*), Z(ldz,*)
            real(ep),  intent(out)   :: WR(*), WI(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dhseqr_quad

        subroutine ztrrfs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, &
                          X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo, trans, diag
            integer,     intent(in)  :: n, nrhs, lda, ldb, ldx
            complex(ep), intent(in)  :: A(lda,*), B(ldb,*), X(ldx,*)
            complex(ep), intent(out) :: work(*)
            real(ep),    intent(out) :: ferr(*), berr(*), rwork(*)
            integer,     intent(out) :: info
        end subroutine ztrrfs_quad

        subroutine dgehrd_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: n, ilo, ihi, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgehrd_quad

        subroutine dsytrd_quad(uplo, n, A, lda, d, e, tau, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: d(*), e(*), tau(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsytrd_quad

        subroutine dsygst_quad(itype, uplo, n, A, lda, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: itype, n, lda, ldb
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(in)    :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dsygst_quad

        subroutine dormtr_quad(side, uplo, trans, m, n, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: side, uplo, trans
            integer,   intent(in)    :: m, n, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormtr_quad

        subroutine dormhr_quad(side, trans, m, n, ilo, ihi, A, lda, tau, &
                          C, ldc, work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, ilo, ihi, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormhr_quad

        ! ── QR factorization — complex ───────────────────────────────
        subroutine zgeqrf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqrf_quad

        subroutine zgeqr2_quad(m, n, A, lda, tau, work, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqr2_quad

        subroutine zungqr_quad(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungqr_quad

        subroutine zunmqr_quad(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmqr_quad

        subroutine zgelqf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgelqf_quad

        subroutine zunglq_quad(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunglq_quad

        subroutine zunmlq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmlq_quad

        subroutine zgeqlf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqlf_quad

        subroutine zungql_quad(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungql_quad

        subroutine zunmql_quad(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmql_quad

        subroutine zgerqf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgerqf_quad

        subroutine zungrq_quad(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungrq_quad

        subroutine zunmrq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmrq_quad

        subroutine zgels_quad(trans, m, n, nrhs, A, lda, B, ldb, &
                         work, lwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgels_quad

        subroutine zposv_quad(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zposv_quad

        subroutine zhegst_quad(itype, uplo, n, A, lda, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: itype, n, lda, ldb
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zhegst_quad

        subroutine zpotrf_quad(uplo, n, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine zpotrf_quad

        subroutine zpocon_quad(uplo, n, A, lda, anorm, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(in)    :: A(lda,*)
            real(ep),    intent(in)    :: anorm
            real(ep),    intent(out)   :: rcond
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zpocon_quad

        subroutine dtrtrs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dtrtrs_quad

        subroutine ztrtrs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo, trans, diag
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine ztrtrs_quad

        subroutine zpotrs_quad(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpotrs_quad

        subroutine zgetrs_quad(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, nrhs, lda, ldb, ipiv(*)
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zgetrs_quad

        subroutine zgetrf_quad(m, n, A, lda, ipiv, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgetrf_quad

        subroutine zpotri_quad(uplo, n, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine zpotri_quad

        subroutine ztrtri_quad(uplo, diag, n, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo, diag
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine ztrtri_quad

        subroutine zgetri_quad(n, A, lda, ipiv, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, lda, lwork
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgetri_quad

        subroutine zgebrd_quad(m, n, A, lda, d, e, tauq, taup, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: d(*), e(*)
            complex(ep), intent(out)   :: tauq(*), taup(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgebrd_quad

        subroutine zgehrd_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, ilo, ihi, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgehrd_quad

        subroutine zhetrd_quad(uplo, n, A, lda, d, e, tau, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: d(*), e(*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhetrd_quad

        subroutine zunmbr_quad(vect, side, trans, m, n, k, A, lda, tau, &
                          C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: vect, side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmbr_quad

        subroutine ztzrzf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine ztzrzf_quad

        subroutine zunmrz_quad(side, trans, m, n, k, l, A, lda, tau, &
                          C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, l, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmrz_quad

        subroutine zggqrf_quad(n, m, p, A, lda, taua, B, ldb, taub, &
                          work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, m, p, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zggqrf_quad

        subroutine zggrqf_quad(m, p, n, A, lda, taua, B, ldb, taub, &
                          work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, p, n, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zggrqf_quad

        subroutine zunmtr_quad(side, uplo, trans, m, n, A, lda, tau, &
                          C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, uplo, trans
            integer,     intent(in)    :: m, n, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmtr_quad

        subroutine zunmhr_quad(side, trans, m, n, ilo, ihi, A, lda, tau, &
                          C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, ilo, ihi, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmhr_quad

        ! ── Symmetric / Hermitian eigenvalue ─────────────────────────
        subroutine dsyev_quad(jobz, uplo, n, A, lda, w, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsyev_quad

        subroutine dsyevd_quad(jobz, uplo, n, A, lda, w, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork, liwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsyevd_quad

        subroutine zheev_quad(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zheev_quad

        subroutine zheevd_quad(jobz, uplo, n, A, lda, w, work, lwork, &
                          rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, lda, lwork, lrwork, liwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zheevd_quad

        ! ── SVD ──────────────────────────────────────────────────────
        subroutine dgesvd_quad(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: jobu, jobvt
            integer,   intent(in)    :: m, n, lda, ldu, ldvt, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dgesvd_quad

        subroutine zgesvd_quad(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, &
                          work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobu, jobvt
            integer,     intent(in)    :: m, n, lda, ldu, ldvt, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: s(*), rwork(*)
            complex(ep), intent(out)   :: U(ldu,*), VT(ldvt,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgesvd_quad

        ! ── Auxiliary ────────────────────────────────────────────────
        function dlange_quad(norm, m, n, A, lda, work) result(r)
            import :: ep
            character, intent(in) :: norm
            integer,   intent(in) :: m, n, lda
            real(ep),  intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function dlange_quad

        function dlansy_quad(norm, uplo, n, A, lda, work) result(r)
            import :: ep
            character, intent(in) :: norm, uplo
            integer,   intent(in) :: n, lda
            real(ep),  intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function dlansy_quad

        function dlanhs_quad(norm, n, A, lda, work) result(r)
            import :: ep
            character, intent(in) :: norm
            integer,   intent(in) :: n, lda
            real(ep),  intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function dlanhs_quad

        function dlantr_quad(norm, uplo, diag, m, n, A, lda, work) result(r)
            import :: ep
            character, intent(in) :: norm, uplo, diag
            integer,   intent(in) :: m, n, lda
            real(ep),  intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function dlantr_quad

        function zlange_quad(norm, m, n, A, lda, work) result(r)
            import :: ep
            character,   intent(in) :: norm
            integer,     intent(in) :: m, n, lda
            complex(ep), intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function zlange_quad

        function zlanhe_quad(norm, uplo, n, A, lda, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo
            integer,     intent(in) :: n, lda
            complex(ep), intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function zlanhe_quad

        function zlansy_quad(norm, uplo, n, A, lda, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo
            integer,     intent(in) :: n, lda
            complex(ep), intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function zlansy_quad

        function zlanhs_quad(norm, n, A, lda, work) result(r)
            import :: ep
            character,   intent(in) :: norm
            integer,     intent(in) :: n, lda
            complex(ep), intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function zlanhs_quad

        function zlantr_quad(norm, uplo, diag, m, n, A, lda, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo, diag
            integer,     intent(in) :: m, n, lda
            complex(ep), intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function zlantr_quad

        subroutine dlacpy_quad(uplo, m, n, A, lda, B, ldb)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: m, n, lda, ldb
            real(ep),  intent(in)  :: A(lda,*)
            real(ep),  intent(out) :: B(ldb,*)
        end subroutine dlacpy_quad

        subroutine dlaset_quad(uplo, m, n, alpha, beta, A, lda)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: m, n, lda
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(inout) :: A(lda,*)
        end subroutine dlaset_quad
    end interface

contains

    subroutine dgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
        integer,  intent(in)    :: n, nrhs, lda, ldb
        real(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,  intent(out)   :: ipiv(*), info
        call dgesv_quad(n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine dgesv

    subroutine dgetrf(m, n, A, lda, ipiv, info)
        integer,  intent(in)    :: m, n, lda
        real(ep), intent(inout) :: A(lda,*)
        integer,  intent(out)   :: ipiv(*), info
        call dgetrf_quad(m, n, A, lda, ipiv, info)
    end subroutine dgetrf

    subroutine dgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dgetrs_quad(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine dgetrs

    subroutine dpotrf(uplo, n, A, lda, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: info
        call dpotrf_quad(uplo, n, A, lda, info)
    end subroutine dpotrf

    subroutine dpotrs(uplo, n, nrhs, A, lda, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dpotrs_quad(uplo, n, nrhs, A, lda, B, ldb, info)
    end subroutine dpotrs

    subroutine dgtsv(n, nrhs, dl, d, du, B, ldb, info)
        integer,  intent(in)    :: n, nrhs, ldb
        real(ep), intent(inout) :: dl(*), d(*), du(*), B(ldb,*)
        integer,  intent(out)   :: info
        call dgtsv_quad(n, nrhs, dl, d, du, B, ldb, info)
    end subroutine dgtsv

    subroutine zgtsv(n, nrhs, dl, d, du, B, ldb, info)
        integer,     intent(in)    :: n, nrhs, ldb
        complex(ep), intent(inout) :: dl(*), d(*), du(*), B(ldb,*)
        integer,     intent(out)   :: info
        call zgtsv_quad(n, nrhs, dl, d, du, B, ldb, info)
    end subroutine zgtsv

    subroutine dptsv(n, nrhs, d, e, B, ldb, info)
        integer,  intent(in)    :: n, nrhs, ldb
        real(ep), intent(inout) :: d(*), e(*), B(ldb,*)
        integer,  intent(out)   :: info
        call dptsv_quad(n, nrhs, d, e, B, ldb, info)
    end subroutine dptsv

    subroutine zptsv(n, nrhs, d, e, B, ldb, info)
        integer,     intent(in)    :: n, nrhs, ldb
        real(ep),    intent(inout) :: d(*)
        complex(ep), intent(inout) :: e(*), B(ldb,*)
        integer,     intent(out)   :: info
        call zptsv_quad(n, nrhs, d, e, B, ldb, info)
    end subroutine zptsv

    subroutine dpttrf(n, d, e, info)
        integer,  intent(in)    :: n
        real(ep), intent(inout) :: d(*), e(*)
        integer,  intent(out)   :: info
        call dpttrf_quad(n, d, e, info)
    end subroutine dpttrf

    subroutine dpttrs(n, nrhs, d, e, B, ldb, info)
        integer,  intent(in)    :: n, nrhs, ldb
        real(ep), intent(in)    :: d(*), e(*)
        real(ep), intent(inout) :: B(ldb,*)
        integer,  intent(out)   :: info
        call dpttrs_quad(n, nrhs, d, e, B, ldb, info)
    end subroutine dpttrs

    subroutine zpttrf(n, d, e, info)
        integer,     intent(in)    :: n
        real(ep),    intent(inout) :: d(*)
        complex(ep), intent(inout) :: e(*)
        integer,     intent(out)   :: info
        call zpttrf_quad(n, d, e, info)
    end subroutine zpttrf

    subroutine zpttrs(uplo, n, nrhs, d, e, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, ldb
        real(ep),    intent(in)    :: d(*)
        complex(ep), intent(in)    :: e(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zpttrs_quad(uplo, n, nrhs, d, e, B, ldb, info)
    end subroutine zpttrs

    subroutine dgecon(norm, n, A, lda, anorm, rcond, work, iwork, info)
        character, intent(in)    :: norm
        integer,   intent(in)    :: n, lda
        real(ep),  intent(in)    :: A(lda,*), anorm
        real(ep),  intent(out)   :: rcond
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgecon_quad(norm, n, A, lda, anorm, rcond, work, iwork, info)
    end subroutine dgecon

    subroutine dgebal(job, n, A, lda, ilo, ihi, scale, info)
        character, intent(in)    :: job
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: ilo, ihi, info
        real(ep),  intent(out)   :: scale(*)
        call dgebal_quad(job, n, A, lda, ilo, ihi, scale, info)
    end subroutine dgebal

    subroutine dpoequ(n, A, lda, S, scond, amax, info)
        integer,  intent(in)    :: n, lda
        real(ep), intent(in)    :: A(lda,*)
        real(ep), intent(out)   :: S(*), scond, amax
        integer,  intent(out)   :: info
        call dpoequ_quad(n, A, lda, S, scond, amax, info)
    end subroutine dpoequ

    subroutine zgeequ(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
        integer,     intent(in)    :: m, n, lda
        complex(ep), intent(in)    :: A(lda,*)
        real(ep),    intent(out)   :: R(*), C(*), rowcnd, colcnd, amax
        integer,     intent(out)   :: info
        call zgeequ_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
    end subroutine zgeequ

    subroutine zpoequ(n, A, lda, S, scond, amax, info)
        integer,     intent(in)    :: n, lda
        complex(ep), intent(in)    :: A(lda,*)
        real(ep),    intent(out)   :: S(*), scond, amax
        integer,     intent(out)   :: info
        call zpoequ_quad(n, A, lda, S, scond, amax, info)
    end subroutine zpoequ

    subroutine dgeequ(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
        integer,  intent(in)    :: m, n, lda
        real(ep), intent(in)    :: A(lda,*)
        real(ep), intent(out)   :: R(*), C(*), rowcnd, colcnd, amax
        integer,  intent(out)   :: info
        call dgeequ_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
    end subroutine dgeequ

    subroutine dpocon(uplo, n, A, lda, anorm, rcond, work, iwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(in)    :: A(lda,*), anorm
        real(ep),  intent(out)   :: rcond
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: iwork(*), info
        call dpocon_quad(uplo, n, A, lda, anorm, rcond, work, iwork, info)
    end subroutine dpocon

    subroutine dstebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, &
                      work, iwork, info)
        character, intent(in)    :: range, order
        integer,   intent(in)    :: n, il, iu
        real(ep),  intent(in)    :: vl, vu, abstol, d(*), e(*)
        integer,   intent(out)   :: m, nsplit
        real(ep),  intent(out)   :: w(*)
        integer,   intent(out)   :: iblock(*), isplit(*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: iwork(*), info
        call dstebz_quad(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit, &
                         work, iwork, info)
    end subroutine dstebz

    subroutine zgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        call zgesv_quad(n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine zgesv

    subroutine dgeqrf(m, n, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dgeqrf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine dgeqrf

    subroutine dorgqr(m, n, k, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, k, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: tau(*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dorgqr_quad(m, n, k, A, lda, tau, work, lwork, info)
    end subroutine dorgqr

    subroutine dgeqr2(m, n, A, lda, tau, work, info)
        integer,  intent(in)    :: m, n, lda
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dgeqr2_quad(m, n, A, lda, tau, work, info)
    end subroutine dgeqr2

    subroutine dormqr(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormqr_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormqr

    subroutine dgelqf(m, n, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dgelqf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine dgelqf

    subroutine dorglq(m, n, k, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, k, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: tau(*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dorglq_quad(m, n, k, A, lda, tau, work, lwork, info)
    end subroutine dorglq

    subroutine dormlq(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormlq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormlq

    subroutine dgeqlf(m, n, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dgeqlf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine dgeqlf

    subroutine dorgql(m, n, k, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, k, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: tau(*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dorgql_quad(m, n, k, A, lda, tau, work, lwork, info)
    end subroutine dorgql

    subroutine dormql(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormql_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormql

    subroutine dgerqf(m, n, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dgerqf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine dgerqf

    subroutine dorgrq(m, n, k, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, k, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: tau(*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dorgrq_quad(m, n, k, A, lda, tau, work, lwork, info)
    end subroutine dorgrq

    subroutine dormrq(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormrq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormrq

    subroutine dgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, nrhs, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dgels_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
    end subroutine dgels

    subroutine dposv(uplo, n, nrhs, A, lda, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        integer,   intent(out)   :: info
        call dposv_quad(uplo, n, nrhs, A, lda, B, ldb, info)
    end subroutine dposv

    subroutine dpotri(uplo, n, A, lda, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: info
        call dpotri_quad(uplo, n, A, lda, info)
    end subroutine dpotri

    subroutine dtrtri(uplo, diag, n, A, lda, info)
        character, intent(in)    :: uplo, diag
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: info
        call dtrtri_quad(uplo, diag, n, A, lda, info)
    end subroutine dtrtri

    subroutine dgetri(n, A, lda, ipiv, work, lwork, info)
        integer,  intent(in)    :: n, lda, lwork
        integer,  intent(in)    :: ipiv(*)
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dgetri_quad(n, A, lda, ipiv, work, lwork, info)
    end subroutine dgetri

    subroutine dgebrd(m, n, A, lda, d, e, tauq, taup, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: d(*), e(*), tauq(*), taup(*), work(*)
        integer,  intent(out)   :: info
        call dgebrd_quad(m, n, A, lda, d, e, tauq, taup, work, lwork, info)
    end subroutine dgebrd

    subroutine dormbr(vect, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: vect, side, trans
        integer,   intent(in)    :: m, n, k, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormbr_quad(vect, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormbr

    subroutine dtzrzf(m, n, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dtzrzf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine dtzrzf

    subroutine dormrz(side, trans, m, n, k, l, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, l, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormrz_quad(side, trans, m, n, k, l, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormrz

    subroutine dggqrf(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
        integer,  intent(in)    :: n, m, p, lda, ldb, lwork
        real(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep), intent(out)   :: taua(*), taub(*), work(*)
        integer,  intent(out)   :: info
        call dggqrf_quad(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
    end subroutine dggqrf

    subroutine dggrqf(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
        integer,  intent(in)    :: m, p, n, lda, ldb, lwork
        real(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep), intent(out)   :: taua(*), taub(*), work(*)
        integer,  intent(out)   :: info
        call dggrqf_quad(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
    end subroutine dggrqf

    subroutine dtrrfs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, X, ldx, ferr, berr, work, iwork, &
                      info)
        character, intent(in)  :: uplo, trans, diag
        integer,   intent(in)  :: n, nrhs, lda, ldb, ldx
        real(ep),  intent(in)  :: A(lda,*), B(ldb,*), X(ldx,*)
        real(ep),  intent(out) :: ferr(*), berr(*), work(*)
        integer,   intent(out) :: iwork(*), info
        call dtrrfs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, X, ldx, ferr, berr, work, iwork, &
                         info)
    end subroutine dtrrfs

    subroutine dtrcon(norm, uplo, diag, n, A, lda, rcond, work, iwork, info)
        character, intent(in)  :: norm, uplo, diag
        integer,   intent(in)  :: n, lda
        real(ep),  intent(in)  :: A(lda,*)
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dtrcon_quad(norm, uplo, diag, n, A, lda, rcond, work, iwork, info)
    end subroutine dtrcon

    subroutine ztrcon(norm, uplo, diag, n, A, lda, rcond, work, rwork, info)
        character,   intent(in)  :: norm, uplo, diag
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        complex(ep), intent(out) :: work(*)
        real(ep),    intent(out) :: rcond, rwork(*)
        integer,     intent(out) :: info
        call ztrcon_quad(norm, uplo, diag, n, A, lda, rcond, work, rwork, info)
    end subroutine ztrcon

    subroutine ztrevc(side, howmny, select, n, T, ldt, VL, ldvl, VR, ldvr, mm, m, work, rwork, info)
        character,   intent(in)    :: side, howmny
        logical,     intent(in)    :: select(*)
        integer,     intent(in)    :: n, ldt, ldvl, ldvr, mm
        complex(ep), intent(in)    :: T(ldt,*)
        complex(ep), intent(inout) :: VL(ldvl,*), VR(ldvr,*)
        integer,     intent(out)   :: m
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: info
        call ztrevc_quad(side, howmny, select, n, T, ldt, VL, ldvl, VR, ldvr, mm, m, work, rwork, &
                         info)
    end subroutine ztrevc

    subroutine dtrsen(job, compq, select, n, T, ldt, Q, ldq, WR, WI, m, S, SEP, work, lwork, iwork, &
                      liwork, info)
        character, intent(in)    :: job, compq
        logical,   intent(in)    :: select(*)
        integer,   intent(in)    :: n, ldt, ldq, lwork, liwork
        real(ep),  intent(inout) :: T(ldt,*), Q(ldq,*)
        real(ep),  intent(out)   :: WR(*), WI(*), S, SEP, work(*)
        integer,   intent(out)   :: m, iwork(*), info
        call dtrsen_quad(job, compq, select, n, T, ldt, Q, ldq, WR, WI, m, S, SEP, work, lwork, &
                         iwork, liwork, info)
    end subroutine dtrsen

    subroutine dhseqr(job, compz, n, ilo, ihi, H, ldh, WR, WI, Z, ldz, work, lwork, info)
        character, intent(in)    :: job, compz
        integer,   intent(in)    :: n, ilo, ihi, ldh, ldz, lwork
        real(ep),  intent(inout) :: H(ldh,*), Z(ldz,*)
        real(ep),  intent(out)   :: WR(*), WI(*), work(*)
        integer,   intent(out)   :: info
        call dhseqr_quad(job, compz, n, ilo, ihi, H, ldh, WR, WI, Z, ldz, work, lwork, info)
    end subroutine dhseqr

    subroutine ztrrfs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, X, ldx, ferr, berr, work, rwork, &
                      info)
        character,   intent(in)  :: uplo, trans, diag
        integer,     intent(in)  :: n, nrhs, lda, ldb, ldx
        complex(ep), intent(in)  :: A(lda,*), B(ldb,*), X(ldx,*)
        complex(ep), intent(out) :: work(*)
        real(ep),    intent(out) :: ferr(*), berr(*), rwork(*)
        integer,     intent(out) :: info
        call ztrrfs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, X, ldx, ferr, berr, work, rwork, &
                         info)
    end subroutine ztrrfs

    subroutine dgehrd(n, ilo, ihi, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: n, ilo, ihi, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dgehrd_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
    end subroutine dgehrd

    subroutine dsytrd(uplo, n, A, lda, d, e, tau, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: d(*), e(*), tau(*), work(*)
        integer,   intent(out)   :: info
        call dsytrd_quad(uplo, n, A, lda, d, e, tau, work, lwork, info)
    end subroutine dsytrd

    subroutine dsygst(itype, uplo, n, A, lda, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: itype, n, lda, ldb
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(in)    :: B(ldb,*)
        integer,   intent(out)   :: info
        call dsygst_quad(itype, uplo, n, A, lda, B, ldb, info)
    end subroutine dsygst

    subroutine dormtr(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, uplo, trans
        integer,   intent(in)    :: m, n, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormtr_quad(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormtr

    subroutine dormhr(side, trans, m, n, ilo, ihi, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, ilo, ihi, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormhr_quad(side, trans, m, n, ilo, ihi, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormhr

    subroutine zgeqrf(m, n, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zgeqrf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine zgeqrf

    subroutine zgeqr2(m, n, A, lda, tau, work, info)
        integer,     intent(in)    :: m, n, lda
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zgeqr2_quad(m, n, A, lda, tau, work, info)
    end subroutine zgeqr2

    subroutine zungqr(m, n, k, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, k, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: tau(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zungqr_quad(m, n, k, A, lda, tau, work, lwork, info)
    end subroutine zungqr

    subroutine zunmqr(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmqr_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmqr

    subroutine zgelqf(m, n, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zgelqf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine zgelqf

    subroutine zunglq(m, n, k, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, k, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: tau(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunglq_quad(m, n, k, A, lda, tau, work, lwork, info)
    end subroutine zunglq

    subroutine zunmlq(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmlq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmlq

    subroutine zgeqlf(m, n, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zgeqlf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine zgeqlf

    subroutine zungql(m, n, k, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, k, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: tau(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zungql_quad(m, n, k, A, lda, tau, work, lwork, info)
    end subroutine zungql

    subroutine zunmql(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmql_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmql

    subroutine zgerqf(m, n, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zgerqf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine zgerqf

    subroutine zungrq(m, n, k, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, k, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: tau(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zungrq_quad(m, n, k, A, lda, tau, work, lwork, info)
    end subroutine zungrq

    subroutine zunmrq(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmrq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmrq

    subroutine zgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgels_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
    end subroutine zgels

    subroutine zposv(uplo, n, nrhs, A, lda, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: info
        call zposv_quad(uplo, n, nrhs, A, lda, B, ldb, info)
    end subroutine zposv

    subroutine zhegst(itype, uplo, n, A, lda, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: itype, n, lda, ldb
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: B(ldb,*)
        integer,     intent(out)   :: info
        call zhegst_quad(itype, uplo, n, A, lda, B, ldb, info)
    end subroutine zhegst

    subroutine zpotrf(uplo, n, A, lda, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: info
        call zpotrf_quad(uplo, n, A, lda, info)
    end subroutine zpotrf

    subroutine zpocon(uplo, n, A, lda, anorm, rcond, work, rwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(in)    :: A(lda,*)
        real(ep),    intent(in)    :: anorm
        real(ep),    intent(out)   :: rcond
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: info
        call zpocon_quad(uplo, n, A, lda, anorm, rcond, work, rwork, info)
    end subroutine zpocon

    subroutine dtrtrs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dtrtrs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
    end subroutine dtrtrs

    subroutine ztrtrs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
        character,   intent(in)    :: uplo, trans, diag
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call ztrtrs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
    end subroutine ztrtrs

    subroutine zpotrs(uplo, n, nrhs, A, lda, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zpotrs_quad(uplo, n, nrhs, A, lda, B, ldb, info)
    end subroutine zpotrs

    subroutine zgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: n, nrhs, lda, ldb, ipiv(*)
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zgetrs_quad(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine zgetrs

    subroutine zgetrf(m, n, A, lda, ipiv, info)
        integer,     intent(in)    :: m, n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ipiv(*), info
        call zgetrf_quad(m, n, A, lda, ipiv, info)
    end subroutine zgetrf

    subroutine zpotri(uplo, n, A, lda, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: info
        call zpotri_quad(uplo, n, A, lda, info)
    end subroutine zpotri

    subroutine ztrtri(uplo, diag, n, A, lda, info)
        character,   intent(in)    :: uplo, diag
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: info
        call ztrtri_quad(uplo, diag, n, A, lda, info)
    end subroutine ztrtri

    subroutine zgetri(n, A, lda, ipiv, work, lwork, info)
        integer,     intent(in)    :: n, lda, lwork
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgetri_quad(n, A, lda, ipiv, work, lwork, info)
    end subroutine zgetri

    subroutine zgebrd(m, n, A, lda, d, e, tauq, taup, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: d(*), e(*)
        complex(ep), intent(out)   :: tauq(*), taup(*), work(*)
        integer,     intent(out)   :: info
        call zgebrd_quad(m, n, A, lda, d, e, tauq, taup, work, lwork, info)
    end subroutine zgebrd

    subroutine zgehrd(n, ilo, ihi, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: n, ilo, ihi, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zgehrd_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
    end subroutine zgehrd

    subroutine zhetrd(uplo, n, A, lda, d, e, tau, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: d(*), e(*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zhetrd_quad(uplo, n, A, lda, d, e, tau, work, lwork, info)
    end subroutine zhetrd

    subroutine zunmbr(vect, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: vect, side, trans
        integer,     intent(in)    :: m, n, k, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmbr_quad(vect, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmbr

    subroutine ztzrzf(m, n, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call ztzrzf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine ztzrzf

    subroutine zunmrz(side, trans, m, n, k, l, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, l, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmrz_quad(side, trans, m, n, k, l, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmrz

    subroutine zggqrf(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
        integer,     intent(in)    :: n, m, p, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: taua(*), taub(*), work(*)
        integer,     intent(out)   :: info
        call zggqrf_quad(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
    end subroutine zggqrf

    subroutine zggrqf(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
        integer,     intent(in)    :: m, p, n, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: taua(*), taub(*), work(*)
        integer,     intent(out)   :: info
        call zggrqf_quad(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
    end subroutine zggrqf

    subroutine zunmtr(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, uplo, trans
        integer,     intent(in)    :: m, n, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmtr_quad(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmtr

    subroutine zunmhr(side, trans, m, n, ilo, ihi, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, ilo, ihi, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmhr_quad(side, trans, m, n, ilo, ihi, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmhr

    subroutine dsyev(jobz, uplo, n, A, lda, w, work, lwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: w(*), work(*)
        integer,   intent(out)   :: info
        call dsyev_quad(jobz, uplo, n, A, lda, w, work, lwork, info)
    end subroutine dsyev

    subroutine dsyevd(jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, lda, lwork, liwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: w(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsyevd_quad(jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info)
    end subroutine dsyevd

    subroutine zheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zheev_quad(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
    end subroutine zheev

    subroutine zheevd(jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, lda, lwork, lrwork, liwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: iwork(*), info
        call zheevd_quad(jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
    end subroutine zheevd

    subroutine dgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, info)
        character, intent(in)    :: jobu, jobvt
        integer,   intent(in)    :: m, n, lda, ldu, ldvt, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*), work(*)
        integer,   intent(out)   :: info
        call dgesvd_quad(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, info)
    end subroutine dgesvd

    subroutine zgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, rwork, info)
        character,   intent(in)    :: jobu, jobvt
        integer,     intent(in)    :: m, n, lda, ldu, ldvt, lwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: s(*), rwork(*)
        complex(ep), intent(out)   :: U(ldu,*), VT(ldvt,*), work(*)
        integer,     intent(out)   :: info
        call zgesvd_quad(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, rwork, info)
    end subroutine zgesvd

    function dlange(norm, m, n, A, lda, work) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: m, n, lda
        real(ep),  intent(in) :: A(lda,*)
        real(ep) :: work(*)
        real(ep) :: r
        r = dlange_quad(norm, m, n, A, lda, work)
    end function dlange

    function dlansy(norm, uplo, n, A, lda, work) result(r)
        character, intent(in) :: norm, uplo
        integer,   intent(in) :: n, lda
        real(ep),  intent(in) :: A(lda,*)
        real(ep) :: work(*)
        real(ep) :: r
        r = dlansy_quad(norm, uplo, n, A, lda, work)
    end function dlansy

    function dlanhs(norm, n, A, lda, work) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: n, lda
        real(ep),  intent(in) :: A(lda,*)
        real(ep) :: work(*)
        real(ep) :: r
        r = dlanhs_quad(norm, n, A, lda, work)
    end function dlanhs

    function dlantr(norm, uplo, diag, m, n, A, lda, work) result(r)
        character, intent(in) :: norm, uplo, diag
        integer,   intent(in) :: m, n, lda
        real(ep),  intent(in) :: A(lda,*)
        real(ep) :: work(*)
        real(ep) :: r
        r = dlantr_quad(norm, uplo, diag, m, n, A, lda, work)
    end function dlantr

    function zlange(norm, m, n, A, lda, work) result(r)
        character,   intent(in) :: norm
        integer,     intent(in) :: m, n, lda
        complex(ep), intent(in) :: A(lda,*)
        real(ep) :: work(*)
        real(ep) :: r
        r = zlange_quad(norm, m, n, A, lda, work)
    end function zlange

    function zlanhe(norm, uplo, n, A, lda, work) result(r)
        character,   intent(in) :: norm, uplo
        integer,     intent(in) :: n, lda
        complex(ep), intent(in) :: A(lda,*)
        real(ep) :: work(*)
        real(ep) :: r
        r = zlanhe_quad(norm, uplo, n, A, lda, work)
    end function zlanhe

    function zlansy(norm, uplo, n, A, lda, work) result(r)
        character,   intent(in) :: norm, uplo
        integer,     intent(in) :: n, lda
        complex(ep), intent(in) :: A(lda,*)
        real(ep) :: work(*)
        real(ep) :: r
        r = zlansy_quad(norm, uplo, n, A, lda, work)
    end function zlansy

    function zlanhs(norm, n, A, lda, work) result(r)
        character,   intent(in) :: norm
        integer,     intent(in) :: n, lda
        complex(ep), intent(in) :: A(lda,*)
        real(ep) :: work(*)
        real(ep) :: r
        r = zlanhs_quad(norm, n, A, lda, work)
    end function zlanhs

    function zlantr(norm, uplo, diag, m, n, A, lda, work) result(r)
        character,   intent(in) :: norm, uplo, diag
        integer,     intent(in) :: m, n, lda
        complex(ep), intent(in) :: A(lda,*)
        real(ep) :: work(*)
        real(ep) :: r
        r = zlantr_quad(norm, uplo, diag, m, n, A, lda, work)
    end function zlantr

    subroutine dlacpy(uplo, m, n, A, lda, B, ldb)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: m, n, lda, ldb
        real(ep),  intent(in)  :: A(lda,*)
        real(ep),  intent(out) :: B(ldb,*)
        call dlacpy_quad(uplo, m, n, A, lda, B, ldb)
    end subroutine dlacpy

    subroutine dlaset(uplo, m, n, alpha, beta, A, lda)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: m, n, lda
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(inout) :: A(lda,*)
        call dlaset_quad(uplo, m, n, alpha, beta, A, lda)
    end subroutine dlaset


end module ref_quad_lapack
