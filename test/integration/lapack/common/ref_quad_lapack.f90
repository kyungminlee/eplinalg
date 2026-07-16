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

        ! ── Linear solve — complex ───────────────────────────────────
        subroutine zgesv_quad(n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgesv_quad

        subroutine zgetrf_quad(m, n, A, lda, ipiv, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgetrf_quad

        subroutine zgetrs_quad(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zgetrs_quad

        subroutine zpotrf_quad(uplo, n, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine zpotrf_quad

        subroutine zpotrs_quad(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpotrs_quad

        subroutine zungqr_quad(m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungqr_quad

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

        subroutine zlacpy_quad(uplo, m, n, A, lda, B, ldb)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: m, n, lda, ldb
            complex(ep), intent(in)  :: A(lda,*)
            complex(ep), intent(out) :: B(ldb,*)
        end subroutine zlacpy_quad

        subroutine zlaset_quad(uplo, m, n, alpha, beta, A, lda)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: m, n, lda
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(inout) :: A(lda,*)
        end subroutine zlaset_quad

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

        ! ── QR factorization — complex ───────────────────────────────
        subroutine zgeqrf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqrf_quad

        ! ── Symmetric / Hermitian eigenvalue ─────────────────────────
        subroutine dsyev_quad(jobz, uplo, n, A, lda, w, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsyev_quad

        subroutine zheev_quad(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zheev_quad

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

        ! ── Auxiliary ────────────────────────────────────────────────
        function dlange_quad(norm, m, n, A, lda, work) result(r)
            import :: ep
            character, intent(in) :: norm
            integer,   intent(in) :: m, n, lda
            real(ep),  intent(in) :: A(lda,*)
            real(ep) :: work(*)
            real(ep) :: r
        end function dlange_quad

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

        ! ── Phase 2 — real user-facing drivers ───────────────────────
        subroutine dposv_quad(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dposv_quad

        subroutine dsysv_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: ipiv(*), info
            real(ep),  intent(out)   :: work(*)
        end subroutine dsysv_quad

        subroutine dsytrf_quad(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: ipiv(*), info
            real(ep),  intent(out)   :: work(*)
        end subroutine dsytrf_quad

        subroutine dsytrs_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dsytrs_quad

        subroutine dgels_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgels_quad

        subroutine dsyevd_quad(jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork, liwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsyevd_quad

        subroutine dgeev_quad(jobvl, jobvr, n, A, lda, wr, wi, &
                         vl, ldvl, vr, ldvr, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobvl, jobvr
            integer,   intent(in)    :: n, lda, ldvl, ldvr, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: wr(*), wi(*), vl(ldvl,*), vr(ldvr,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dgeev_quad

        subroutine dgesdd_quad(jobz, m, n, A, lda, s, U, ldu, VT, ldvt, &
                          work, lwork, iwork, info)
            import :: ep
            character, intent(in)    :: jobz
            integer,   intent(in)    :: m, n, lda, ldu, ldvt, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgesdd_quad

        subroutine dpotri_quad(uplo, n, A, lda, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine dpotri_quad

        subroutine dgetri_quad(n, A, lda, ipiv, work, lwork, info)
            import :: ep
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgetri_quad

        ! ── Phase 3 — triangular + complex driver mirrors ────────────
        subroutine dtrtri_quad(uplo, diag, n, A, lda, info)
            import :: ep
            character, intent(in)    :: uplo, diag
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine dtrtri_quad

        subroutine dtrtrs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dtrtrs_quad

        subroutine ztrtri_quad(uplo, diag, n, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo, diag
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine ztrtri_quad

        subroutine ztrtrs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo, trans, diag
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine ztrtrs_quad

        subroutine zposv_quad(uplo, n, nrhs, A, lda, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zposv_quad

        subroutine zhesv_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zhesv_quad

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

        subroutine zgeev_quad(jobvl, jobvr, n, A, lda, w, vl, ldvl, vr, ldvr, &
                         work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobvl, jobvr
            integer,     intent(in)    :: n, lda, ldvl, ldvr, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: w(*), vl(ldvl,*), vr(ldvr,*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zgeev_quad

        subroutine zgesdd_quad(jobz, m, n, A, lda, s, U, ldu, VT, ldvt, &
                          work, lwork, rwork, iwork, info)
            import :: ep
            character,   intent(in)    :: jobz
            integer,     intent(in)    :: m, n, lda, ldu, ldvt, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: s(*), rwork(*)
            complex(ep), intent(out)   :: U(ldu,*), VT(ldvt,*), work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zgesdd_quad

        subroutine zgels_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgels_quad

        ! ── Phase 4 — selected/generalized eig + QR/LQ family ────────
        subroutine dsyevr_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, isuppz, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, lda, il, iu, ldz, lwork, liwork
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: m, isuppz(*), iwork(*), info
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        end subroutine dsyevr_quad

        subroutine zheevr_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, isuppz, work, lwork, &
                          rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, lda, il, iu, ldz, lwork, lrwork, liwork
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: m, isuppz(*), iwork(*), info
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
        end subroutine zheevr_quad

        subroutine dsygv_quad(itype, jobz, uplo, n, A, lda, B, ldb, w, &
                         work, lwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: itype, n, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsygv_quad

        subroutine zhegv_quad(itype, jobz, uplo, n, A, lda, B, ldb, w, &
                         work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: itype, n, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhegv_quad

        subroutine dgehrd_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: n, ilo, ihi, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgehrd_quad

        subroutine dorghr_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: n, ilo, ihi, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorghr_quad

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

        ! ── Phase 5 — complex mirrors + QR helpers ───────────────────
        subroutine zgetri_quad(n, A, lda, ipiv, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgetri_quad

        subroutine zpotri_quad(uplo, n, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine zpotri_quad

        subroutine zhetrf_quad(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zhetrf_quad

        subroutine zhetrs_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zhetrs_quad

        subroutine zgehrd_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, ilo, ihi, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgehrd_quad

        subroutine zunghr_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, ilo, ihi, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunghr_quad

        subroutine dgeqp3_quad(m, n, A, lda, jpvt, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            integer,  intent(inout) :: jpvt(*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqp3_quad

        subroutine zgeqp3_quad(m, n, A, lda, jpvt, tau, work, lwork, rwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(inout) :: jpvt(*)
            complex(ep), intent(out)   :: tau(*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zgeqp3_quad

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

        ! ── Phase 6 — banded + tridiagonal ───────────────────────────
        subroutine dgbtrf_quad(m, n, kl, ku, AB, ldab, ipiv, info)
            import :: ep
            integer,  intent(in)    :: m, n, kl, ku, ldab
            real(ep), intent(inout) :: AB(ldab,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine dgbtrf_quad

        subroutine dgbtrs_quad(trans, n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldb
            real(ep),  intent(in)    :: AB(ldab,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dgbtrs_quad

        subroutine dgbsv_quad(n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, kl, ku, nrhs, ldab, ldb
            real(ep), intent(inout) :: AB(ldab,*), B(ldb,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine dgbsv_quad

        subroutine dpbtrf_quad(uplo, n, kd, AB, ldab, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, kd, ldab
            real(ep),  intent(inout) :: AB(ldab,*)
            integer,   intent(out)   :: info
        end subroutine dpbtrf_quad

        subroutine dpbtrs_quad(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, kd, nrhs, ldab, ldb
            real(ep),  intent(in)    :: AB(ldab,*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dpbtrs_quad

        subroutine dpbsv_quad(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, kd, nrhs, ldab, ldb
            real(ep),  intent(inout) :: AB(ldab,*), B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dpbsv_quad

        subroutine dgttrf_quad(n, dl, d, du, du2, ipiv, info)
            import :: ep
            integer,  intent(in)    :: n
            real(ep), intent(inout) :: dl(*), d(*), du(*)
            real(ep), intent(out)   :: du2(*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine dgttrf_quad

        subroutine dgttrs_quad(trans, n, nrhs, dl, d, du, du2, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(in)    :: dl(*), d(*), du(*), du2(*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dgttrs_quad

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
        ! ── Phase 7 — symmetric / Hermitian eigenvalue family ────────
        subroutine dsyevx_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, work, lwork, &
                          iwork, ifail, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, lda, il, iu, ldz, lwork
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        end subroutine dsyevx_quad

        subroutine zheevx_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, work, lwork, &
                          rwork, iwork, ifail, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, lda, il, iu, ldz, lwork
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: m, iwork(*), ifail(*), info
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
        end subroutine zheevx_quad

        subroutine dstev_quad(jobz, n, d, e, z, ldz, work, info)
            import :: ep
            character, intent(in)    :: jobz
            integer,   intent(in)    :: n, ldz
            real(ep),  intent(inout) :: d(*), e(*)
            real(ep),  intent(out)   :: z(ldz,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dstev_quad

        subroutine dstevd_quad(jobz, n, d, e, z, ldz, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz
            integer,   intent(in)    :: n, ldz, lwork, liwork
            real(ep),  intent(inout) :: d(*), e(*)
            real(ep),  intent(out)   :: z(ldz,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dstevd_quad

        subroutine dstevx_quad(jobz, range, n, d, e, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, work, iwork, ifail, info)
            import :: ep
            character, intent(in)    :: jobz, range
            integer,   intent(in)    :: n, il, iu, ldz
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: d(*), e(*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        end subroutine dstevx_quad

        subroutine dstevr_quad(jobz, range, n, d, e, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, isuppz, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, range
            integer,   intent(in)    :: n, il, iu, ldz, lwork, liwork
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: d(*), e(*)
            integer,   intent(out)   :: m, isuppz(*), iwork(*), info
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        end subroutine dstevr_quad

        subroutine dsbev_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, kd, ldab, ldz
            real(ep),  intent(inout) :: AB(ldab,*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsbev_quad

        subroutine dsbevd_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, &
                          work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, kd, ldab, ldz, lwork, liwork
            real(ep),  intent(inout) :: AB(ldab,*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsbevd_quad

        subroutine dsbevx_quad(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, &
                          vl, vu, il, iu, abstol, m, w, z, ldz, &
                          work, iwork, ifail, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, kd, ldab, ldq, il, iu, ldz
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: AB(ldab,*)
            real(ep),  intent(out)   :: Q(ldq,*), w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine dsbevx_quad

        subroutine zhbev_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, &
                         work, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, kd, ldab, ldz
            complex(ep), intent(inout) :: AB(ldab,*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zhbev_quad

        subroutine zhbevd_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, &
                          work, lwork, rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, kd, ldab, ldz, lwork, lrwork, liwork
            complex(ep), intent(inout) :: AB(ldab,*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zhbevd_quad

        subroutine zhbevx_quad(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, &
                          vl, vu, il, iu, abstol, m, w, z, ldz, &
                          work, rwork, iwork, ifail, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, kd, ldab, ldq, il, iu, ldz
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: AB(ldab,*)
            complex(ep), intent(out)   :: Q(ldq,*), z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine zhbevx_quad

        subroutine dspev_quad(jobz, uplo, n, AP, w, z, ldz, work, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ldz
            real(ep),  intent(inout) :: AP(*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dspev_quad

        subroutine dspevd_quad(jobz, uplo, n, AP, w, z, ldz, &
                          work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ldz, lwork, liwork
            real(ep),  intent(inout) :: AP(*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dspevd_quad

        subroutine dspevx_quad(jobz, range, uplo, n, AP, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, work, iwork, ifail, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, il, iu, ldz
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        end subroutine dspevx_quad

        subroutine zhpev_quad(jobz, uplo, n, AP, w, z, ldz, work, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, ldz
            complex(ep), intent(inout) :: AP(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zhpev_quad

        subroutine zhpevd_quad(jobz, uplo, n, AP, w, z, ldz, &
                          work, lwork, rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, ldz, lwork, lrwork, liwork
            complex(ep), intent(inout) :: AP(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zhpevd_quad

        subroutine zhpevx_quad(jobz, range, uplo, n, AP, vl, vu, il, iu, &
                          abstol, m, w, z, ldz, work, rwork, &
                          iwork, ifail, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, il, iu, ldz
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: AP(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine zhpevx_quad

        ! ── Phase 8 — QL / RQ family ─────────────────────────────────
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

        ! ── Phase 9 — real packed sym factor / solve / inverse ───────
        subroutine dsptrf_quad(uplo, n, AP, ipiv, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(out)   :: ipiv(*), info
        end subroutine dsptrf_quad

        subroutine dsptrs_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(in)    :: AP(*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dsptrs_quad

        subroutine dspsv_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(inout) :: AP(*), B(ldb,*)
            integer,   intent(out)   :: ipiv(*), info
        end subroutine dspsv_quad

        subroutine dsptri_quad(uplo, n, AP, ipiv, work, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dsptri_quad

        subroutine dpptrf_quad(uplo, n, AP, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(out)   :: info
        end subroutine dpptrf_quad

        subroutine dpptrs_quad(uplo, n, nrhs, AP, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(in)    :: AP(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dpptrs_quad

        subroutine dppsv_quad(uplo, n, nrhs, AP, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(inout) :: AP(*), B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dppsv_quad

        subroutine dpptri_quad(uplo, n, AP, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(out)   :: info
        end subroutine dpptri_quad

        ! ── Phase 10 — complex packed factor / solve / inverse ───────
        subroutine zhptrf_quad(uplo, n, AP, ipiv, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zhptrf_quad

        subroutine zhptrs_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(in)    :: AP(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zhptrs_quad

        subroutine zhpsv_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(inout) :: AP(*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zhpsv_quad

        subroutine zhptri_quad(uplo, n, AP, ipiv, work, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhptri_quad

        subroutine zsptrf_quad(uplo, n, AP, ipiv, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zsptrf_quad

        subroutine zsptrs_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(in)    :: AP(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zsptrs_quad

        subroutine zspsv_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(inout) :: AP(*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zspsv_quad

        subroutine zsptri_quad(uplo, n, AP, ipiv, work, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsptri_quad

        subroutine zpptrf_quad(uplo, n, AP, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(out)   :: info
        end subroutine zpptrf_quad

        subroutine zpptrs_quad(uplo, n, nrhs, AP, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(in)    :: AP(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpptrs_quad

        subroutine zppsv_quad(uplo, n, nrhs, AP, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(inout) :: AP(*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zppsv_quad

        subroutine zpptri_quad(uplo, n, AP, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(out)   :: info
        end subroutine zpptri_quad

        ! ── Phase 11 — complex banded / tridiagonal ──────────────────
        subroutine zgbtrf_quad(m, n, kl, ku, AB, ldab, ipiv, info)
            import :: ep
            integer,     intent(in)    :: m, n, kl, ku, ldab
            complex(ep), intent(inout) :: AB(ldab,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgbtrf_quad

        subroutine zgbtrs_quad(trans, n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldb
            complex(ep), intent(in)    :: AB(ldab,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zgbtrs_quad

        subroutine zgbsv_quad(n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldb
            complex(ep), intent(inout) :: AB(ldab,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgbsv_quad

        subroutine zpbtrf_quad(uplo, n, kd, AB, ldab, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, kd, ldab
            complex(ep), intent(inout) :: AB(ldab,*)
            integer,     intent(out)   :: info
        end subroutine zpbtrf_quad

        subroutine zpbtrs_quad(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, kd, nrhs, ldab, ldb
            complex(ep), intent(in)    :: AB(ldab,*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpbtrs_quad

        subroutine zpbsv_quad(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, kd, nrhs, ldab, ldb
            complex(ep), intent(inout) :: AB(ldab,*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpbsv_quad

        subroutine zgttrf_quad(n, dl, d, du, du2, ipiv, info)
            import :: ep
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: dl(*), d(*), du(*)
            complex(ep), intent(out)   :: du2(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgttrf_quad

        subroutine zgttrs_quad(trans, n, nrhs, dl, d, du, du2, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(in)    :: dl(*), d(*), du(*), du2(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zgttrs_quad

        subroutine zgtsv_quad(n, nrhs, dl, d, du, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(inout) :: dl(*), d(*), du(*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zgtsv_quad

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

        subroutine zptsv_quad(n, nrhs, d, e, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, ldb
            real(ep),    intent(inout) :: d(*)
            complex(ep), intent(inout) :: e(*), B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zptsv_quad

        ! ── Phase 12 — triangular packed / banded ────────────────────
        subroutine dtbtrs_quad(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, kd, nrhs, ldab, ldb
            real(ep),  intent(in)    :: AB(ldab,*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dtbtrs_quad

        subroutine ztbtrs_quad(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo, trans, diag
            integer,     intent(in)    :: n, kd, nrhs, ldab, ldb
            complex(ep), intent(in)    :: AB(ldab,*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine ztbtrs_quad

        subroutine dtptrs_quad(uplo, trans, diag, n, nrhs, AP, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(in)    :: AP(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dtptrs_quad

        subroutine ztptrs_quad(uplo, trans, diag, n, nrhs, AP, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo, trans, diag
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(in)    :: AP(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine ztptrs_quad

        subroutine dtptri_quad(uplo, diag, n, AP, info)
            import :: ep
            character, intent(in)    :: uplo, diag
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            integer,   intent(out)   :: info
        end subroutine dtptri_quad

        subroutine ztptri_quad(uplo, diag, n, AP, info)
            import :: ep
            character,   intent(in)    :: uplo, diag
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            integer,     intent(out)   :: info
        end subroutine ztptri_quad

        ! ── Phase 13 — condition number estimators ───────────────────
        subroutine dgecon_quad(norm, n, A, lda, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: norm
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*), anorm
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dgecon_quad

        subroutine zgecon_quad(norm, n, A, lda, anorm, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: norm
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zgecon_quad

        subroutine dpocon_quad(uplo, n, A, lda, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*), anorm
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dpocon_quad

        subroutine zpocon_quad(uplo, n, A, lda, anorm, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zpocon_quad

        subroutine dgbcon_quad(norm, n, kl, ku, AB, ldab, ipiv, anorm, rcond, &
                          work, iwork, info)
            import :: ep
            character, intent(in)  :: norm
            integer,   intent(in)  :: n, kl, ku, ldab
            real(ep),  intent(in)  :: AB(ldab,*), anorm
            integer,   intent(in)  :: ipiv(*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dgbcon_quad

        subroutine zgbcon_quad(norm, n, kl, ku, AB, ldab, ipiv, anorm, rcond, &
                          work, rwork, info)
            import :: ep
            character,   intent(in)  :: norm
            integer,     intent(in)  :: n, kl, ku, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zgbcon_quad

        subroutine dgtcon_quad(norm, n, dl, d, du, du2, ipiv, anorm, rcond, &
                          work, iwork, info)
            import :: ep
            character, intent(in)  :: norm
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: dl(*), d(*), du(*), du2(*), anorm
            integer,   intent(in)  :: ipiv(*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dgtcon_quad

        subroutine zgtcon_quad(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: norm
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: dl(*), d(*), du(*), du2(*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zgtcon_quad

        subroutine dsycon_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*), anorm
            integer,   intent(in)  :: ipiv(*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dsycon_quad

        subroutine zhecon_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zhecon_quad

        ! ── Phase 14 — condition estimators part 2 ───────────────────
        subroutine dpbcon_quad(uplo, n, kd, AB, ldab, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, kd, ldab
            real(ep),  intent(in)  :: AB(ldab,*), anorm
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dpbcon_quad

        subroutine zpbcon_quad(uplo, n, kd, AB, ldab, anorm, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, kd, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zpbcon_quad

        subroutine dppcon_quad(uplo, n, AP, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: AP(*), anorm
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dppcon_quad

        subroutine zppcon_quad(uplo, n, AP, anorm, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: AP(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zppcon_quad

        subroutine dptcon_quad(n, d, e, anorm, rcond, work, info)
            import :: ep
            integer,  intent(in)  :: n
            real(ep), intent(in)  :: d(*), e(*), anorm
            real(ep), intent(out) :: rcond, work(*)
            integer,  intent(out) :: info
        end subroutine dptcon_quad

        subroutine zptcon_quad(n, d, e, anorm, rcond, rwork, info)
            import :: ep
            integer,     intent(in)  :: n
            real(ep),    intent(in)  :: d(*), anorm
            complex(ep), intent(in)  :: e(*)
            real(ep),    intent(out) :: rcond, rwork(*)
            integer,     intent(out) :: info
        end subroutine zptcon_quad

        subroutine dspcon_quad(uplo, n, AP, ipiv, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: AP(*), anorm
            integer,   intent(in)  :: ipiv(*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dspcon_quad

        subroutine zhpcon_quad(uplo, n, AP, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: AP(*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zhpcon_quad

        subroutine zspcon_quad(uplo, n, AP, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: AP(*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zspcon_quad

        subroutine dtrcon_quad(norm, uplo, diag, n, A, lda, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: norm, uplo, diag
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtrcon_quad

        subroutine ztrcon_quad(norm, uplo, diag, n, A, lda, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: norm, uplo, diag
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztrcon_quad

        subroutine dtpcon_quad(norm, uplo, diag, n, AP, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: norm, uplo, diag
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: AP(*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtpcon_quad

        subroutine ztpcon_quad(norm, uplo, diag, n, AP, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: norm, uplo, diag
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: AP(*)
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztpcon_quad

        subroutine dtbcon_quad(norm, uplo, diag, n, kd, AB, ldab, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: norm, uplo, diag
            integer,   intent(in)  :: n, kd, ldab
            real(ep),  intent(in)  :: AB(ldab,*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtbcon_quad

        subroutine ztbcon_quad(norm, uplo, diag, n, kd, AB, ldab, rcond, work, rwork, info)
            import :: ep
            character,   intent(in)  :: norm, uplo, diag
            integer,     intent(in)  :: n, kd, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            real(ep),    intent(out) :: rcond, rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztbcon_quad

        ! ── Phase 15 — iterative refinement ──────────────────────────
        subroutine dgerfs_quad(trans, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgerfs_quad

        subroutine zgerfs_quad(trans, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgerfs_quad

        subroutine dporfs_quad(uplo, n, nrhs, A, lda, AF, ldaf, B, ldb, &
                          X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dporfs_quad

        subroutine zporfs_quad(uplo, n, nrhs, A, lda, AF, ldaf, B, ldb, &
                          X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zporfs_quad

        ! ── Phase 16 — banded / tridiag refinement ───────────────────
        subroutine dgbrfs_quad(trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, &
                          B, ldb, X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx
            real(ep),  intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgbrfs_quad

        subroutine zgbrfs_quad(trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, &
                          B, ldb, X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx
            complex(ep), intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgbrfs_quad

        subroutine dpbrfs_quad(uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, B, ldb, &
                          X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, kd, nrhs, ldab, ldafb, ldb, ldx
            real(ep),  intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dpbrfs_quad

        subroutine zpbrfs_quad(uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, B, ldb, &
                          X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, kd, nrhs, ldab, ldafb, ldb, ldx
            complex(ep), intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zpbrfs_quad

        subroutine dgtrfs_quad(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, &
                          B, ldb, X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: n, nrhs, ldb, ldx
            real(ep),  intent(in)    :: dl(*), d(*), du(*), dlf(*), df(*), duf(*), du2(*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(in)    :: B(ldb,*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgtrfs_quad

        subroutine zgtrfs_quad(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, &
                          B, ldb, X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            complex(ep), intent(in)    :: dl(*), d(*), du(*), dlf(*), df(*), duf(*), du2(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(in)    :: B(ldb,*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgtrfs_quad

        subroutine dptrfs_quad(n, nrhs, d, e, df, ef, B, ldb, X, ldx, &
                          ferr, berr, work, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, ldb, ldx
            real(ep), intent(in)    :: d(*), e(*), df(*), ef(*), B(ldb,*)
            real(ep), intent(inout) :: X(ldx,*)
            real(ep), intent(out)   :: ferr(*), berr(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dptrfs_quad

        subroutine zptrfs_quad(uplo, n, nrhs, d, e, df, ef, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            real(ep),    intent(in)    :: d(*), df(*)
            complex(ep), intent(in)    :: e(*), ef(*), B(ldb,*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zptrfs_quad

        ! ── Phase 17 — sym/Herm + packed refinement ──────────────────
        subroutine dsyrfs_quad(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsyrfs_quad

        subroutine zsyrfs_quad(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsyrfs_quad

        subroutine zherfs_quad(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zherfs_quad

        subroutine dsprfs_quad(uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, &
                          ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb, ldx
            real(ep),  intent(in)    :: AP(*), AFP(*), B(ldb,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsprfs_quad

        subroutine zsprfs_quad(uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            complex(ep), intent(in)    :: AP(*), AFP(*), B(ldb,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsprfs_quad

        subroutine dpprfs_quad(uplo, n, nrhs, AP, AFP, B, ldb, X, ldx, &
                          ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, ldb, ldx
            real(ep),  intent(in)    :: AP(*), AFP(*), B(ldb,*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dpprfs_quad

        subroutine zpprfs_quad(uplo, n, nrhs, AP, AFP, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            complex(ep), intent(in)    :: AP(*), AFP(*), B(ldb,*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zpprfs_quad

        ! ── Phase 18 — triangular refinement ─────────────────────────
        subroutine dtrrfs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, X, ldx, &
                          ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo, trans, diag
            integer,   intent(in)  :: n, nrhs, lda, ldb, ldx
            real(ep),  intent(in)  :: A(lda,*), B(ldb,*), X(ldx,*)
            real(ep),  intent(out) :: ferr(*), berr(*), work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtrrfs_quad

        subroutine ztrrfs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo, trans, diag
            integer,     intent(in)  :: n, nrhs, lda, ldb, ldx
            complex(ep), intent(in)  :: A(lda,*), B(ldb,*), X(ldx,*)
            real(ep),    intent(out) :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztrrfs_quad

        subroutine dtprfs_quad(uplo, trans, diag, n, nrhs, AP, B, ldb, X, ldx, &
                          ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo, trans, diag
            integer,   intent(in)  :: n, nrhs, ldb, ldx
            real(ep),  intent(in)  :: AP(*), B(ldb,*), X(ldx,*)
            real(ep),  intent(out) :: ferr(*), berr(*), work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtprfs_quad

        subroutine ztprfs_quad(uplo, trans, diag, n, nrhs, AP, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo, trans, diag
            integer,     intent(in)  :: n, nrhs, ldb, ldx
            complex(ep), intent(in)  :: AP(*), B(ldb,*), X(ldx,*)
            real(ep),    intent(out) :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztprfs_quad

        subroutine dtbrfs_quad(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, X, ldx, &
                          ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo, trans, diag
            integer,   intent(in)  :: n, kd, nrhs, ldab, ldb, ldx
            real(ep),  intent(in)  :: AB(ldab,*), B(ldb,*), X(ldx,*)
            real(ep),  intent(out) :: ferr(*), berr(*), work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dtbrfs_quad

        subroutine ztbrfs_quad(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, X, ldx, &
                          ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)  :: uplo, trans, diag
            integer,     intent(in)  :: n, kd, nrhs, ldab, ldb, ldx
            complex(ep), intent(in)  :: AB(ldab,*), B(ldb,*), X(ldx,*)
            real(ep),    intent(out) :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine ztbrfs_quad

        ! ── Equilibration — real ─────────────────────────────────────
        subroutine dgeequ_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,  intent(in)  :: m, n, lda
            real(ep), intent(in)  :: A(lda,*)
            real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,  intent(out) :: info
        end subroutine dgeequ_quad

        subroutine dgbequ_quad(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,  intent(in)  :: m, n, kl, ku, ldab
            real(ep), intent(in)  :: AB(ldab,*)
            real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,  intent(out) :: info
        end subroutine dgbequ_quad

        subroutine dpoequ_quad(n, A, lda, S, scond, amax, info)
            import :: ep
            integer,  intent(in)  :: n, lda
            real(ep), intent(in)  :: A(lda,*)
            real(ep), intent(out) :: S(*), scond, amax
            integer,  intent(out) :: info
        end subroutine dpoequ_quad

        subroutine dppequ_quad(uplo, n, AP, S, scond, amax, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: AP(*)
            real(ep),  intent(out) :: S(*), scond, amax
            integer,   intent(out) :: info
        end subroutine dppequ_quad

        subroutine dpbequ_quad(uplo, n, kd, AB, ldab, S, scond, amax, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, kd, ldab
            real(ep),  intent(in)  :: AB(ldab,*)
            real(ep),  intent(out) :: S(*), scond, amax
            integer,   intent(out) :: info
        end subroutine dpbequ_quad

        ! ── Equilibration — complex ──────────────────────────────────
        subroutine zgeequ_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,     intent(in)  :: m, n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,     intent(out) :: info
        end subroutine zgeequ_quad

        subroutine zgbequ_quad(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,     intent(in)  :: m, n, kl, ku, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,     intent(out) :: info
        end subroutine zgbequ_quad

        subroutine zpoequ_quad(n, A, lda, S, scond, amax, info)
            import :: ep
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: S(*), scond, amax
            integer,     intent(out) :: info
        end subroutine zpoequ_quad

        subroutine zppequ_quad(uplo, n, AP, S, scond, amax, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: AP(*)
            real(ep),    intent(out) :: S(*), scond, amax
            integer,     intent(out) :: info
        end subroutine zppequ_quad

        subroutine zpbequ_quad(uplo, n, kd, AB, ldab, S, scond, amax, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, kd, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            real(ep),    intent(out) :: S(*), scond, amax
            integer,     intent(out) :: info
        end subroutine zpbequ_quad

        ! ── Improved equilibration (equb) ────────────────────────────
        subroutine dgeequb_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,  intent(in)  :: m, n, lda
            real(ep), intent(in)  :: A(lda,*)
            real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,  intent(out) :: info
        end subroutine dgeequb_quad

        subroutine dgbequb_quad(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,  intent(in)  :: m, n, kl, ku, ldab
            real(ep), intent(in)  :: AB(ldab,*)
            real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,  intent(out) :: info
        end subroutine dgbequb_quad

        subroutine dpoequb_quad(n, A, lda, S, scond, amax, info)
            import :: ep
            integer,  intent(in)  :: n, lda
            real(ep), intent(in)  :: A(lda,*)
            real(ep), intent(out) :: S(*), scond, amax
            integer,  intent(out) :: info
        end subroutine dpoequb_quad

        subroutine dsyequb_quad(uplo, n, A, lda, S, scond, amax, work, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*)
            real(ep),  intent(out) :: S(*), scond, amax, work(*)
            integer,   intent(out) :: info
        end subroutine dsyequb_quad

        subroutine zgeequb_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,     intent(in)  :: m, n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,     intent(out) :: info
        end subroutine zgeequb_quad

        subroutine zgbequb_quad(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
            import :: ep
            integer,     intent(in)  :: m, n, kl, ku, ldab
            complex(ep), intent(in)  :: AB(ldab,*)
            real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
            integer,     intent(out) :: info
        end subroutine zgbequb_quad

        subroutine zpoequb_quad(n, A, lda, S, scond, amax, info)
            import :: ep
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: S(*), scond, amax
            integer,     intent(out) :: info
        end subroutine zpoequb_quad

        subroutine zsyequb_quad(uplo, n, A, lda, S, scond, amax, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: S(*), scond, amax
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zsyequb_quad

        subroutine zheequb_quad(uplo, n, A, lda, S, scond, amax, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            real(ep),    intent(out) :: S(*), scond, amax
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zheequb_quad

        ! ── Tridiagonal reduction ────────────────────────────────────
        subroutine dsytrd_quad(uplo, n, A, lda, D, E, tau, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: D(*), E(*), tau(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsytrd_quad

        subroutine zhetrd_quad(uplo, n, A, lda, D, E, tau, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: D(*), E(*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhetrd_quad

        subroutine dsptrd_quad(uplo, n, AP, D, E, tau, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: AP(*)
            real(ep),  intent(out)   :: D(*), E(*), tau(*)
            integer,   intent(out)   :: info
        end subroutine dsptrd_quad

        subroutine zhptrd_quad(uplo, n, AP, D, E, tau, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: AP(*)
            real(ep),    intent(out)   :: D(*), E(*)
            complex(ep), intent(out)   :: tau(*)
            integer,     intent(out)   :: info
        end subroutine zhptrd_quad

        subroutine dsbtrd_quad(vect, uplo, n, kd, AB, ldab, D, E, Q, ldq, work, info)
            import :: ep
            character, intent(in)    :: vect, uplo
            integer,   intent(in)    :: n, kd, ldab, ldq
            real(ep),  intent(inout) :: AB(ldab,*), Q(ldq,*)
            real(ep),  intent(out)   :: D(*), E(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsbtrd_quad

        subroutine zhbtrd_quad(vect, uplo, n, kd, AB, ldab, D, E, Q, ldq, work, info)
            import :: ep
            character,   intent(in)    :: vect, uplo
            integer,     intent(in)    :: n, kd, ldab, ldq
            complex(ep), intent(inout) :: AB(ldab,*), Q(ldq,*)
            real(ep),    intent(out)   :: D(*), E(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhbtrd_quad

        ! ── Orthogonal/unitary Q generation/application for *trd ─────
        subroutine dorgtr_quad(uplo, n, A, lda, tau, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(in)    :: tau(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dorgtr_quad

        subroutine zungtr_quad(uplo, n, A, lda, tau, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungtr_quad

        subroutine dormtr_quad(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
            import :: ep
            character, intent(in)    :: side, uplo, trans
            integer,   intent(in)    :: m, n, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormtr_quad

        subroutine zunmtr_quad(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, uplo, trans
            integer,     intent(in)    :: m, n, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmtr_quad

        subroutine dopgtr_quad(uplo, n, AP, tau, Q, ldq, work, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, ldq
            real(ep),  intent(in)  :: AP(*), tau(*)
            real(ep),  intent(out) :: Q(ldq,*), work(*)
            integer,   intent(out) :: info
        end subroutine dopgtr_quad

        subroutine zupgtr_quad(uplo, n, AP, tau, Q, ldq, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, ldq
            complex(ep), intent(in)  :: AP(*), tau(*)
            complex(ep), intent(out) :: Q(ldq,*), work(*)
            integer,     intent(out) :: info
        end subroutine zupgtr_quad

        subroutine dopmtr_quad(side, uplo, trans, m, n, AP, tau, C, ldc, work, info)
            import :: ep
            character, intent(in)    :: side, uplo, trans
            integer,   intent(in)    :: m, n, ldc
            real(ep),  intent(in)    :: AP(*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dopmtr_quad

        subroutine zupmtr_quad(side, uplo, trans, m, n, AP, tau, C, ldc, work, info)
            import :: ep
            character,   intent(in)    :: side, uplo, trans
            integer,     intent(in)    :: m, n, ldc
            complex(ep), intent(in)    :: AP(*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zupmtr_quad

        ! ── Bidiagonal reduction ─────────────────────────────────────
        subroutine dgebrd_quad(m, n, A, lda, D, E, tauq, taup, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: D(*), E(*), tauq(*), taup(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgebrd_quad

        subroutine zgebrd_quad(m, n, A, lda, D, E, tauq, taup, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: D(*), E(*)
            complex(ep), intent(out)   :: tauq(*), taup(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgebrd_quad

        subroutine dorgbr_quad(vect, m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            character, intent(in)    :: vect
            integer,   intent(in)    :: m, n, k, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(in)    :: tau(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dorgbr_quad

        subroutine zungbr_quad(vect, m, n, k, A, lda, tau, work, lwork, info)
            import :: ep
            character,   intent(in)    :: vect
            integer,     intent(in)    :: m, n, k, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungbr_quad

        subroutine dormbr_quad(vect, side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: vect, side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormbr_quad

        subroutine zunmbr_quad(vect, side, trans, m, n, k, A, lda, tau, C, ldc, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: vect, side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmbr_quad

        ! ── Generalized eigenvalue (packed / banded) ─────────────────
        subroutine dspgv_quad(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: itype, n, ldz
            real(ep),  intent(inout) :: AP(*), BP(*)
            real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dspgv_quad

        subroutine zhpgv_quad(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: itype, n, ldz
            complex(ep), intent(inout) :: AP(*), BP(*)
            real(ep),    intent(out)   :: W(*), rwork(*)
            complex(ep), intent(out)   :: Z(ldz,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhpgv_quad

        subroutine dsbgv_quad(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, &
                         work, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ka, kb, ldab, ldbb, ldz
            real(ep),  intent(inout) :: AB(ldab,*), BB(ldbb,*)
            real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsbgv_quad

        subroutine zhbgv_quad(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, &
                         work, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, ka, kb, ldab, ldbb, ldz
            complex(ep), intent(inout) :: AB(ldab,*), BB(ldbb,*)
            real(ep),    intent(out)   :: W(*), rwork(*)
            complex(ep), intent(out)   :: Z(ldz,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhbgv_quad

        ! ── Generalized eigenvalue D&C variants ──────────────────────
        subroutine dsygvd_quad(itype, jobz, uplo, n, A, lda, B, ldb, W, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: itype, n, lda, ldb, lwork, liwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: W(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsygvd_quad

        subroutine zhegvd_quad(itype, jobz, uplo, n, A, lda, B, ldb, W, work, lwork, &
                          rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: itype, n, lda, ldb, lwork, lrwork, liwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),    intent(out)   :: W(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zhegvd_quad

        subroutine dspgvd_quad(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, lwork, &
                          iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: itype, n, ldz, lwork, liwork
            real(ep),  intent(inout) :: AP(*), BP(*)
            real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dspgvd_quad

        subroutine zhpgvd_quad(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, lwork, &
                          rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: itype, n, ldz, lwork, lrwork, liwork
            complex(ep), intent(inout) :: AP(*), BP(*)
            real(ep),    intent(out)   :: W(*), rwork(*)
            complex(ep), intent(out)   :: Z(ldz,*), work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zhpgvd_quad

        subroutine dsbgvd_quad(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, &
                          work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, ka, kb, ldab, ldbb, ldz, lwork, liwork
            real(ep),  intent(inout) :: AB(ldab,*), BB(ldbb,*)
            real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsbgvd_quad

        subroutine zhbgvd_quad(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, &
                          work, lwork, rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, ka, kb, ldab, ldbb, ldz, lwork, &
                                          lrwork, liwork
            complex(ep), intent(inout) :: AB(ldab,*), BB(ldbb,*)
            real(ep),    intent(out)   :: W(*), rwork(*)
            complex(ep), intent(out)   :: Z(ldz,*), work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zhbgvd_quad

        ! ── Tridiagonal eigenvalue solvers ───────────────────────────
        subroutine dsterf_quad(n, D, E, info)
            import :: ep
            integer,  intent(in)    :: n
            real(ep), intent(inout) :: D(*), E(*)
            integer,  intent(out)   :: info
        end subroutine dsterf_quad

        subroutine dsteqr_quad(compz, n, D, E, Z, ldz, work, info)
            import :: ep
            character, intent(in)    :: compz
            integer,   intent(in)    :: n, ldz
            real(ep),  intent(inout) :: D(*), E(*), Z(ldz,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dsteqr_quad

        subroutine zsteqr_quad(compz, n, D, E, Z, ldz, work, info)
            import :: ep
            character,   intent(in)    :: compz
            integer,     intent(in)    :: n, ldz
            real(ep),    intent(inout) :: D(*), E(*)
            complex(ep), intent(inout) :: Z(ldz,*)
            real(ep),    intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsteqr_quad

        subroutine dstedc_quad(compz, n, D, E, Z, ldz, work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: compz
            integer,   intent(in)    :: n, ldz, lwork, liwork
            real(ep),  intent(inout) :: D(*), E(*), Z(ldz,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dstedc_quad

        subroutine zstedc_quad(compz, n, D, E, Z, ldz, work, lwork, rwork, lrwork, &
                          iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: compz
            integer,     intent(in)    :: n, ldz, lwork, lrwork, liwork
            real(ep),    intent(inout) :: D(*), E(*)
            complex(ep), intent(inout) :: Z(ldz,*)
            real(ep),    intent(out)   :: rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zstedc_quad

        subroutine dstebz_quad(range, order, n, vl, vu, il, iu, abstol, D, E, &
                          m, nsplit, W, iblock, isplit, work, iwork, info)
            import :: ep
            character, intent(in)  :: range, order
            integer,   intent(in)  :: n, il, iu
            real(ep),  intent(in)  :: vl, vu, abstol, D(*), E(*)
            integer,   intent(out) :: m, nsplit
            real(ep),  intent(out) :: W(*), work(*)
            integer,   intent(out) :: iblock(*), isplit(*), iwork(*), info
        end subroutine dstebz_quad

        ! ── Bidiagonal SVD ───────────────────────────────────────────
        subroutine dbdsqr_quad(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, &
                          C, ldc, work, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, ncvt, nru, ncc, ldvt, ldu, ldc
            real(ep),  intent(inout) :: D(*), E(*), VT(ldvt,*), U(ldu,*), C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dbdsqr_quad

        subroutine zbdsqr_quad(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, &
                          C, ldc, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, ncvt, nru, ncc, ldvt, ldu, ldc
            real(ep),    intent(inout) :: D(*), E(*)
            complex(ep), intent(inout) :: VT(ldvt,*), U(ldu,*), C(ldc,*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zbdsqr_quad

        subroutine dbdsdc_quad(uplo, compq, n, D, E, U, ldu, VT, ldvt, Q, IQ, &
                          work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo, compq
            integer,   intent(in)    :: n, ldu, ldvt
            real(ep),  intent(inout) :: D(*), E(*)
            real(ep),  intent(out)   :: U(ldu,*), VT(ldvt,*), Q(*), work(*)
            integer,   intent(out)   :: IQ(*), iwork(*), info
        end subroutine dbdsdc_quad

        ! ── ddisna + balancing/back-transform ────────────────────────
        subroutine ddisna_quad(job, m, n, D, sep, info)
            import :: ep
            character, intent(in)  :: job
            integer,   intent(in)  :: m, n
            real(ep),  intent(in)  :: D(*)
            real(ep),  intent(out) :: sep(*)
            integer,   intent(out) :: info
        end subroutine ddisna_quad

        subroutine dgebal_quad(job, n, A, lda, ilo, ihi, scale, info)
            import :: ep
            character, intent(in)    :: job
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: ilo, ihi
            real(ep),  intent(out)   :: scale(*)
            integer,   intent(out)   :: info
        end subroutine dgebal_quad

        subroutine zgebal_quad(job, n, A, lda, ilo, ihi, scale, info)
            import :: ep
            character,   intent(in)    :: job
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ilo, ihi
            real(ep),    intent(out)   :: scale(*)
            integer,     intent(out)   :: info
        end subroutine zgebal_quad

        subroutine dgebak_quad(job, side, n, ilo, ihi, scale, m, V, ldv, info)
            import :: ep
            character, intent(in)    :: job, side
            integer,   intent(in)    :: n, ilo, ihi, m, ldv
            real(ep),  intent(in)    :: scale(*)
            real(ep),  intent(inout) :: V(ldv,*)
            integer,   intent(out)   :: info
        end subroutine dgebak_quad

        subroutine zgebak_quad(job, side, n, ilo, ihi, scale, m, V, ldv, info)
            import :: ep
            character,   intent(in)    :: job, side
            integer,     intent(in)    :: n, ilo, ihi, m, ldv
            real(ep),    intent(in)    :: scale(*)
            complex(ep), intent(inout) :: V(ldv,*)
            integer,     intent(out)   :: info
        end subroutine zgebak_quad

        ! ── Generalized balance/back-transform/Hessenberg ────────────
        subroutine dggbal_quad(job, n, A, lda, B, ldb, ilo, ihi, lscale, rscale, &
                          work, info)
            import :: ep
            character, intent(in)    :: job
            integer,   intent(in)    :: n, lda, ldb
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: ilo, ihi
            real(ep),  intent(out)   :: lscale(*), rscale(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dggbal_quad

        subroutine zggbal_quad(job, n, A, lda, B, ldb, ilo, ihi, lscale, rscale, &
                          work, info)
            import :: ep
            character,   intent(in)    :: job
            integer,     intent(in)    :: n, lda, ldb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ilo, ihi
            real(ep),    intent(out)   :: lscale(*), rscale(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zggbal_quad

        subroutine dggbak_quad(job, side, n, ilo, ihi, lscale, rscale, m, V, ldv, info)
            import :: ep
            character, intent(in)    :: job, side
            integer,   intent(in)    :: n, ilo, ihi, m, ldv
            real(ep),  intent(in)    :: lscale(*), rscale(*)
            real(ep),  intent(inout) :: V(ldv,*)
            integer,   intent(out)   :: info
        end subroutine dggbak_quad

        subroutine zggbak_quad(job, side, n, ilo, ihi, lscale, rscale, m, V, ldv, info)
            import :: ep
            character,   intent(in)    :: job, side
            integer,     intent(in)    :: n, ilo, ihi, m, ldv
            real(ep),    intent(in)    :: lscale(*), rscale(*)
            complex(ep), intent(inout) :: V(ldv,*)
            integer,     intent(out)   :: info
        end subroutine zggbak_quad

        subroutine dgghrd_quad(compq, compz, n, ilo, ihi, A, lda, B, ldb, &
                          Q, ldq, Z, ldz, info)
            import :: ep
            character, intent(in)    :: compq, compz
            integer,   intent(in)    :: n, ilo, ihi, lda, ldb, ldq, ldz
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
            integer,   intent(out)   :: info
        end subroutine dgghrd_quad

        subroutine zgghrd_quad(compq, compz, n, ilo, ihi, A, lda, B, ldb, &
                          Q, ldq, Z, ldz, info)
            import :: ep
            character,   intent(in)    :: compq, compz
            integer,     intent(in)    :: n, ilo, ihi, lda, ldb, ldq, ldz
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
            integer,     intent(out)   :: info
        end subroutine zgghrd_quad

        ! ── Generalized QR / RQ factorizations ───────────────────────
        subroutine dggqrf_quad(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: n, m, p, lda, ldb, lwork
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dggqrf_quad

        subroutine zggqrf_quad(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, m, p, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zggqrf_quad

        subroutine dggrqf_quad(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, p, n, lda, ldb, lwork
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dggrqf_quad

        subroutine zggrqf_quad(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, p, n, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: taua(*), taub(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zggrqf_quad

        ! ── Sylvester equation ───────────────────────────────────────
        subroutine dtrsyl_quad(trana, tranb, isgn, m, n, A, lda, B, ldb, C, ldc, &
                          scale, info)
            import :: ep
            character, intent(in)    :: trana, tranb
            integer,   intent(in)    :: isgn, m, n, lda, ldb, ldc
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: scale
            integer,   intent(out)   :: info
        end subroutine dtrsyl_quad

        subroutine ztrsyl_quad(trana, tranb, isgn, m, n, A, lda, B, ldb, C, ldc, &
                          scale, info)
            import :: ep
            character,   intent(in)    :: trana, tranb
            integer,     intent(in)    :: isgn, m, n, lda, ldb, ldc
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: C(ldc,*)
            real(ep),    intent(out)   :: scale
            integer,     intent(out)   :: info
        end subroutine ztrsyl_quad

        ! ── Hessenberg Schur factorization ───────────────────────────
        subroutine dhseqr_quad(job, compz, n, ilo, ihi, H, ldh, WR, WI, Z, ldz, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: job, compz
            integer,   intent(in)    :: n, ilo, ihi, ldh, ldz, lwork
            real(ep),  intent(inout) :: H(ldh,*), Z(ldz,*)
            real(ep),  intent(out)   :: WR(*), WI(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dhseqr_quad

        subroutine zhseqr_quad(job, compz, n, ilo, ihi, H, ldh, W, Z, ldz, &
                          work, lwork, info)
            import :: ep
            character,   intent(in)    :: job, compz
            integer,     intent(in)    :: n, ilo, ihi, ldh, ldz, lwork
            complex(ep), intent(inout) :: H(ldh,*), Z(ldz,*)
            complex(ep), intent(out)   :: W(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhseqr_quad

        ! ── Eigenvectors of (quasi-)triangular Schur form ────────────
        subroutine dtrevc_quad(side, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, &
                          mm, m, work, info)
            import :: ep
            character, intent(in)    :: side, howmny
            integer,   intent(in)    :: n, ldt, ldvl, ldvr, mm
            logical,   intent(inout) :: sel(*)
            real(ep),  intent(in)    :: T(ldt,*)
            real(ep),  intent(inout) :: VL(ldvl,*), VR(ldvr,*)
            integer,   intent(out)   :: m
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dtrevc_quad

        subroutine ztrevc_quad(side, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, &
                          mm, m, work, rwork, info)
            import :: ep
            character,   intent(in)    :: side, howmny
            integer,     intent(in)    :: n, ldt, ldvl, ldvr, mm
            logical,     intent(in)    :: sel(*)
            complex(ep), intent(inout) :: T(ldt,*), VL(ldvl,*), VR(ldvr,*)
            integer,     intent(out)   :: m
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine ztrevc_quad

        ! ── Generalized non-symmetric eigenvalue ─────────────────────
        subroutine dggev_quad(jobvl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, &
                         VL, ldvl, VR, ldvr, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobvl, jobvr
            integer,   intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
            real(ep),  intent(out)   :: VL(ldvl,*), VR(ldvr,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dggev_quad

        subroutine zggev_quad(jobvl, jobvr, n, A, lda, B, ldb, alpha, beta, &
                         VL, ldvl, VR, ldvr, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobvl, jobvr
            integer,     intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: alpha(*), beta(*)
            complex(ep), intent(out)   :: VL(ldvl,*), VR(ldvr,*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zggev_quad

        ! ── Modern (LAPACK 3.7+) blocked QR/LQ ───────────────────────
        subroutine dgeqr_quad(m, n, A, lda, T, tsize, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, tsize, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: T(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqr_quad

        subroutine zgeqr_quad(m, n, A, lda, T, tsize, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, tsize, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: T(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqr_quad

        subroutine dgelq_quad(m, n, A, lda, T, tsize, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, tsize, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: T(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgelq_quad

        subroutine zgelq_quad(m, n, A, lda, T, tsize, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, tsize, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: T(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgelq_quad

        subroutine dgemqr_quad(side, trans, m, n, k, A, lda, T, tsize, &
                          C, ldc, work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, tsize, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), T(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgemqr_quad

        subroutine zgemqr_quad(side, trans, m, n, k, A, lda, T, tsize, &
                          C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, tsize, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), T(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgemqr_quad

        subroutine dgemlq_quad(side, trans, m, n, k, A, lda, T, tsize, &
                          C, ldc, work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, tsize, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), T(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgemlq_quad

        subroutine zgemlq_quad(side, trans, m, n, k, A, lda, T, tsize, &
                          C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, tsize, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), T(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgemlq_quad

        ! ── Recursive blocked QR/LQ with explicit T (LAPACK 3.6+) ────
        subroutine dgeqrt_quad(m, n, nb, A, lda, T, ldt, work, info)
            import :: ep
            integer,  intent(in)    :: m, n, nb, lda, ldt
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: T(ldt,*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqrt_quad

        subroutine zgeqrt_quad(m, n, nb, A, lda, T, ldt, work, info)
            import :: ep
            integer,     intent(in)    :: m, n, nb, lda, ldt
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: T(ldt,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqrt_quad

        subroutine dgelqt_quad(m, n, mb, A, lda, T, ldt, work, info)
            import :: ep
            integer,  intent(in)    :: m, n, mb, lda, ldt
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: T(ldt,*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgelqt_quad

        subroutine zgelqt_quad(m, n, mb, A, lda, T, ldt, work, info)
            import :: ep
            integer,     intent(in)    :: m, n, mb, lda, ldt
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: T(ldt,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgelqt_quad

        subroutine dgemqrt_quad(side, trans, m, n, k, nb, V, ldv, T, ldt, &
                           C, ldc, work, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, nb, ldv, ldt, ldc
            real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgemqrt_quad

        subroutine zgemqrt_quad(side, trans, m, n, k, nb, V, ldv, T, ldt, &
                           C, ldc, work, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, nb, ldv, ldt, ldc
            complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgemqrt_quad

        subroutine dgemlqt_quad(side, trans, m, n, k, mb, V, ldv, T, ldt, &
                           C, ldc, work, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, mb, ldv, ldt, ldc
            real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgemlqt_quad

        subroutine zgemlqt_quad(side, trans, m, n, k, mb, V, ldv, T, ldt, &
                           C, ldc, work, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, mb, ldv, ldt, ldc
            complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgemlqt_quad

        ! ── Tall-skinny least squares ────────────────────────────────
        subroutine dgetsls_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgetsls_quad

        subroutine zgetsls_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgetsls_quad

        ! ── Expert linear solvers (SVX) ──────────────────────────────
        subroutine dgesvx_quad(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, &
                          R, C, B, ldb, X, ldx, rcond, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: fact, trans
            character, intent(inout) :: equed
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            real(ep),  intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*), R(*), C(*)
            integer,   intent(inout) :: ipiv(*)
            real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgesvx_quad

        subroutine zgesvx_quad(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, &
                          R, C, B, ldb, X, ldx, rcond, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, trans
            character,   intent(inout) :: equed
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),    intent(inout) :: R(*), C(*)
            integer,     intent(inout) :: ipiv(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zgesvx_quad

        subroutine dposvx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S, &
                          B, ldb, X, ldx, rcond, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: fact, uplo
            character, intent(inout) :: equed
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            real(ep),  intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*), S(*)
            real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dposvx_quad

        subroutine zposvx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S, &
                          B, ldb, X, ldx, rcond, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, uplo
            character,   intent(inout) :: equed
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
            complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),    intent(inout) :: S(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zposvx_quad

        ! ── Banded/packed positive-definite expert solvers ───────────
        subroutine dpbsvx_quad(fact, uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, equed, &
                          S, B, ldb, X, ldx, rcond, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: fact, uplo
            character, intent(inout) :: equed
            integer,   intent(in)    :: n, kd, nrhs, ldab, ldafb, ldb, ldx
            real(ep),  intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*), S(*)
            real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dpbsvx_quad

        subroutine zpbsvx_quad(fact, uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, equed, &
                          S, B, ldb, X, ldx, rcond, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, uplo
            character,   intent(inout) :: equed
            integer,     intent(in)    :: n, kd, nrhs, ldab, ldafb, ldb, ldx
            complex(ep), intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            real(ep),    intent(inout) :: S(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zpbsvx_quad

        subroutine dppsvx_quad(fact, uplo, n, nrhs, AP, AFP, equed, S, B, ldb, &
                          X, ldx, rcond, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: fact, uplo
            character, intent(inout) :: equed
            integer,   intent(in)    :: n, nrhs, ldb, ldx
            real(ep),  intent(inout) :: AP(*), AFP(*), B(ldb,*), S(*)
            real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dppsvx_quad

        subroutine zppsvx_quad(fact, uplo, n, nrhs, AP, AFP, equed, S, B, ldb, &
                          X, ldx, rcond, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, uplo
            character,   intent(inout) :: equed
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            complex(ep), intent(inout) :: AP(*), AFP(*), B(ldb,*)
            real(ep),    intent(inout) :: S(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zppsvx_quad

        ! ── Banded general / symmetric / Hermitian expert solvers ────
        subroutine dgbsvx_quad(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, &
                          equed, R, C, B, ldb, X, ldx, rcond, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: fact, trans
            character, intent(inout) :: equed
            integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx
            real(ep),  intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*), R(*), C(*)
            integer,   intent(inout) :: ipiv(*)
            real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgbsvx_quad

        subroutine zgbsvx_quad(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, &
                          equed, R, C, B, ldb, X, ldx, rcond, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, trans
            character,   intent(inout) :: equed
            integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx
            complex(ep), intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            real(ep),    intent(inout) :: R(*), C(*)
            integer,     intent(inout) :: ipiv(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zgbsvx_quad

        subroutine dsysvx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, rcond, ferr, berr, work, lwork, iwork, info)
            import :: ep
            character, intent(in)    :: fact, uplo
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, lwork
            real(ep),  intent(inout) :: AF(ldaf,*)
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            integer,   intent(inout) :: ipiv(*)
            real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsysvx_quad

        subroutine zsysvx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, rcond, ferr, berr, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, uplo
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, lwork
            complex(ep), intent(inout) :: AF(ldaf,*)
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            integer,     intent(inout) :: ipiv(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zsysvx_quad

        subroutine zhesvx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, &
                          X, ldx, rcond, ferr, berr, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, uplo
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, lwork
            complex(ep), intent(inout) :: AF(ldaf,*)
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            integer,     intent(inout) :: ipiv(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zhesvx_quad

        ! ── Phase L1 — 2-stage symmetric/Hermitian eigensolvers ──────
        subroutine dsyev_2stage_quad(jobz, uplo, n, A, lda, w, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsyev_2stage_quad

        subroutine zheev_2stage_quad(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zheev_2stage_quad

        subroutine dsyevd_2stage_quad(jobz, uplo, n, A, lda, w, work, lwork, &
                                 iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, lda, lwork, liwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsyevd_2stage_quad

        subroutine zheevd_2stage_quad(jobz, uplo, n, A, lda, w, work, lwork, &
                                 rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, lda, lwork, lrwork, liwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zheevd_2stage_quad

        subroutine dsyevr_2stage_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                                 abstol, m, w, z, ldz, isuppz, &
                                 work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, lda, il, iu, ldz, lwork, liwork
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: m, isuppz(*), iwork(*), info
        end subroutine dsyevr_2stage_quad

        subroutine zheevr_2stage_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                                 abstol, m, w, z, ldz, isuppz, &
                                 work, lwork, rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, lda, il, iu, ldz, lwork, lrwork, liwork
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            integer,     intent(out)   :: m, isuppz(*), iwork(*), info
        end subroutine zheevr_2stage_quad

        subroutine dsyevx_2stage_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                                 abstol, m, w, z, ldz, work, lwork, iwork, &
                                 ifail, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, lda, il, iu, ldz, lwork
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine dsyevx_2stage_quad

        subroutine zheevx_2stage_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, &
                                 abstol, m, w, z, ldz, work, lwork, rwork, &
                                 iwork, ifail, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, lda, il, iu, ldz, lwork
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            integer,     intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine zheevx_2stage_quad

        subroutine dsbev_2stage_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, &
                                work, lwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, kd, ldab, ldz, lwork
            real(ep),  intent(inout) :: AB(ldab,*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsbev_2stage_quad

        subroutine zhbev_2stage_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, &
                                work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, kd, ldab, ldz, lwork
            complex(ep), intent(inout) :: AB(ldab,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhbev_2stage_quad

        subroutine dsbevd_2stage_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, &
                                 work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: n, kd, ldab, ldz, lwork, liwork
            real(ep),  intent(inout) :: AB(ldab,*)
            real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsbevd_2stage_quad

        subroutine zhbevd_2stage_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, &
                                 work, lwork, rwork, lrwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: n, kd, ldab, ldz, lwork, lrwork, liwork
            complex(ep), intent(inout) :: AB(ldab,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: z(ldz,*), work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zhbevd_2stage_quad

        subroutine dsbevx_2stage_quad(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, &
                                 vl, vu, il, iu, abstol, m, w, z, ldz, &
                                 work, lwork, iwork, ifail, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, kd, ldab, ldq, il, iu, ldz, lwork
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: AB(ldab,*)
            real(ep),  intent(out)   :: Q(ldq,*), w(*), z(ldz,*), work(*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine dsbevx_2stage_quad

        subroutine zhbevx_2stage_quad(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, &
                                 vl, vu, il, iu, abstol, m, w, z, ldz, &
                                 work, lwork, rwork, iwork, ifail, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, kd, ldab, ldq, il, iu, ldz, lwork
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: AB(ldab,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: Q(ldq,*), z(ldz,*), work(*)
            integer,     intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine zhbevx_2stage_quad

        subroutine dsygv_2stage_quad(itype, jobz, uplo, n, A, lda, B, ldb, w, &
                                work, lwork, info)
            import :: ep
            character, intent(in)    :: jobz, uplo
            integer,   intent(in)    :: itype, n, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: w(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsygv_2stage_quad

        subroutine zhegv_2stage_quad(itype, jobz, uplo, n, A, lda, B, ldb, w, &
                                work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobz, uplo
            integer,     intent(in)    :: itype, n, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhegv_2stage_quad

        ! ── Phase L2 — Hessenberg/Schur expert drivers ───────────────
        subroutine dgees_quad(jobvs, sort, select, n, A, lda, sdim, wr, wi, &
                         vs, ldvs, work, lwork, bwork, info)
            import :: ep
            character, intent(in)    :: jobvs, sort
            interface
                logical function select(re, im)
                    import :: ep
                    real(ep), intent(in) :: re, im
                end function select
            end interface
            integer,   intent(in)    :: n, lda, ldvs, lwork
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: sdim, info
            real(ep),  intent(out)   :: wr(*), wi(*), vs(ldvs,*), work(*)
            logical,   intent(out)   :: bwork(*)
        end subroutine dgees_quad

        subroutine zgees_quad(jobvs, sort, select, n, A, lda, sdim, w, &
                         vs, ldvs, work, lwork, rwork, bwork, info)
            import :: ep
            character,   intent(in)    :: jobvs, sort
            interface
                logical function select(z)
                    import :: ep
                    complex(ep), intent(in) :: z
                end function select
            end interface
            integer,     intent(in)    :: n, lda, ldvs, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: sdim, info
            complex(ep), intent(out)   :: w(*), vs(ldvs,*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            logical,     intent(out)   :: bwork(*)
        end subroutine zgees_quad

        subroutine dgeesx_quad(jobvs, sort, select, sense, n, A, lda, sdim, &
                          wr, wi, vs, ldvs, rconde, rcondv, work, lwork, &
                          iwork, liwork, bwork, info)
            import :: ep
            character, intent(in)    :: jobvs, sort, sense
            interface
                logical function select(re, im)
                    import :: ep
                    real(ep), intent(in) :: re, im
                end function select
            end interface
            integer,   intent(in)    :: n, lda, ldvs, lwork, liwork
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: sdim, iwork(*), info
            real(ep),  intent(out)   :: wr(*), wi(*), vs(ldvs,*), work(*)
            real(ep),  intent(out)   :: rconde, rcondv
            logical,   intent(out)   :: bwork(*)
        end subroutine dgeesx_quad

        subroutine zgeesx_quad(jobvs, sort, select, sense, n, A, lda, sdim, &
                          w, vs, ldvs, rconde, rcondv, work, lwork, &
                          rwork, bwork, info)
            import :: ep
            character,   intent(in)    :: jobvs, sort, sense
            interface
                logical function select(z)
                    import :: ep
                    complex(ep), intent(in) :: z
                end function select
            end interface
            integer,     intent(in)    :: n, lda, ldvs, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: sdim, info
            complex(ep), intent(out)   :: w(*), vs(ldvs,*), work(*)
            real(ep),    intent(out)   :: rconde, rcondv, rwork(*)
            logical,     intent(out)   :: bwork(*)
        end subroutine zgeesx_quad

        subroutine dgeevx_quad(balanc, jobvl, jobvr, sense, n, A, lda, &
                          wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, &
                          scale, abnrm, rconde, rcondv, work, lwork, &
                          iwork, info)
            import :: ep
            character, intent(in)    :: balanc, jobvl, jobvr, sense
            integer,   intent(in)    :: n, lda, ldvl, ldvr, lwork
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: ilo, ihi, iwork(*), info
            real(ep),  intent(out)   :: wr(*), wi(*), vl(ldvl,*), vr(ldvr,*)
            real(ep),  intent(out)   :: scale(*), abnrm
            real(ep),  intent(out)   :: rconde(*), rcondv(*), work(*)
        end subroutine dgeevx_quad

        subroutine zgeevx_quad(balanc, jobvl, jobvr, sense, n, A, lda, w, &
                          vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, &
                          rconde, rcondv, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: balanc, jobvl, jobvr, sense
            integer,     intent(in)    :: n, lda, ldvl, ldvr, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ilo, ihi, info
            complex(ep), intent(out)   :: w(*), vl(ldvl,*), vr(ldvr,*), work(*)
            real(ep),    intent(out)   :: scale(*), abnrm
            real(ep),    intent(out)   :: rconde(*), rcondv(*), rwork(*)
        end subroutine zgeevx_quad

        subroutine dhsein_quad(side, eigsrc, initv, select, n, H, ldh, &
                          wr, wi, vl, ldvl, vr, ldvr, mm, m, work, &
                          ifaill, ifailr, info)
            import :: ep
            character, intent(in)    :: side, eigsrc, initv
            integer,   intent(in)    :: n, ldh, ldvl, ldvr, mm
            logical,   intent(inout) :: select(*)
            real(ep),  intent(in)    :: H(ldh,*)
            real(ep),  intent(inout) :: wr(*)
            real(ep),  intent(in)    :: wi(*)
            real(ep),  intent(inout) :: vl(ldvl,*), vr(ldvr,*)
            integer,   intent(out)   :: m, ifaill(*), ifailr(*), info
            real(ep),  intent(out)   :: work(*)
        end subroutine dhsein_quad

        subroutine zhsein_quad(side, eigsrc, initv, select, n, H, ldh, &
                          w, vl, ldvl, vr, ldvr, mm, m, work, rwork, &
                          ifaill, ifailr, info)
            import :: ep
            character,   intent(in)    :: side, eigsrc, initv
            integer,     intent(in)    :: n, ldh, ldvl, ldvr, mm
            logical,     intent(in)    :: select(*)
            complex(ep), intent(in)    :: H(ldh,*)
            complex(ep), intent(inout) :: w(*), vl(ldvl,*), vr(ldvr,*)
            integer,     intent(out)   :: m, ifaill(*), ifailr(*), info
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rwork(*)
        end subroutine zhsein_quad

        subroutine dtrevc3_quad(side, howmny, select, n, T, ldt, vl, ldvl, &
                           vr, ldvr, mm, m, work, lwork, info)
            import :: ep
            character, intent(in)    :: side, howmny
            integer,   intent(in)    :: n, ldt, ldvl, ldvr, mm, lwork
            logical,   intent(inout) :: select(*)
            real(ep),  intent(in)    :: T(ldt,*)
            real(ep),  intent(inout) :: vl(ldvl,*), vr(ldvr,*)
            integer,   intent(out)   :: m, info
            real(ep),  intent(out)   :: work(*)
        end subroutine dtrevc3_quad

        subroutine ztrevc3_quad(side, howmny, select, n, T, ldt, vl, ldvl, &
                           vr, ldvr, mm, m, work, lwork, rwork, lrwork, info)
            import :: ep
            character,   intent(in)    :: side, howmny
            integer,     intent(in)    :: n, ldt, ldvl, ldvr, mm, lwork, lrwork
            logical,     intent(in)    :: select(*)
            complex(ep), intent(inout) :: T(ldt,*), vl(ldvl,*), vr(ldvr,*)
            integer,     intent(out)   :: m, info
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rwork(*)
        end subroutine ztrevc3_quad

        subroutine dtrexc_quad(compq, n, T, ldt, Q, ldq, ifst, ilst, work, info)
            import :: ep
            character, intent(in)    :: compq
            integer,   intent(in)    :: n, ldt, ldq
            integer,   intent(inout) :: ifst, ilst
            real(ep),  intent(inout) :: T(ldt,*), Q(ldq,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dtrexc_quad

        subroutine ztrexc_quad(compq, n, T, ldt, Q, ldq, ifst, ilst, info)
            import :: ep
            character,   intent(in)    :: compq
            integer,     intent(in)    :: n, ldt, ldq, ifst, ilst
            complex(ep), intent(inout) :: T(ldt,*), Q(ldq,*)
            integer,     intent(out)   :: info
        end subroutine ztrexc_quad

        subroutine dtrsen_quad(job, compq, select, n, T, ldt, Q, ldq, &
                          wr, wi, m, s, sep, work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: job, compq
            integer,   intent(in)    :: n, ldt, ldq, lwork, liwork
            logical,   intent(in)    :: select(*)
            real(ep),  intent(inout) :: T(ldt,*), Q(ldq,*)
            real(ep),  intent(out)   :: wr(*), wi(*), s, sep, work(*)
            integer,   intent(out)   :: m, iwork(*), info
        end subroutine dtrsen_quad

        subroutine ztrsen_quad(job, compq, select, n, T, ldt, Q, ldq, &
                          w, m, s, sep, work, lwork, info)
            import :: ep
            character,   intent(in)    :: job, compq
            integer,     intent(in)    :: n, ldt, ldq, lwork
            logical,     intent(in)    :: select(*)
            complex(ep), intent(inout) :: T(ldt,*), Q(ldq,*)
            complex(ep), intent(out)   :: w(*), work(*)
            real(ep),    intent(out)   :: s, sep
            integer,     intent(out)   :: m, info
        end subroutine ztrsen_quad

        ! ── Phase L3 — Generalized non-symmetric expert drivers + QZ ──
        subroutine dgges_quad(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, &
                         sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, &
                         work, lwork, bwork, info)
            import :: ep
            character, intent(in) :: jobvsl, jobvsr, sort
            interface
                logical function selctg(re, im, b)
                    import :: ep
                    real(ep), intent(in) :: re, im, b
                end function selctg
            end interface
            integer,   intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: sdim, info
            real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
            real(ep),  intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
            logical,   intent(out)   :: bwork(*)
        end subroutine dgges_quad

        subroutine zgges_quad(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, &
                         sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, &
                         work, lwork, rwork, bwork, info)
            import :: ep
            character,   intent(in) :: jobvsl, jobvsr, sort
            interface
                logical function selctg(a, b)
                    import :: ep
                    complex(ep), intent(in) :: a, b
                end function selctg
            end interface
            integer,     intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: sdim, info
            complex(ep), intent(out)   :: alpha(*), beta(*)
            complex(ep), intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            logical,     intent(out)   :: bwork(*)
        end subroutine zgges_quad

        subroutine dgges3_quad(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, &
                          sdim, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, &
                          work, lwork, bwork, info)
            import :: ep
            character, intent(in) :: jobvsl, jobvsr, sort
            interface
                logical function selctg(re, im, b)
                    import :: ep
                    real(ep), intent(in) :: re, im, b
                end function selctg
            end interface
            integer,   intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: sdim, info
            real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
            real(ep),  intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
            logical,   intent(out)   :: bwork(*)
        end subroutine dgges3_quad

        subroutine zgges3_quad(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, &
                          sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, &
                          work, lwork, rwork, bwork, info)
            import :: ep
            character,   intent(in) :: jobvsl, jobvsr, sort
            interface
                logical function selctg(a, b)
                    import :: ep
                    complex(ep), intent(in) :: a, b
                end function selctg
            end interface
            integer,     intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: sdim, info
            complex(ep), intent(out)   :: alpha(*), beta(*)
            complex(ep), intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            logical,     intent(out)   :: bwork(*)
        end subroutine zgges3_quad

        subroutine dggesx_quad(jobvsl, jobvsr, sort, selctg, sense, n, A, lda, &
                          B, ldb, sdim, alphar, alphai, beta, vsl, ldvsl, &
                          vsr, ldvsr, rconde, rcondv, work, lwork, &
                          iwork, liwork, bwork, info)
            import :: ep
            character, intent(in) :: jobvsl, jobvsr, sort, sense
            interface
                logical function selctg(re, im, b)
                    import :: ep
                    real(ep), intent(in) :: re, im, b
                end function selctg
            end interface
            integer,   intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork, liwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: sdim, iwork(*), info
            real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
            real(ep),  intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
            real(ep),  intent(out)   :: rconde(2), rcondv(2)
            logical,   intent(out)   :: bwork(*)
        end subroutine dggesx_quad

        subroutine zggesx_quad(jobvsl, jobvsr, sort, selctg, sense, n, A, lda, &
                          B, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, &
                          rconde, rcondv, work, lwork, rwork, iwork, liwork, &
                          bwork, info)
            import :: ep
            character,   intent(in) :: jobvsl, jobvsr, sort, sense
            interface
                logical function selctg(a, b)
                    import :: ep
                    complex(ep), intent(in) :: a, b
                end function selctg
            end interface
            integer,     intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork, liwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: sdim, iwork(*), info
            complex(ep), intent(out)   :: alpha(*), beta(*)
            complex(ep), intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
            real(ep),    intent(out)   :: rconde(2), rcondv(2), rwork(*)
            logical,     intent(out)   :: bwork(*)
        end subroutine zggesx_quad

        subroutine dggev3_quad(jobvl, jobvr, n, A, lda, B, ldb, alphar, alphai, &
                          beta, vl, ldvl, vr, ldvr, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobvl, jobvr
            integer,   intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: info
            real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
            real(ep),  intent(out)   :: vl(ldvl,*), vr(ldvr,*), work(*)
        end subroutine dggev3_quad

        subroutine zggev3_quad(jobvl, jobvr, n, A, lda, B, ldb, alpha, beta, &
                          vl, ldvl, vr, ldvr, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: jobvl, jobvr
            integer,     intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: info
            complex(ep), intent(out)   :: alpha(*), beta(*)
            complex(ep), intent(out)   :: vl(ldvl,*), vr(ldvr,*), work(*)
            real(ep),    intent(out)   :: rwork(*)
        end subroutine zggev3_quad

        subroutine dggevx_quad(balanc, jobvl, jobvr, sense, n, A, lda, B, ldb, &
                          alphar, alphai, beta, vl, ldvl, vr, ldvr, &
                          ilo, ihi, lscale, rscale, abnrm, bbnrm, &
                          rconde, rcondv, work, lwork, iwork, bwork, info)
            import :: ep
            character, intent(in)    :: balanc, jobvl, jobvr, sense
            integer,   intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: ilo, ihi, iwork(*), info
            real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
            real(ep),  intent(out)   :: vl(ldvl,*), vr(ldvr,*)
            real(ep),  intent(out)   :: lscale(*), rscale(*), abnrm, bbnrm
            real(ep),  intent(out)   :: rconde(*), rcondv(*), work(*)
            logical,   intent(out)   :: bwork(*)
        end subroutine dggevx_quad

        subroutine zggevx_quad(balanc, jobvl, jobvr, sense, n, A, lda, B, ldb, &
                          alpha, beta, vl, ldvl, vr, ldvr, &
                          ilo, ihi, lscale, rscale, abnrm, bbnrm, &
                          rconde, rcondv, work, lwork, rwork, iwork, bwork, info)
            import :: ep
            character,   intent(in)    :: balanc, jobvl, jobvr, sense
            integer,     intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ilo, ihi, iwork(*), info
            complex(ep), intent(out)   :: alpha(*), beta(*)
            complex(ep), intent(out)   :: vl(ldvl,*), vr(ldvr,*), work(*)
            real(ep),    intent(out)   :: lscale(*), rscale(*), abnrm, bbnrm
            real(ep),    intent(out)   :: rconde(*), rcondv(*), rwork(*)
            logical,     intent(out)   :: bwork(*)
        end subroutine zggevx_quad

        subroutine dhgeqz_quad(job, compq, compz, n, ilo, ihi, H, ldh, T, ldt, &
                          alphar, alphai, beta, Q, ldq, Z, ldz, work, lwork, info)
            import :: ep
            character, intent(in)    :: job, compq, compz
            integer,   intent(in)    :: n, ilo, ihi, ldh, ldt, ldq, ldz, lwork
            real(ep),  intent(inout) :: H(ldh,*), T(ldt,*), Q(ldq,*), Z(ldz,*)
            real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dhgeqz_quad

        subroutine zhgeqz_quad(job, compq, compz, n, ilo, ihi, H, ldh, T, ldt, &
                          alpha, beta, Q, ldq, Z, ldz, work, lwork, rwork, info)
            import :: ep
            character,   intent(in)    :: job, compq, compz
            integer,     intent(in)    :: n, ilo, ihi, ldh, ldt, ldq, ldz, lwork
            complex(ep), intent(inout) :: H(ldh,*), T(ldt,*), Q(ldq,*), Z(ldz,*)
            complex(ep), intent(out)   :: alpha(*), beta(*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zhgeqz_quad

        subroutine dtgevc_quad(side, howmny, select, n, S, lds, P, ldp, &
                          vl, ldvl, vr, ldvr, mm, m, work, info)
            import :: ep
            character, intent(in)    :: side, howmny
            integer,   intent(in)    :: n, lds, ldp, ldvl, ldvr, mm
            logical,   intent(in)    :: select(*)
            real(ep),  intent(in)    :: S(lds,*), P(ldp,*)
            real(ep),  intent(inout) :: vl(ldvl,*), vr(ldvr,*)
            integer,   intent(out)   :: m, info
            real(ep),  intent(out)   :: work(*)
        end subroutine dtgevc_quad

        subroutine ztgevc_quad(side, howmny, select, n, S, lds, P, ldp, &
                          vl, ldvl, vr, ldvr, mm, m, work, rwork, info)
            import :: ep
            character,   intent(in)    :: side, howmny
            integer,     intent(in)    :: n, lds, ldp, ldvl, ldvr, mm
            logical,     intent(in)    :: select(*)
            complex(ep), intent(in)    :: S(lds,*), P(ldp,*)
            complex(ep), intent(inout) :: vl(ldvl,*), vr(ldvr,*)
            integer,     intent(out)   :: m, info
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rwork(*)
        end subroutine ztgevc_quad

        ! ── Phase L4 — Generalized post-processing + constrained LS ──
        subroutine dtgexc_quad(wantq, wantz, n, A, lda, B, ldb, Q, ldq, Z, ldz, &
                          ifst, ilst, work, lwork, info)
            import :: ep
            logical,   intent(in)    :: wantq, wantz
            integer,   intent(in)    :: n, lda, ldb, ldq, ldz, lwork
            integer,   intent(inout) :: ifst, ilst
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dtgexc_quad

        subroutine ztgexc_quad(wantq, wantz, n, A, lda, B, ldb, Q, ldq, Z, ldz, &
                          ifst, ilst, info)
            import :: ep
            logical,     intent(in)    :: wantq, wantz
            integer,     intent(in)    :: n, lda, ldb, ldq, ldz, ifst, ilst
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
            integer,     intent(out)   :: info
        end subroutine ztgexc_quad

        subroutine dtgsen_quad(ijob, wantq, wantz, select, n, A, lda, B, ldb, &
                          alphar, alphai, beta, Q, ldq, Z, ldz, m, pl, pr, &
                          dif, work, lwork, iwork, liwork, info)
            import :: ep
            integer,   intent(in)    :: ijob, n, lda, ldb, ldq, ldz, lwork, liwork
            logical,   intent(in)    :: wantq, wantz, select(*)
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
            real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*), pl, pr
            real(ep),  intent(out)   :: dif(2), work(*)
            integer,   intent(out)   :: m, iwork(*), info
        end subroutine dtgsen_quad

        subroutine ztgsen_quad(ijob, wantq, wantz, select, n, A, lda, B, ldb, &
                          alpha, beta, Q, ldq, Z, ldz, m, pl, pr, dif, &
                          work, lwork, iwork, liwork, info)
            import :: ep
            integer,     intent(in)    :: ijob, n, lda, ldb, ldq, ldz, lwork, liwork
            logical,     intent(in)    :: wantq, wantz, select(*)
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
            complex(ep), intent(out)   :: alpha(*), beta(*), work(*)
            real(ep),    intent(out)   :: pl, pr, dif(2)
            integer,     intent(out)   :: m, iwork(*), info
        end subroutine ztgsen_quad

        subroutine dtgsna_quad(job, howmny, select, n, A, lda, B, ldb, vl, ldvl, &
                          vr, ldvr, s, dif, mm, m, work, lwork, iwork, info)
            import :: ep
            character, intent(in)    :: job, howmny
            integer,   intent(in)    :: n, lda, ldb, ldvl, ldvr, mm, lwork
            logical,   intent(in)    :: select(*)
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*), vl(ldvl,*), vr(ldvr,*)
            real(ep),  intent(out)   :: s(*), dif(*), work(*)
            integer,   intent(out)   :: m, iwork(*), info
        end subroutine dtgsna_quad

        subroutine ztgsna_quad(job, howmny, select, n, A, lda, B, ldb, vl, ldvl, &
                          vr, ldvr, s, dif, mm, m, work, lwork, iwork, info)
            import :: ep
            character,   intent(in)    :: job, howmny
            integer,     intent(in)    :: n, lda, ldb, ldvl, ldvr, mm, lwork
            logical,     intent(in)    :: select(*)
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*), vl(ldvl,*), vr(ldvr,*)
            real(ep),    intent(out)   :: s(*), dif(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: m, iwork(*), info
        end subroutine ztgsna_quad

        subroutine dtgsja_quad(jobu, jobv, jobq, m, p, n, k, l, A, lda, B, ldb, &
                          tola, tolb, alpha, beta, U, ldu, V, ldv, Q, ldq, &
                          work, ncycle, info)
            import :: ep
            character, intent(in)    :: jobu, jobv, jobq
            integer,   intent(in)    :: m, p, n, k, l, lda, ldb, ldu, ldv, ldq
            real(ep),  intent(in)    :: tola, tolb
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: U(ldu,*), V(ldv,*), Q(ldq,*)
            real(ep),  intent(out)   :: alpha(*), beta(*), work(*)
            integer,   intent(out)   :: ncycle, info
        end subroutine dtgsja_quad

        subroutine ztgsja_quad(jobu, jobv, jobq, m, p, n, k, l, A, lda, B, ldb, &
                          tola, tolb, alpha, beta, U, ldu, V, ldv, Q, ldq, &
                          work, ncycle, info)
            import :: ep
            character,   intent(in)    :: jobu, jobv, jobq
            integer,     intent(in)    :: m, p, n, k, l, lda, ldb, ldu, ldv, ldq
            real(ep),    intent(in)    :: tola, tolb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: U(ldu,*), V(ldv,*), Q(ldq,*)
            real(ep),    intent(out)   :: alpha(*), beta(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: ncycle, info
        end subroutine ztgsja_quad

        subroutine dtgsyl_quad(trans, ijob, m, n, A, lda, B, ldb, C, ldc, &
                          D, ldd, E, lde, F, ldf, scale, dif, work, lwork, &
                          iwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: ijob, m, n, lda, ldb, ldc, ldd, lde, ldf, lwork
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*), D(ldd,*), E(lde,*)
            real(ep),  intent(inout) :: C(ldc,*), F(ldf,*)
            real(ep),  intent(out)   :: scale, dif, work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dtgsyl_quad

        subroutine ztgsyl_quad(trans, ijob, m, n, A, lda, B, ldb, C, ldc, &
                          D, ldd, E, lde, F, ldf, scale, dif, work, lwork, &
                          iwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: ijob, m, n, lda, ldb, ldc, ldd, lde, ldf, lwork
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*), D(ldd,*), E(lde,*)
            complex(ep), intent(inout) :: C(ldc,*), F(ldf,*)
            real(ep),    intent(out)   :: scale, dif
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine ztgsyl_quad

        subroutine dggglm_quad(n, m, p, A, lda, B, ldb, d, x, y, work, lwork, info)
            import :: ep
            integer,   intent(in)    :: n, m, p, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*), d(*)
            real(ep),  intent(out)   :: x(*), y(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dggglm_quad

        subroutine zggglm_quad(n, m, p, A, lda, B, ldb, d, x, y, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: n, m, p, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), d(*)
            complex(ep), intent(out)   :: x(*), y(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zggglm_quad

        subroutine dgglse_quad(m, n, p, A, lda, B, ldb, c, d, x, work, lwork, info)
            import :: ep
            integer,   intent(in)    :: m, n, p, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*), c(*), d(*)
            real(ep),  intent(out)   :: x(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dgglse_quad

        subroutine zgglse_quad(m, n, p, A, lda, B, ldb, c, d, x, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, p, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), c(*), d(*)
            complex(ep), intent(out)   :: x(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgglse_quad

        subroutine dgelsd_quad(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, &
                          work, lwork, iwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep), intent(in)    :: rcond
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep), intent(out)   :: S(*), work(*)
            integer,  intent(out)   :: rank, iwork(*), info
        end subroutine dgelsd_quad

        subroutine zgelsd_quad(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, &
                          work, lwork, rwork, iwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep),    intent(in)    :: rcond
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),    intent(out)   :: S(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: rank, iwork(*), info
        end subroutine zgelsd_quad

        subroutine dgelss_quad(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, &
                          work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep), intent(in)    :: rcond
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep), intent(out)   :: S(*), work(*)
            integer,  intent(out)   :: rank, info
        end subroutine dgelss_quad

        subroutine zgelss_quad(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, &
                          work, lwork, rwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep),    intent(in)    :: rcond
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),    intent(out)   :: S(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: rank, info
        end subroutine zgelss_quad

        subroutine dgelsy_quad(m, n, nrhs, A, lda, B, ldb, jpvt, rcond, rank, &
                          work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep), intent(in)    :: rcond
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(inout) :: jpvt(*)
            integer,  intent(out)   :: rank, info
        end subroutine dgelsy_quad

        subroutine zgelsy_quad(m, n, nrhs, A, lda, B, ldb, jpvt, rcond, rank, &
                          work, lwork, rwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep),    intent(in)    :: rcond
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(inout) :: jpvt(*)
            integer,     intent(out)   :: rank, info
        end subroutine zgelsy_quad

        subroutine dgelst_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dgelst_quad

        subroutine zgelst_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgelst_quad

        subroutine dgejsv_quad(joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, &
                          sva, U, ldu, V, ldv, work, lwork, iwork, info)
            import :: ep
            character, intent(in)    :: joba, jobu, jobv, jobr, jobt, jobp
            integer,   intent(in)    :: m, n, lda, ldu, ldv, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: sva(*), U(ldu,*), V(ldv,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgejsv_quad

        subroutine zgejsv_quad(joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, &
                          sva, U, ldu, V, ldv, cwork, lwork, rwork, lrwork, &
                          iwork, info)
            import :: ep
            character,   intent(in)    :: joba, jobu, jobv, jobr, jobt, jobp
            integer,     intent(in)    :: m, n, lda, ldu, ldv, lwork, lrwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: sva(*), rwork(*)
            complex(ep), intent(out)   :: U(ldu,*), V(ldv,*), cwork(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zgejsv_quad

        subroutine dgesvj_quad(joba, jobu, jobv, m, n, A, lda, sva, mv, V, ldv, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: joba, jobu, jobv
            integer,   intent(in)    :: m, n, lda, mv, ldv, lwork
            real(ep),  intent(inout) :: A(lda,*), V(ldv,*), work(*)
            real(ep),  intent(out)   :: sva(*)
            integer,   intent(out)   :: info
        end subroutine dgesvj_quad

        subroutine zgesvj_quad(joba, jobu, jobv, m, n, A, lda, sva, mv, V, ldv, &
                          cwork, lwork, rwork, lrwork, info)
            import :: ep
            character,   intent(in)    :: joba, jobu, jobv
            integer,     intent(in)    :: m, n, lda, mv, ldv, lwork, lrwork
            complex(ep), intent(inout) :: A(lda,*), V(ldv,*), cwork(*)
            real(ep),    intent(inout) :: rwork(*)
            real(ep),    intent(out)   :: sva(*)
            integer,     intent(out)   :: info
        end subroutine zgesvj_quad

        ! Phase L18 — Auxiliary matrix norm utilities.
        function dlangb_quad(norm, n, kl, ku, AB, ldab, work) result(r)
            import :: ep
            character, intent(in) :: norm
            integer,   intent(in) :: n, kl, ku, ldab
            real(ep),  intent(in) :: AB(ldab,*)
            real(ep),  intent(out) :: work(*)
            real(ep) :: r
        end function dlangb_quad

        function zlangb_quad(norm, n, kl, ku, AB, ldab, work) result(r)
            import :: ep
            character,   intent(in) :: norm
            integer,     intent(in) :: n, kl, ku, ldab
            complex(ep), intent(in) :: AB(ldab,*)
            real(ep),    intent(out) :: work(*)
            real(ep) :: r
        end function zlangb_quad

        function dlangt_quad(norm, n, dl, d, du) result(r)
            import :: ep
            character, intent(in) :: norm
            integer,   intent(in) :: n
            real(ep),  intent(in) :: dl(*), d(*), du(*)
            real(ep) :: r
        end function dlangt_quad

        function zlangt_quad(norm, n, dl, d, du) result(r)
            import :: ep
            character,   intent(in) :: norm
            integer,     intent(in) :: n
            complex(ep), intent(in) :: dl(*), d(*), du(*)
            real(ep) :: r
        end function zlangt_quad

        function dlanhs_quad(norm, n, A, lda, work) result(r)
            import :: ep
            character, intent(in) :: norm
            integer,   intent(in) :: n, lda
            real(ep),  intent(in) :: A(lda,*)
            real(ep),  intent(out) :: work(*)
            real(ep) :: r
        end function dlanhs_quad

        function zlanhs_quad(norm, n, A, lda, work) result(r)
            import :: ep
            character,   intent(in) :: norm
            integer,     intent(in) :: n, lda
            complex(ep), intent(in) :: A(lda,*)
            real(ep),    intent(out) :: work(*)
            real(ep) :: r
        end function zlanhs_quad

        function dlansb_quad(norm, uplo, n, k, AB, ldab, work) result(r)
            import :: ep
            character, intent(in) :: norm, uplo
            integer,   intent(in) :: n, k, ldab
            real(ep),  intent(in) :: AB(ldab,*)
            real(ep),  intent(out) :: work(*)
            real(ep) :: r
        end function dlansb_quad

        function zlanhb_quad(norm, uplo, n, k, AB, ldab, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo
            integer,     intent(in) :: n, k, ldab
            complex(ep), intent(in) :: AB(ldab,*)
            real(ep),    intent(out) :: work(*)
            real(ep) :: r
        end function zlanhb_quad

        function dlansp_quad(norm, uplo, n, AP, work) result(r)
            import :: ep
            character, intent(in) :: norm, uplo
            integer,   intent(in) :: n
            real(ep),  intent(in) :: AP(*)
            real(ep),  intent(out) :: work(*)
            real(ep) :: r
        end function dlansp_quad

        function zlanhp_quad(norm, uplo, n, AP, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo
            integer,     intent(in) :: n
            complex(ep), intent(in) :: AP(*)
            real(ep),    intent(out) :: work(*)
            real(ep) :: r
        end function zlanhp_quad

        function dlanst_quad(norm, n, d, e) result(r)
            import :: ep
            character, intent(in) :: norm
            integer,   intent(in) :: n
            real(ep),  intent(in) :: d(*), e(*)
            real(ep) :: r
        end function dlanst_quad

        function zlanht_quad(norm, n, d, e) result(r)
            import :: ep
            character,   intent(in) :: norm
            integer,     intent(in) :: n
            real(ep),    intent(in) :: d(*)
            complex(ep), intent(in) :: e(*)
            real(ep) :: r
        end function zlanht_quad

        function dlansy_quad(norm, uplo, n, A, lda, work) result(r)
            import :: ep
            character, intent(in) :: norm, uplo
            integer,   intent(in) :: n, lda
            real(ep),  intent(in) :: A(lda,*)
            real(ep),  intent(out) :: work(*)
            real(ep) :: r
        end function dlansy_quad

        function zlansy_quad(norm, uplo, n, A, lda, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo
            integer,     intent(in) :: n, lda
            complex(ep), intent(in) :: A(lda,*)
            real(ep),    intent(out) :: work(*)
            real(ep) :: r
        end function zlansy_quad

        function zlansb_quad(norm, uplo, n, k, AB, ldab, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo
            integer,     intent(in) :: n, k, ldab
            complex(ep), intent(in) :: AB(ldab,*)
            real(ep),    intent(out) :: work(*)
            real(ep) :: r
        end function zlansb_quad

        function zlansp_quad(norm, uplo, n, AP, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo
            integer,     intent(in) :: n
            complex(ep), intent(in) :: AP(*)
            real(ep),    intent(out) :: work(*)
            real(ep) :: r
        end function zlansp_quad

        function dlansf_quad(norm, transr, uplo, n, A, work) result(r)
            import :: ep
            character, intent(in) :: norm, transr, uplo
            integer,   intent(in) :: n
            real(ep),  intent(in) :: A(0:*)
            real(ep),  intent(out) :: work(0:*)
            real(ep) :: r
        end function dlansf_quad

        function dlantb_quad(norm, uplo, diag, n, k, AB, ldab, work) result(r)
            import :: ep
            character, intent(in) :: norm, uplo, diag
            integer,   intent(in) :: n, k, ldab
            real(ep),  intent(in) :: AB(ldab,*)
            real(ep),  intent(out) :: work(*)
            real(ep) :: r
        end function dlantb_quad

        function zlantb_quad(norm, uplo, diag, n, k, AB, ldab, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo, diag
            integer,     intent(in) :: n, k, ldab
            complex(ep), intent(in) :: AB(ldab,*)
            real(ep),    intent(out) :: work(*)
            real(ep) :: r
        end function zlantb_quad

        function dlantp_quad(norm, uplo, diag, n, AP, work) result(r)
            import :: ep
            character, intent(in) :: norm, uplo, diag
            integer,   intent(in) :: n
            real(ep),  intent(in) :: AP(*)
            real(ep),  intent(out) :: work(*)
            real(ep) :: r
        end function dlantp_quad

        function zlantp_quad(norm, uplo, diag, n, AP, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo, diag
            integer,     intent(in) :: n
            complex(ep), intent(in) :: AP(*)
            real(ep),    intent(out) :: work(*)
            real(ep) :: r
        end function zlantp_quad

        function dlantr_quad(norm, uplo, diag, m, n, A, lda, work) result(r)
            import :: ep
            character, intent(in) :: norm, uplo, diag
            integer,   intent(in) :: m, n, lda
            real(ep),  intent(in) :: A(lda,*)
            real(ep),  intent(out) :: work(*)
            real(ep) :: r
        end function dlantr_quad

        function zlantr_quad(norm, uplo, diag, m, n, A, lda, work) result(r)
            import :: ep
            character,   intent(in) :: norm, uplo, diag
            integer,     intent(in) :: m, n, lda
            complex(ep), intent(in) :: A(lda,*)
            real(ep),    intent(out) :: work(*)
            real(ep) :: r
        end function zlantr_quad

        ! Phase L19 — Permutation/norm helpers.
        subroutine dlapmr_quad(forwrd, m, n, X, ldx, K)
            import :: ep
            logical,  intent(in)    :: forwrd
            integer,  intent(in)    :: m, n, ldx
            real(ep), intent(inout) :: X(ldx,*)
            integer,  intent(inout) :: K(*)
        end subroutine dlapmr_quad

        subroutine zlapmr_quad(forwrd, m, n, X, ldx, K)
            import :: ep
            logical,     intent(in)    :: forwrd
            integer,     intent(in)    :: m, n, ldx
            complex(ep), intent(inout) :: X(ldx,*)
            integer,     intent(inout) :: K(*)
        end subroutine zlapmr_quad

        subroutine dlapmt_quad(forwrd, m, n, X, ldx, K)
            import :: ep
            logical,  intent(in)    :: forwrd
            integer,  intent(in)    :: m, n, ldx
            real(ep), intent(inout) :: X(ldx,*)
            integer,  intent(inout) :: K(*)
        end subroutine dlapmt_quad

        subroutine zlapmt_quad(forwrd, m, n, X, ldx, K)
            import :: ep
            logical,     intent(in)    :: forwrd
            integer,     intent(in)    :: m, n, ldx
            complex(ep), intent(inout) :: X(ldx,*)
            integer,     intent(inout) :: K(*)
        end subroutine zlapmt_quad

        subroutine dlapll_quad(n, X, incx, Y, incy, ssmin)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(inout) :: X(*), Y(*)
            real(ep), intent(out)   :: ssmin
        end subroutine dlapll_quad

        subroutine zlapll_quad(n, X, incx, Y, incy, ssmin)
            import :: ep
            integer,     intent(in)    :: n, incx, incy
            complex(ep), intent(inout) :: X(*), Y(*)
            real(ep),    intent(out)   :: ssmin
        end subroutine zlapll_quad

        subroutine dlacn2_quad(n, V, X, isgn, est, kase, isave)
            import :: ep
            integer,  intent(in)    :: n
            real(ep), intent(out)   :: V(*)
            real(ep), intent(inout) :: X(*), est
            integer,  intent(out)   :: isgn(*)
            integer,  intent(inout) :: kase, isave(3)
        end subroutine dlacn2_quad

        subroutine zlacn2_quad(n, V, X, est, kase, isave)
            import :: ep
            integer,     intent(in)    :: n
            complex(ep), intent(out)   :: V(*)
            complex(ep), intent(inout) :: X(*)
            real(ep),    intent(inout) :: est
            integer,     intent(inout) :: kase, isave(3)
        end subroutine zlacn2_quad

        subroutine dlacon_quad(n, V, X, isgn, est, kase)
            import :: ep
            integer,  intent(in)    :: n
            real(ep), intent(out)   :: V(*)
            real(ep), intent(inout) :: X(*), est
            integer,  intent(out)   :: isgn(*)
            integer,  intent(inout) :: kase
        end subroutine dlacon_quad

        subroutine zlacon_quad(n, V, X, est, kase)
            import :: ep
            integer,     intent(in)    :: n
            complex(ep), intent(out)   :: V(*)
            complex(ep), intent(inout) :: X(*)
            real(ep),    intent(inout) :: est
            integer,     intent(inout) :: kase
        end subroutine zlacon_quad

        subroutine dlartg_quad(f, g, c, s, r)
            import :: ep
            real(ep), intent(in)  :: f, g
            real(ep), intent(out) :: c, s, r
        end subroutine dlartg_quad

        subroutine zlartg_quad(f, g, c, s, r)
            import :: ep
            complex(ep), intent(in)  :: f, g
            real(ep),    intent(out) :: c
            complex(ep), intent(out) :: s, r
        end subroutine zlartg_quad

        subroutine dlartgp_quad(f, g, c, s, r)
            import :: ep
            real(ep), intent(in)  :: f, g
            real(ep), intent(out) :: c, s, r
        end subroutine dlartgp_quad

        subroutine dlartgs_quad(x, y, sigma, c, s)
            import :: ep
            real(ep), intent(in)  :: x, y, sigma
            real(ep), intent(out) :: c, s
        end subroutine dlartgs_quad

        ! Phase L20 — Small public scalar utilities.
        function dlapy2_quad(x, y) result(r)
            import :: ep
            real(ep), intent(in) :: x, y
            real(ep) :: r
        end function dlapy2_quad

        function dlapy3_quad(x, y, z) result(r)
            import :: ep
            real(ep), intent(in) :: x, y, z
            real(ep) :: r
        end function dlapy3_quad

        subroutine dladiv_quad(a, b, c, d, p, q)
            import :: ep
            real(ep), intent(in)  :: a, b, c, d
            real(ep), intent(out) :: p, q
        end subroutine dladiv_quad

        function zladiv_quad(x, y) result(r)
            import :: ep
            complex(ep), intent(in) :: x, y
            complex(ep) :: r
        end function zladiv_quad

        subroutine dlamrg_quad(n1, n2, A, dtrd1, dtrd2, index)
            import :: ep
            integer,  intent(in)  :: n1, n2, dtrd1, dtrd2
            real(ep), intent(in)  :: A(*)
            integer,  intent(out) :: index(*)
        end subroutine dlamrg_quad

        subroutine dlarnv_quad(idist, iseed, n, X)
            import :: ep
            integer,  intent(in)    :: idist, n
            integer,  intent(inout) :: iseed(4)
            real(ep), intent(out)   :: X(*)
        end subroutine dlarnv_quad

        subroutine zlarnv_quad(idist, iseed, n, X)
            import :: ep
            integer,     intent(in)    :: idist, n
            integer,     intent(inout) :: iseed(4)
            complex(ep), intent(out)   :: X(*)
        end subroutine zlarnv_quad

        ! Phase L21 — Generalized sym/Hermitian glue.
        subroutine dsbgst_quad(vect, uplo, n, ka, kb, AB, ldab, BB, ldbb, &
                          X, ldx, work, info)
            import :: ep
            character, intent(in)    :: vect, uplo
            integer,   intent(in)    :: n, ka, kb, ldab, ldbb, ldx
            real(ep),  intent(inout) :: AB(ldab,*)
            real(ep),  intent(in)    :: BB(ldbb,*)
            real(ep),  intent(out)   :: X(ldx,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsbgst_quad

        subroutine zhbgst_quad(vect, uplo, n, ka, kb, AB, ldab, BB, ldbb, &
                          X, ldx, work, rwork, info)
            import :: ep
            character,   intent(in)    :: vect, uplo
            integer,     intent(in)    :: n, ka, kb, ldab, ldbb, ldx
            complex(ep), intent(inout) :: AB(ldab,*)
            complex(ep), intent(in)    :: BB(ldbb,*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zhbgst_quad

        subroutine dspgst_quad(itype, uplo, n, AP, BP, info)
            import :: ep
            integer,   intent(in)    :: itype, n
            character, intent(in)    :: uplo
            real(ep),  intent(inout) :: AP(*)
            real(ep),  intent(in)    :: BP(*)
            integer,   intent(out)   :: info
        end subroutine dspgst_quad

        subroutine zhpgst_quad(itype, uplo, n, AP, BP, info)
            import :: ep
            integer,     intent(in)    :: itype, n
            character,   intent(in)    :: uplo
            complex(ep), intent(inout) :: AP(*)
            complex(ep), intent(in)    :: BP(*)
            integer,     intent(out)   :: info
        end subroutine zhpgst_quad

        subroutine dsygst_quad(itype, uplo, n, A, lda, B, ldb, info)
            import :: ep
            integer,   intent(in)    :: itype, n, lda, ldb
            character, intent(in)    :: uplo
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(in)    :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dsygst_quad

        subroutine zhegst_quad(itype, uplo, n, A, lda, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: itype, n, lda, ldb
            character,   intent(in)    :: uplo
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zhegst_quad

        subroutine dsbgvx_quad(jobz, range, uplo, n, ka, kb, AB, ldab, BB, ldbb, &
                          Q, ldq, vl, vu, il, iu, abstol, m, w, Z, ldz, &
                          work, iwork, ifail, info)
            import :: ep
            character, intent(in)    :: jobz, range, uplo
            integer,   intent(in)    :: n, ka, kb, ldab, ldbb, ldq, il, iu, ldz
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: AB(ldab,*), BB(ldbb,*)
            real(ep),  intent(out)   :: Q(ldq,*), w(*), Z(ldz,*), work(*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine dsbgvx_quad

        subroutine zhbgvx_quad(jobz, range, uplo, n, ka, kb, AB, ldab, BB, ldbb, &
                          Q, ldq, vl, vu, il, iu, abstol, m, w, Z, ldz, &
                          work, rwork, iwork, ifail, info)
            import :: ep
            character,   intent(in)    :: jobz, range, uplo
            integer,     intent(in)    :: n, ka, kb, ldab, ldbb, ldq, il, iu, ldz
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: AB(ldab,*), BB(ldbb,*)
            complex(ep), intent(out)   :: Q(ldq,*), Z(ldz,*), work(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            integer,     intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine zhbgvx_quad

        subroutine dspgvx_quad(itype, jobz, range, uplo, n, AP, BP, vl, vu, il, iu, &
                          abstol, m, w, Z, ldz, work, iwork, ifail, info)
            import :: ep
            integer,   intent(in)    :: itype, n, il, iu, ldz
            character, intent(in)    :: jobz, range, uplo
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: AP(*), BP(*)
            real(ep),  intent(out)   :: w(*), Z(ldz,*), work(*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine dspgvx_quad

        subroutine zhpgvx_quad(itype, jobz, range, uplo, n, AP, BP, vl, vu, il, iu, &
                          abstol, m, w, Z, ldz, work, rwork, iwork, ifail, info)
            import :: ep
            integer,     intent(in)    :: itype, n, il, iu, ldz
            character,   intent(in)    :: jobz, range, uplo
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: AP(*), BP(*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: Z(ldz,*), work(*)
            integer,     intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine zhpgvx_quad

        subroutine dsygvx_quad(itype, jobz, range, uplo, n, A, lda, B, ldb, &
                          vl, vu, il, iu, abstol, m, w, Z, ldz, &
                          work, lwork, iwork, ifail, info)
            import :: ep
            integer,   intent(in)    :: itype, n, lda, ldb, il, iu, ldz, lwork
            character, intent(in)    :: jobz, range, uplo
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: w(*), Z(ldz,*), work(*)
            integer,   intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine dsygvx_quad

        subroutine zhegvx_quad(itype, jobz, range, uplo, n, A, lda, B, ldb, &
                          vl, vu, il, iu, abstol, m, w, Z, ldz, &
                          work, lwork, rwork, iwork, ifail, info)
            import :: ep
            integer,     intent(in)    :: itype, n, lda, ldb, il, iu, ldz, lwork
            character,   intent(in)    :: jobz, range, uplo
            real(ep),    intent(in)    :: vl, vu, abstol
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),    intent(out)   :: w(*), rwork(*)
            complex(ep), intent(out)   :: Z(ldz,*), work(*)
            integer,     intent(out)   :: m, iwork(*), ifail(*), info
        end subroutine zhegvx_quad

        ! P12 — RFP / packed / triangular-full storage conversion
        subroutine dtfttp_quad(transr, uplo, n, ARF, AP, info)
            import :: ep
            character, intent(in)  :: transr, uplo
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: ARF(0:*)
            real(ep),  intent(out) :: AP(0:*)
            integer,   intent(out) :: info
        end subroutine dtfttp_quad

        subroutine ztfttp_quad(transr, uplo, n, ARF, AP, info)
            import :: ep
            character,   intent(in)  :: transr, uplo
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: ARF(0:*)
            complex(ep), intent(out) :: AP(0:*)
            integer,     intent(out) :: info
        end subroutine ztfttp_quad

        subroutine dtfttr_quad(transr, uplo, n, ARF, A, lda, info)
            import :: ep
            character, intent(in)    :: transr, uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(in)    :: ARF(0:*)
            real(ep),  intent(inout) :: A(0:lda-1, 0:*)
            integer,   intent(out)   :: info
        end subroutine dtfttr_quad

        subroutine ztfttr_quad(transr, uplo, n, ARF, A, lda, info)
            import :: ep
            character,   intent(in)    :: transr, uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(in)    :: ARF(0:*)
            complex(ep), intent(inout) :: A(0:lda-1, 0:*)
            integer,     intent(out)   :: info
        end subroutine ztfttr_quad

        subroutine dtpttf_quad(transr, uplo, n, AP, ARF, info)
            import :: ep
            character, intent(in)  :: transr, uplo
            integer,   intent(in)  :: n
            real(ep),  intent(in)  :: AP(0:*)
            real(ep),  intent(out) :: ARF(0:*)
            integer,   intent(out) :: info
        end subroutine dtpttf_quad

        subroutine ztpttf_quad(transr, uplo, n, AP, ARF, info)
            import :: ep
            character,   intent(in)  :: transr, uplo
            integer,     intent(in)  :: n
            complex(ep), intent(in)  :: AP(0:*)
            complex(ep), intent(out) :: ARF(0:*)
            integer,     intent(out) :: info
        end subroutine ztpttf_quad

        subroutine dtpttr_quad(uplo, n, AP, A, lda, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(in)    :: AP(*)
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: info
        end subroutine dtpttr_quad

        subroutine ztpttr_quad(uplo, n, AP, A, lda, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(in)    :: AP(*)
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: info
        end subroutine ztpttr_quad

        subroutine dtrttf_quad(transr, uplo, n, A, lda, ARF, info)
            import :: ep
            character, intent(in)  :: transr, uplo
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(0:lda-1, 0:*)
            real(ep),  intent(out) :: ARF(0:*)
            integer,   intent(out) :: info
        end subroutine dtrttf_quad

        subroutine ztrttf_quad(transr, uplo, n, A, lda, ARF, info)
            import :: ep
            character,   intent(in)  :: transr, uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(0:lda-1, 0:*)
            complex(ep), intent(out) :: ARF(0:*)
            integer,     intent(out) :: info
        end subroutine ztrttf_quad

        subroutine dtrttp_quad(uplo, n, A, lda, AP, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*)
            real(ep),  intent(out) :: AP(*)
            integer,   intent(out) :: info
        end subroutine dtrttp_quad

        subroutine ztrttp_quad(uplo, n, A, lda, AP, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            complex(ep), intent(out) :: AP(*)
            integer,     intent(out) :: info
        end subroutine ztrttp_quad

        ! P-misc — small utilities & unblocked apply/generate Householder

        logical function disnan_quad(din)
            import :: ep
            real(ep), intent(in) :: din
        end function disnan_quad

        subroutine drscl_quad(n, sa, sx, incx)
            import :: ep
            integer,  intent(in)    :: n, incx
            real(ep), intent(in)    :: sa
            real(ep), intent(inout) :: sx(*)
        end subroutine drscl_quad

        subroutine zdrscl_quad(n, sa, sx, incx)
            import :: ep
            integer,     intent(in)    :: n, incx
            real(ep),    intent(in)    :: sa
            complex(ep), intent(inout) :: sx(*)
        end subroutine zdrscl_quad

        subroutine zrscl_quad(n, a, x, incx)
            import :: ep
            integer,     intent(in)    :: n, incx
            complex(ep), intent(in)    :: a
            complex(ep), intent(inout) :: x(*)
        end subroutine zrscl_quad

        subroutine zrot_quad(n, cx, incx, cy, incy, c, s)
            import :: ep
            integer,     intent(in)    :: n, incx, incy
            real(ep),    intent(in)    :: c
            complex(ep), intent(in)    :: s
            complex(ep), intent(inout) :: cx(*), cy(*)
        end subroutine zrot_quad

        subroutine dorg2l_quad(m, n, k, A, lda, tau, work, info)
            import :: ep
            integer,  intent(in)    :: m, n, k, lda
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorg2l_quad

        subroutine zung2l_quad(m, n, k, A, lda, tau, work, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zung2l_quad

        subroutine dorg2r_quad(m, n, k, A, lda, tau, work, info)
            import :: ep
            integer,  intent(in)    :: m, n, k, lda
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: tau(*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorg2r_quad

        subroutine zung2r_quad(m, n, k, A, lda, tau, work, info)
            import :: ep
            integer,     intent(in)    :: m, n, k, lda
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: tau(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zung2r_quad

        subroutine dorm2l_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dorm2l_quad

        subroutine zunm2l_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunm2l_quad

        subroutine dorm2r_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dorm2r_quad

        subroutine zunm2r_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunm2r_quad

        subroutine dormlq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, lda, ldc, lwork
            real(ep),  intent(in)    :: A(lda,*), tau(*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dormlq_quad

        subroutine zunmlq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, lda, ldc, lwork
            complex(ep), intent(in)    :: A(lda,*), tau(*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunmlq_quad

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

        subroutine dtzrzf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dtzrzf_quad

        subroutine ztzrzf_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine ztzrzf_quad

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

        ! P11 — RFP factor / solve / triangular ops
        subroutine dpftrf_quad(transr, uplo, n, A, info)
            import :: ep
            character, intent(in)    :: transr, uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: A(0:*)
            integer,   intent(out)   :: info
        end subroutine dpftrf_quad

        subroutine zpftrf_quad(transr, uplo, n, A, info)
            import :: ep
            character,   intent(in)    :: transr, uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: A(0:*)
            integer,     intent(out)   :: info
        end subroutine zpftrf_quad

        subroutine dpftri_quad(transr, uplo, n, A, info)
            import :: ep
            character, intent(in)    :: transr, uplo
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: A(0:*)
            integer,   intent(out)   :: info
        end subroutine dpftri_quad

        subroutine zpftri_quad(transr, uplo, n, A, info)
            import :: ep
            character,   intent(in)    :: transr, uplo
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: A(0:*)
            integer,     intent(out)   :: info
        end subroutine zpftri_quad

        subroutine dpftrs_quad(transr, uplo, n, nrhs, A, B, ldb, info)
            import :: ep
            character, intent(in)    :: transr, uplo
            integer,   intent(in)    :: n, nrhs, ldb
            real(ep),  intent(in)    :: A(0:*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dpftrs_quad

        subroutine zpftrs_quad(transr, uplo, n, nrhs, A, B, ldb, info)
            import :: ep
            character,   intent(in)    :: transr, uplo
            integer,     intent(in)    :: n, nrhs, ldb
            complex(ep), intent(in)    :: A(0:*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zpftrs_quad

        subroutine dpstrf_quad(uplo, n, A, lda, piv, rank, tol, work, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: piv(*), rank, info
            real(ep),  intent(in)    :: tol
            real(ep),  intent(out)   :: work(*)
        end subroutine dpstrf_quad

        subroutine zpstrf_quad(uplo, n, A, lda, piv, rank, tol, work, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: piv(*), rank, info
            real(ep),    intent(in)    :: tol
            real(ep),    intent(out)   :: work(*)
        end subroutine zpstrf_quad

        subroutine dpteqr_quad(compz, n, D, E, Z, ldz, work, info)
            import :: ep
            character, intent(in)    :: compz
            integer,   intent(in)    :: n, ldz
            real(ep),  intent(inout) :: D(*), E(*), Z(ldz,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dpteqr_quad

        subroutine zpteqr_quad(compz, n, D, E, Z, ldz, work, info)
            import :: ep
            character,   intent(in)    :: compz
            integer,     intent(in)    :: n, ldz
            real(ep),    intent(inout) :: D(*), E(*)
            complex(ep), intent(inout) :: Z(ldz,*)
            real(ep),    intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zpteqr_quad

        subroutine dpbstf_quad(uplo, n, kd, AB, ldab, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, kd, ldab
            real(ep),  intent(inout) :: AB(ldab,*)
            integer,   intent(out)   :: info
        end subroutine dpbstf_quad

        subroutine zpbstf_quad(uplo, n, kd, AB, ldab, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, kd, ldab
            complex(ep), intent(inout) :: AB(ldab,*)
            integer,     intent(out)   :: info
        end subroutine zpbstf_quad

        subroutine dsfrk_quad(transr, uplo, trans, n, k, alpha, A, lda, beta, C)
            import :: ep
            character, intent(in)    :: transr, uplo, trans
            integer,   intent(in)    :: n, k, lda
            real(ep),  intent(in)    :: alpha, beta, A(lda,*)
            real(ep),  intent(inout) :: C(0:*)
        end subroutine dsfrk_quad

        subroutine zhfrk_quad(transr, uplo, trans, n, k, alpha, A, lda, beta, C)
            import :: ep
            character,   intent(in)    :: transr, uplo, trans
            integer,     intent(in)    :: n, k, lda
            real(ep),    intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: C(0:*)
        end subroutine zhfrk_quad

        subroutine dtfsm_quad(transr, side, uplo, trans, diag, m, n, alpha, A, B, ldb)
            import :: ep
            character, intent(in)    :: transr, side, uplo, trans, diag
            integer,   intent(in)    :: m, n, ldb
            real(ep),  intent(in)    :: alpha, A(0:*)
            real(ep),  intent(inout) :: B(0:ldb-1, 0:*)
        end subroutine dtfsm_quad

        subroutine ztfsm_quad(transr, side, uplo, trans, diag, m, n, alpha, A, B, ldb)
            import :: ep
            character,   intent(in)    :: transr, side, uplo, trans, diag
            integer,     intent(in)    :: m, n, ldb
            complex(ep), intent(in)    :: alpha, A(0:*)
            complex(ep), intent(inout) :: B(0:ldb-1, 0:*)
        end subroutine ztfsm_quad

        subroutine dtftri_quad(transr, uplo, diag, n, A, info)
            import :: ep
            character, intent(in)    :: transr, uplo, diag
            integer,   intent(in)    :: n
            real(ep),  intent(inout) :: A(0:*)
            integer,   intent(out)   :: info
        end subroutine dtftri_quad

        subroutine ztftri_quad(transr, uplo, diag, n, A, info)
            import :: ep
            character,   intent(in)    :: transr, uplo, diag
            integer,     intent(in)    :: n
            complex(ep), intent(inout) :: A(0:*)
            integer,     intent(out)   :: info
        end subroutine ztftri_quad

        ! P13 — pentagonal QR/LQ
        subroutine dtpqrt_quad(m, n, l, nb, A, lda, B, ldb, T, ldt, work, info)
            import :: ep
            integer,  intent(in)    :: m, n, l, nb, lda, ldb, ldt
            real(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dtpqrt_quad

        subroutine ztpqrt_quad(m, n, l, nb, A, lda, B, ldb, T, ldt, work, info)
            import :: ep
            integer,     intent(in)    :: m, n, l, nb, lda, ldb, ldt
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine ztpqrt_quad

        subroutine dtpqrt2_quad(m, n, l, A, lda, B, ldb, T, ldt, info)
            import :: ep
            integer,  intent(in)    :: m, n, l, lda, ldb, ldt
            real(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
            integer,  intent(out)   :: info
        end subroutine dtpqrt2_quad

        subroutine ztpqrt2_quad(m, n, l, A, lda, B, ldb, T, ldt, info)
            import :: ep
            integer,     intent(in)    :: m, n, l, lda, ldb, ldt
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
            integer,     intent(out)   :: info
        end subroutine ztpqrt2_quad

        subroutine dtpmqrt_quad(side, trans, m, n, k, l, nb, V, ldv, T, ldt, &
                           A, lda, B, ldb, work, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, l, nb, ldv, ldt, lda, ldb
            real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dtpmqrt_quad

        subroutine ztpmqrt_quad(side, trans, m, n, k, l, nb, V, ldv, T, ldt, &
                           A, lda, B, ldb, work, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, l, nb, ldv, ldt, lda, ldb
            complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine ztpmqrt_quad

        subroutine dtplqt_quad(m, n, l, mb, A, lda, B, ldb, T, ldt, work, info)
            import :: ep
            integer,  intent(in)    :: m, n, l, mb, lda, ldb, ldt
            real(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dtplqt_quad

        subroutine ztplqt_quad(m, n, l, mb, A, lda, B, ldb, T, ldt, work, info)
            import :: ep
            integer,     intent(in)    :: m, n, l, mb, lda, ldb, ldt
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine ztplqt_quad

        subroutine dtplqt2_quad(m, n, l, A, lda, B, ldb, T, ldt, info)
            import :: ep
            integer,  intent(in)    :: m, n, l, lda, ldb, ldt
            real(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
            integer,  intent(out)   :: info
        end subroutine dtplqt2_quad

        subroutine ztplqt2_quad(m, n, l, A, lda, B, ldb, T, ldt, info)
            import :: ep
            integer,     intent(in)    :: m, n, l, lda, ldb, ldt
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
            integer,     intent(out)   :: info
        end subroutine ztplqt2_quad

        subroutine dtpmlqt_quad(side, trans, m, n, k, l, mb, V, ldv, T, ldt, &
                           A, lda, B, ldb, work, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, k, l, mb, ldv, ldt, lda, ldb
            real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dtpmlqt_quad

        subroutine ztpmlqt_quad(side, trans, m, n, k, l, mb, V, ldv, T, ldt, &
                           A, lda, B, ldb, work, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, k, l, mb, ldv, ldt, lda, ldb
            complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine ztpmlqt_quad

        subroutine dtprfb_quad(side, trans, direct, storev, m, n, k, l, &
                          V, ldv, T, ldt, A, lda, B, ldb, work, ldwork)
            import :: ep
            character, intent(in)    :: side, trans, direct, storev
            integer,   intent(in)    :: m, n, k, l, ldv, ldt, lda, ldb, ldwork
            real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: work(ldwork,*)
        end subroutine dtprfb_quad

        subroutine ztprfb_quad(side, trans, direct, storev, m, n, k, l, &
                          V, ldv, T, ldt, A, lda, B, ldb, work, ldwork)
            import :: ep
            character,   intent(in)    :: side, trans, direct, storev
            integer,     intent(in)    :: m, n, k, l, ldv, ldt, lda, ldb, ldwork
            complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: work(ldwork,*)
        end subroutine ztprfb_quad

        ! P14 — TSQR + Householder reconstruction
        subroutine dgeqr2p_quad(m, n, A, lda, tau, work, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqr2p_quad

        subroutine zgeqr2p_quad(m, n, A, lda, tau, work, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqr2p_quad

        subroutine dgeqrfp_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, lda, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: tau(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgeqrfp_quad

        subroutine zgeqrfp_quad(m, n, A, lda, tau, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgeqrfp_quad

        subroutine dlatsqr_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, mb, nb, lda, ldt, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: T(ldt,*), work(*)
            integer,  intent(out)   :: info
        end subroutine dlatsqr_quad

        subroutine zlatsqr_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, mb, nb, lda, ldt, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: T(ldt,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zlatsqr_quad

        subroutine dorgtsqr_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, mb, nb, lda, ldt, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: T(ldt,*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorgtsqr_quad

        subroutine zungtsqr_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, mb, nb, lda, ldt, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: T(ldt,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungtsqr_quad

        subroutine dorgtsqr_row_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, mb, nb, lda, ldt, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(in)    :: T(ldt,*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorgtsqr_row_quad

        subroutine zungtsqr_row_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, mb, nb, lda, ldt, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: T(ldt,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zungtsqr_row_quad

        subroutine dorhr_col_quad(m, n, nb, A, lda, T, ldt, D, info)
            import :: ep
            integer,  intent(in)    :: m, n, nb, lda, ldt
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: T(ldt,*), D(*)
            integer,  intent(out)   :: info
        end subroutine dorhr_col_quad

        subroutine zunhr_col_quad(m, n, nb, A, lda, T, ldt, D, info)
            import :: ep
            integer,     intent(in)    :: m, n, nb, lda, ldt
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: T(ldt,*), D(*)
            integer,     intent(out)   :: info
        end subroutine zunhr_col_quad

        subroutine dgetsqrhrt_quad(m, n, mb1, nb1, nb2, A, lda, T, ldt, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, mb1, nb1, nb2, lda, ldt, lwork
            real(ep), intent(inout) :: A(lda,*)
            real(ep), intent(out)   :: T(ldt,*), work(*)
            integer,  intent(out)   :: info
        end subroutine dgetsqrhrt_quad

        subroutine zgetsqrhrt_quad(m, n, mb1, nb1, nb2, A, lda, T, ldt, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, mb1, nb1, nb2, lda, ldt, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: T(ldt,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgetsqrhrt_quad

        subroutine dgeqp3rk_quad(m, n, nrhs, kmax, abstol, reltol, A, lda, &
                            K, maxc2nrmk, relmaxc2nrmk, jpiv, tau, &
                            work, lwork, iwork, info)
            import :: ep
            integer,  intent(in)    :: m, n, nrhs, kmax, lda, lwork
            real(ep), intent(in)    :: abstol, reltol
            real(ep), intent(inout) :: A(lda,*)
            integer,  intent(out)   :: K, jpiv(*), iwork(*), info
            real(ep), intent(out)   :: maxc2nrmk, relmaxc2nrmk, tau(*), work(*)
        end subroutine dgeqp3rk_quad

        subroutine zgeqp3rk_quad(m, n, nrhs, kmax, abstol, reltol, A, lda, &
                            K, maxc2nrmk, relmaxc2nrmk, jpiv, tau, &
                            work, lwork, rwork, iwork, info)
            import :: ep
            integer,     intent(in)    :: m, n, nrhs, kmax, lda, lwork
            real(ep),    intent(in)    :: abstol, reltol
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: K, jpiv(*), iwork(*), info
            real(ep),    intent(out)   :: maxc2nrmk, relmaxc2nrmk, rwork(*)
            complex(ep), intent(out)   :: tau(*), work(*)
        end subroutine zgeqp3rk_quad

        ! P15 — CS decomposition + bdsvdx + bbcsd
        subroutine dbdsvdx_quad(uplo, jobz, range, n, D, E, vl, vu, il, iu, &
                           ns, S, Z, ldz, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo, jobz, range
            integer,   intent(in)  :: n, il, iu, ldz
            real(ep),  intent(in)  :: D(*), E(*), vl, vu
            integer,   intent(out) :: ns, iwork(*), info
            real(ep),  intent(out) :: S(*), Z(ldz,*), work(*)
        end subroutine dbdsvdx_quad

        subroutine dbbcsd_quad(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, &
                          theta, phi, U1, ldu1, U2, ldu2, V1t, ldv1t, V2t, ldv2t, &
                          B11d, B11e, B12d, B12e, B21d, B21e, B22d, B22e, &
                          work, lwork, info)
            import :: ep
            character, intent(in)    :: jobu1, jobu2, jobv1t, jobv2t, trans
            integer,   intent(in)    :: m, p, q, ldu1, ldu2, ldv1t, ldv2t, lwork
            real(ep),  intent(inout) :: theta(*), phi(*)
            real(ep),  intent(inout) :: U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*), V2t(ldv2t,*)
            real(ep),  intent(out)   :: B11d(*), B11e(*), B12d(*), B12e(*)
            real(ep),  intent(out)   :: B21d(*), B21e(*), B22d(*), B22e(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dbbcsd_quad

        subroutine zbbcsd_quad(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, &
                          theta, phi, U1, ldu1, U2, ldu2, V1t, ldv1t, V2t, ldv2t, &
                          B11d, B11e, B12d, B12e, B21d, B21e, B22d, B22e, &
                          rwork, lrwork, info)
            import :: ep
            character,   intent(in)    :: jobu1, jobu2, jobv1t, jobv2t, trans
            integer,     intent(in)    :: m, p, q, ldu1, ldu2, ldv1t, ldv2t, lrwork
            real(ep),    intent(inout) :: theta(*), phi(*)
            complex(ep), intent(inout) :: U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*), V2t(ldv2t,*)
            real(ep),    intent(out)   :: B11d(*), B11e(*), B12d(*), B12e(*)
            real(ep),    intent(out)   :: B21d(*), B21e(*), B22d(*), B22e(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zbbcsd_quad

        subroutine dorbdb_quad(trans, signs, m, p, q, X11, ldx11, X12, ldx12, &
                          X21, ldx21, X22, ldx22, theta, phi, &
                          taup1, taup2, tauq1, tauq2, work, lwork, info)
            import :: ep
            character, intent(in)    :: trans, signs
            integer,   intent(in)    :: m, p, q, ldx11, ldx12, ldx21, ldx22, lwork
            real(ep),  intent(inout) :: X11(ldx11,*), X12(ldx12,*), X21(ldx21,*), X22(ldx22,*)
            real(ep),  intent(out)   :: theta(*), phi(*), taup1(*), taup2(*), tauq1(*), tauq2(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dorbdb_quad

        subroutine zunbdb_quad(trans, signs, m, p, q, X11, ldx11, X12, ldx12, &
                          X21, ldx21, X22, ldx22, theta, phi, &
                          taup1, taup2, tauq1, tauq2, work, lwork, info)
            import :: ep
            character,   intent(in)    :: trans, signs
            integer,     intent(in)    :: m, p, q, ldx11, ldx12, ldx21, ldx22, lwork
            complex(ep), intent(inout) :: X11(ldx11,*), X12(ldx12,*), X21(ldx21,*), X22(ldx22,*)
            real(ep),    intent(out)   :: theta(*), phi(*)
            complex(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), tauq2(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zunbdb_quad

        recursive subroutine dorcsd_quad(jobu1, jobu2, jobv1t, jobv2t, trans, signs, &
                                    m, p, q, X11, ldx11, X12, ldx12, X21, ldx21, &
                                    X22, ldx22, theta, U1, ldu1, U2, ldu2, &
                                    V1t, ldv1t, V2t, ldv2t, work, lwork, iwork, info)
            import :: ep
            character, intent(in)    :: jobu1, jobu2, jobv1t, jobv2t, trans, signs
            integer,   intent(in)    :: m, p, q, ldx11, ldx12, ldx21, ldx22
            integer,   intent(in)    :: ldu1, ldu2, ldv1t, ldv2t, lwork
            real(ep),  intent(inout) :: X11(ldx11,*), X12(ldx12,*), X21(ldx21,*), X22(ldx22,*)
            real(ep),  intent(out)   :: theta(*)
            real(ep),  intent(out)   :: U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*), V2t(ldv2t,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dorcsd_quad

        recursive subroutine zuncsd_quad(jobu1, jobu2, jobv1t, jobv2t, trans, signs, &
                                    m, p, q, X11, ldx11, X12, ldx12, X21, ldx21, &
                                    X22, ldx22, theta, U1, ldu1, U2, ldu2, &
                                    V1t, ldv1t, V2t, ldv2t, &
                                    work, lwork, rwork, lrwork, iwork, info)
            import :: ep
            character,   intent(in)    :: jobu1, jobu2, jobv1t, jobv2t, trans, signs
            integer,     intent(in)    :: m, p, q, ldx11, ldx12, ldx21, ldx22
            integer,     intent(in)    :: ldu1, ldu2, ldv1t, ldv2t, lwork, lrwork
            complex(ep), intent(inout) :: X11(ldx11,*), X12(ldx12,*), X21(ldx21,*), X22(ldx22,*)
            real(ep),    intent(out)   :: theta(*)
            complex(ep), intent(out)   :: U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*), V2t(ldv2t,*)
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zuncsd_quad

        subroutine dorcsd2by1_quad(jobu1, jobu2, jobv1t, m, p, q, X11, ldx11, &
                              X21, ldx21, theta, U1, ldu1, U2, ldu2, V1t, ldv1t, &
                              work, lwork, iwork, info)
            import :: ep
            character, intent(in)    :: jobu1, jobu2, jobv1t
            integer,   intent(in)    :: m, p, q, ldx11, ldx21, ldu1, ldu2, ldv1t, lwork
            real(ep),  intent(inout) :: X11(ldx11,*), X21(ldx21,*)
            real(ep),  intent(out)   :: theta(*), U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dorcsd2by1_quad

        subroutine zuncsd2by1_quad(jobu1, jobu2, jobv1t, m, p, q, X11, ldx11, &
                              X21, ldx21, theta, U1, ldu1, U2, ldu2, V1t, ldv1t, &
                              work, lwork, rwork, lrwork, iwork, info)
            import :: ep
            character,   intent(in)    :: jobu1, jobu2, jobv1t
            integer,     intent(in)    :: m, p, q, ldx11, ldx21, ldu1, ldu2, ldv1t, lwork, lrwork
            complex(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
            real(ep),    intent(out)   :: theta(*)
            complex(ep), intent(out)   :: U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*)
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zuncsd2by1_quad

        subroutine dggsvd3_quad(jobu, jobv, jobq, m, n, p, k, l, A, lda, B, ldb, &
                           alpha, beta, U, ldu, V, ldv, Q, ldq, work, lwork, &
                           iwork, info)
            import :: ep
            character, intent(in)    :: jobu, jobv, jobq
            integer,   intent(in)    :: m, n, p, lda, ldb, ldu, ldv, ldq, lwork
            integer,   intent(out)   :: k, l
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: alpha(*), beta(*)
            real(ep),  intent(out)   :: U(ldu,*), V(ldv,*), Q(ldq,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dggsvd3_quad

        subroutine zggsvd3_quad(jobu, jobv, jobq, m, n, p, k, l, A, lda, B, ldb, &
                           alpha, beta, U, ldu, V, ldv, Q, ldq, work, lwork, &
                           rwork, iwork, info)
            import :: ep
            character,   intent(in)    :: jobu, jobv, jobq
            integer,     intent(in)    :: m, n, p, lda, ldb, ldu, ldv, ldq, lwork
            integer,     intent(out)   :: k, l
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),    intent(out)   :: alpha(*), beta(*)
            complex(ep), intent(out)   :: U(ldu,*), V(ldv,*), Q(ldq,*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zggsvd3_quad

        subroutine dggsvp3_quad(jobu, jobv, jobq, m, p, n, A, lda, B, ldb, tola, tolb, &
                           k, l, U, ldu, V, ldv, Q, ldq, iwork, tau, work, lwork, info)
            import :: ep
            character, intent(in)    :: jobu, jobv, jobq
            integer,   intent(in)    :: m, p, n, lda, ldb, ldu, ldv, ldq, lwork
            real(ep),  intent(in)    :: tola, tolb
            integer,   intent(out)   :: k, l
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: U(ldu,*), V(ldv,*), Q(ldq,*), tau(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dggsvp3_quad

        subroutine zggsvp3_quad(jobu, jobv, jobq, m, p, n, A, lda, B, ldb, tola, tolb, &
                           k, l, U, ldu, V, ldv, Q, ldq, iwork, rwork, tau, work, lwork, info)
            import :: ep
            character,   intent(in)    :: jobu, jobv, jobq
            integer,     intent(in)    :: m, p, n, lda, ldb, ldu, ldv, ldq, lwork
            real(ep),    intent(in)    :: tola, tolb
            integer,     intent(out)   :: k, l
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: U(ldu,*), V(ldv,*), Q(ldq,*), tau(*), work(*)
            real(ep),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: iwork(*), info
        end subroutine zggsvp3_quad

        subroutine dorbdb1_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, &
                           taup1, taup2, tauq1, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, p, q, ldx11, ldx21, lwork
            real(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
            real(ep), intent(out)   :: theta(*), phi(*)
            real(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dorbdb1_quad

        subroutine dorbdb2_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, &
                           taup1, taup2, tauq1, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, p, q, ldx11, ldx21, lwork
            real(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
            real(ep), intent(out)   :: theta(*), phi(*)
            real(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dorbdb2_quad

        subroutine dorbdb3_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, &
                           taup1, taup2, tauq1, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, p, q, ldx11, ldx21, lwork
            real(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
            real(ep), intent(out)   :: theta(*), phi(*)
            real(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dorbdb3_quad

        subroutine dorbdb4_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, &
                           taup1, taup2, tauq1, phantom, work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m, p, q, ldx11, ldx21, lwork
            real(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
            real(ep), intent(out)   :: theta(*), phi(*)
            real(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), phantom(*), work(*)
            integer,  intent(out)   :: info
        end subroutine dorbdb4_quad

        subroutine dorbdb5_quad(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, &
                           work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m1, m2, n, incx1, incx2, ldq1, ldq2, lwork
            real(ep), intent(inout) :: X1(*), X2(*)
            real(ep), intent(in)    :: Q1(ldq1,*), Q2(ldq2,*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorbdb5_quad

        subroutine dorbdb6_quad(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, &
                           work, lwork, info)
            import :: ep
            integer,  intent(in)    :: m1, m2, n, incx1, incx2, ldq1, ldq2, lwork
            real(ep), intent(inout) :: X1(*), X2(*)
            real(ep), intent(in)    :: Q1(ldq1,*), Q2(ldq2,*)
            real(ep), intent(out)   :: work(*)
            integer,  intent(out)   :: info
        end subroutine dorbdb6_quad

        subroutine zunbdb1_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, &
                           taup1, taup2, tauq1, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, p, q, ldx11, ldx21, lwork
            complex(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
            real(ep),    intent(out)   :: theta(*), phi(*)
            complex(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zunbdb1_quad

        subroutine zunbdb2_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, &
                           taup1, taup2, tauq1, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, p, q, ldx11, ldx21, lwork
            complex(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
            real(ep),    intent(out)   :: theta(*), phi(*)
            complex(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zunbdb2_quad

        subroutine zunbdb3_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, &
                           taup1, taup2, tauq1, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, p, q, ldx11, ldx21, lwork
            complex(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
            real(ep),    intent(out)   :: theta(*), phi(*)
            complex(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zunbdb3_quad

        subroutine zunbdb4_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, &
                           taup1, taup2, tauq1, phantom, work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m, p, q, ldx11, ldx21, lwork
            complex(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
            real(ep),    intent(out)   :: theta(*), phi(*)
            complex(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), phantom(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zunbdb4_quad

        subroutine zunbdb5_quad(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, &
                           work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m1, m2, n, incx1, incx2, ldq1, ldq2, lwork
            complex(ep), intent(inout) :: X1(*), X2(*)
            complex(ep), intent(in)    :: Q1(ldq1,*), Q2(ldq2,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunbdb5_quad

        subroutine zunbdb6_quad(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, &
                           work, lwork, info)
            import :: ep
            integer,     intent(in)    :: m1, m2, n, incx1, incx2, ldq1, ldq2, lwork
            complex(ep), intent(inout) :: X1(*), X2(*)
            complex(ep), intent(in)    :: Q1(ldq1,*), Q2(ldq2,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunbdb6_quad

        subroutine dorm22_quad(side, trans, m, n, n1, n2, Q, ldq, C, ldc, work, lwork, info)
            import :: ep
            character, intent(in)    :: side, trans
            integer,   intent(in)    :: m, n, n1, n2, ldq, ldc, lwork
            real(ep),  intent(in)    :: Q(ldq,*)
            real(ep),  intent(inout) :: C(ldc,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dorm22_quad

        subroutine zunm22_quad(side, trans, m, n, n1, n2, Q, ldq, C, ldc, work, lwork, info)
            import :: ep
            character,   intent(in)    :: side, trans
            integer,     intent(in)    :: m, n, n1, n2, ldq, ldc, lwork
            complex(ep), intent(in)    :: Q(ldq,*)
            complex(ep), intent(inout) :: C(ldc,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zunm22_quad

        ! P6/P7 — modern SVD + banded bidiag reduction
        subroutine dgesvdq_quad(joba, jobp, jobr, jobu, jobv, m, n, A, lda, &
                           S, U, ldu, V, ldv, numrank, iwork, liwork, &
                           work, lwork, rwork, lrwork, info)
            import :: ep
            character, intent(in)    :: joba, jobp, jobr, jobu, jobv
            integer,   intent(in)    :: m, n, lda, ldu, ldv, liwork, lwork, lrwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: S(*), U(ldu,*), V(ldv,*), work(*), rwork(*)
            integer,   intent(out)   :: numrank, iwork(*), info
        end subroutine dgesvdq_quad

        subroutine zgesvdq_quad(joba, jobp, jobr, jobu, jobv, m, n, A, lda, &
                           S, U, ldu, V, ldv, numrank, iwork, liwork, &
                           cwork, lcwork, rwork, lrwork, info)
            import :: ep
            character,   intent(in)    :: joba, jobp, jobr, jobu, jobv
            integer,     intent(in)    :: m, n, lda, ldu, ldv, liwork, lcwork, lrwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: S(*), rwork(*)
            complex(ep), intent(out)   :: U(ldu,*), V(ldv,*), cwork(*)
            integer,     intent(out)   :: numrank, iwork(*), info
        end subroutine zgesvdq_quad

        subroutine dgesvdx_quad(jobu, jobvt, range, m, n, A, lda, vl, vu, il, iu, &
                           ns, S, U, ldu, Vt, ldvt, work, lwork, iwork, info)
            import :: ep
            character, intent(in)    :: jobu, jobvt, range
            integer,   intent(in)    :: m, n, il, iu, lda, ldu, ldvt, lwork
            real(ep),  intent(in)    :: vl, vu
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: ns, iwork(*), info
            real(ep),  intent(out)   :: S(*), U(ldu,*), Vt(ldvt,*), work(*)
        end subroutine dgesvdx_quad

        subroutine zgesvdx_quad(jobu, jobvt, range, m, n, A, lda, vl, vu, il, iu, &
                           ns, S, U, ldu, Vt, ldvt, work, lwork, rwork, iwork, info)
            import :: ep
            character,   intent(in)    :: jobu, jobvt, range
            integer,     intent(in)    :: m, n, il, iu, lda, ldu, ldvt, lwork
            real(ep),    intent(in)    :: vl, vu
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ns, iwork(*), info
            real(ep),    intent(out)   :: S(*), rwork(*)
            complex(ep), intent(out)   :: U(ldu,*), Vt(ldvt,*), work(*)
        end subroutine zgesvdx_quad

        subroutine dgbbrd_quad(vect, m, n, ncc, kl, ku, AB, ldab, D, E, &
                          Q, ldq, Pt, ldpt, C, ldc, work, info)
            import :: ep
            character, intent(in)    :: vect
            integer,   intent(in)    :: m, n, ncc, kl, ku, ldab, ldq, ldpt, ldc
            real(ep),  intent(inout) :: AB(ldab,*), C(ldc,*)
            real(ep),  intent(out)   :: D(*), E(*), Q(ldq,*), Pt(ldpt,*), work(*)
            integer,   intent(out)   :: info
        end subroutine dgbbrd_quad

        subroutine zgbbrd_quad(vect, m, n, ncc, kl, ku, AB, ldab, D, E, &
                          Q, ldq, Pt, ldpt, C, ldc, work, rwork, info)
            import :: ep
            character,   intent(in)    :: vect
            integer,     intent(in)    :: m, n, ncc, kl, ku, ldab, ldq, ldpt, ldc
            complex(ep), intent(inout) :: AB(ldab,*), C(ldc,*)
            real(ep),    intent(out)   :: D(*), E(*), rwork(*)
            complex(ep), intent(out)   :: Q(ldq,*), Pt(ldpt,*), work(*)
            integer,     intent(out)   :: info
        end subroutine zgbbrd_quad

        ! P16 — MRRR / inverse-iteration tridiagonal eigensolvers
        subroutine dstegr_quad(jobz, range, n, D, E, vl, vu, il, iu, abstol, &
                          m, W, Z, ldz, isuppz, work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, range
            integer,   intent(in)    :: n, il, iu, ldz, lwork, liwork
            real(ep),  intent(in)    :: vl, vu, abstol
            real(ep),  intent(inout) :: D(*), E(*)
            integer,   intent(out)   :: m, isuppz(*), iwork(*), info
            real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
        end subroutine dstegr_quad

        subroutine zstegr_quad(jobz, range, n, D, E, vl, vu, il, iu, abstol, &
                          m, W, Z, ldz, isuppz, work, lwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, range
            integer,     intent(in)    :: n, il, iu, ldz, lwork, liwork
            real(ep),    intent(in)    :: vl, vu, abstol
            real(ep),    intent(inout) :: D(*), E(*)
            integer,     intent(out)   :: m, isuppz(*), iwork(*), info
            real(ep),    intent(out)   :: W(*), work(*)
            complex(ep), intent(out)   :: Z(ldz,*)
        end subroutine zstegr_quad

        subroutine dstemr_quad(jobz, range, n, D, E, vl, vu, il, iu, &
                          m, W, Z, ldz, nzc, isuppz, tryrac, &
                          work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobz, range
            integer,   intent(in)    :: n, il, iu, ldz, nzc, lwork, liwork
            logical,   intent(inout) :: tryrac
            real(ep),  intent(in)    :: vl, vu
            real(ep),  intent(inout) :: D(*), E(*)
            integer,   intent(out)   :: m, isuppz(*), iwork(*), info
            real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
        end subroutine dstemr_quad

        subroutine zstemr_quad(jobz, range, n, D, E, vl, vu, il, iu, &
                          m, W, Z, ldz, nzc, isuppz, tryrac, &
                          work, lwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobz, range
            integer,     intent(in)    :: n, il, iu, ldz, nzc, lwork, liwork
            logical,     intent(inout) :: tryrac
            real(ep),    intent(in)    :: vl, vu
            real(ep),    intent(inout) :: D(*), E(*)
            integer,     intent(out)   :: m, isuppz(*), iwork(*), info
            real(ep),    intent(out)   :: W(*), work(*)
            complex(ep), intent(out)   :: Z(ldz,*)
        end subroutine zstemr_quad

        subroutine dstein_quad(n, D, E, m, W, iblock, isplit, Z, ldz, work, iwork, ifail, info)
            import :: ep
            integer,  intent(in)    :: n, m, ldz
            real(ep), intent(in)    :: D(*), E(*), W(*)
            integer,  intent(in)    :: iblock(*), isplit(*)
            real(ep), intent(out)   :: Z(ldz,*), work(*)
            integer,  intent(out)   :: iwork(*), ifail(*), info
        end subroutine dstein_quad

        subroutine zstein_quad(n, D, E, m, W, iblock, isplit, Z, ldz, work, iwork, ifail, info)
            import :: ep
            integer,     intent(in)    :: n, m, ldz
            real(ep),    intent(in)    :: D(*), E(*), W(*)
            integer,     intent(in)    :: iblock(*), isplit(*)
            complex(ep), intent(out)   :: Z(ldz,*)
            real(ep),    intent(out)   :: work(*)
            integer,     intent(out)   :: iwork(*), ifail(*), info
        end subroutine zstein_quad

        ! ── Bunch–Kaufman _rook variants ───────────────────────────────
        subroutine dsytrf_rook_quad(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: ipiv(*), info
            real(ep),  intent(out)   :: work(*)
        end subroutine dsytrf_rook_quad

        subroutine dsytf2_rook_quad(uplo, n, A, lda, ipiv, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: ipiv(*), info
        end subroutine dsytf2_rook_quad

        subroutine dsytrs_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dsytrs_rook_quad

        subroutine dsytri_rook_quad(uplo, n, A, lda, ipiv, work, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dsytri_rook_quad

        subroutine dsycon_rook_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*), anorm
            integer,   intent(in)  :: ipiv(*)
            real(ep),  intent(out) :: rcond
            real(ep),  intent(out) :: work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dsycon_rook_quad

        subroutine dsysv_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: ipiv(*), info
            real(ep),  intent(out)   :: work(*)
        end subroutine dsysv_rook_quad

        subroutine zhetrf_rook_quad(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zhetrf_rook_quad

        subroutine zhetf2_rook_quad(uplo, n, A, lda, ipiv, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zhetf2_rook_quad

        subroutine zhetrs_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zhetrs_rook_quad

        subroutine zhetri_rook_quad(uplo, n, A, lda, ipiv, work, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhetri_rook_quad

        subroutine zhecon_rook_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zhecon_rook_quad

        subroutine zhesv_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zhesv_rook_quad

        ! ── Bunch–Kaufman _rk / _3x family ─────────────────────────────
        subroutine dsytrf_rk_quad(uplo, n, A, lda, e, ipiv, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: e(*), work(*)
            integer,   intent(out)   :: ipiv(*), info
        end subroutine dsytrf_rk_quad

        subroutine dsytf2_rk_quad(uplo, n, A, lda, e, ipiv, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: e(*)
            integer,   intent(out)   :: ipiv(*), info
        end subroutine dsytf2_rk_quad

        subroutine dsysv_rk_quad(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            real(ep),  intent(out)   :: e(*), work(*)
            integer,   intent(out)   :: ipiv(*), info
        end subroutine dsysv_rk_quad

        subroutine dsytri_3x_quad(uplo, n, A, lda, e, ipiv, work, nb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, nb
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(in)    :: e(*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(out)   :: work(n+nb+1, nb+3)
            integer,   intent(out)   :: info
        end subroutine dsytri_3x_quad

        subroutine dsytri_3_quad(uplo, n, A, lda, e, ipiv, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(in)    :: e(*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dsytri_3_quad

        subroutine zhetrf_rk_quad(uplo, n, A, lda, e, ipiv, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: e(*), work(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zhetrf_rk_quad

        subroutine zhetf2_rk_quad(uplo, n, A, lda, e, ipiv, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: e(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zhetf2_rk_quad

        subroutine zhesv_rk_quad(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: e(*), work(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zhesv_rk_quad

        subroutine zhetri_3x_quad(uplo, n, A, lda, e, ipiv, work, nb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, nb
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: e(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(n+nb+1, nb+3)
            integer,     intent(out)   :: info
        end subroutine zhetri_3x_quad

        subroutine zsytrf_rk_quad(uplo, n, A, lda, e, ipiv, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: e(*), work(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zsytrf_rk_quad

        subroutine zsytf2_rk_quad(uplo, n, A, lda, e, ipiv, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: e(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zsytf2_rk_quad

        subroutine zsysv_rk_quad(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            complex(ep), intent(out)   :: e(*), work(*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zsysv_rk_quad

        subroutine zsytri_3x_quad(uplo, n, A, lda, e, ipiv, work, nb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, nb
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: e(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(n+nb+1, nb+3)
            integer,     intent(out)   :: info
        end subroutine zsytri_3x_quad

        subroutine zsytri_3_quad(uplo, n, A, lda, e, ipiv, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(in)    :: e(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsytri_3_quad

        subroutine dsytri2_quad(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dsytri2_quad

        subroutine zsytri2_quad(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsytri2_quad

        subroutine dsytrs2_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dsytrs2_quad

        subroutine zsytrs2_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsytrs2_quad

        subroutine dsytrs_3_quad(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb
            real(ep),  intent(in)    :: A(lda,*), e(*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dsytrs_3_quad

        subroutine zsytrs_3_quad(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*), e(*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zsytrs_3_quad

        subroutine dsycon_3_quad(uplo, n, A, lda, e, ipiv, anorm, rcond, work, iwork, info)
            import :: ep
            character, intent(in)  :: uplo
            integer,   intent(in)  :: n, lda
            real(ep),  intent(in)  :: A(lda,*), e(*), anorm
            integer,   intent(in)  :: ipiv(*)
            real(ep),  intent(out) :: rcond, work(*)
            integer,   intent(out) :: iwork(*), info
        end subroutine dsycon_3_quad

        subroutine zsycon_3_quad(uplo, n, A, lda, e, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*), e(*)
            real(ep),    intent(in)  :: anorm
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zsycon_3_quad

        ! ── Aasen _aa / _aa_2stage family ──────────────────────────────
        subroutine dsytrf_aa_quad(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, lwork
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(out)   :: ipiv(*), info
            real(ep),  intent(out)   :: work(*)
        end subroutine dsytrf_aa_quad

        subroutine dsytrs_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb, lwork
            real(ep),  intent(in)    :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(inout) :: B(ldb,*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dsytrs_aa_quad

        subroutine dsysv_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
            integer,   intent(out)   :: ipiv(*), info
            real(ep),  intent(out)   :: work(*)
        end subroutine dsysv_aa_quad

        subroutine dsytrf_aa_2stage_quad(uplo, n, A, lda, TB, ltb, ipiv, ipiv2, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, ltb, lwork
            real(ep),  intent(inout) :: A(lda,*), TB(*)
            integer,   intent(out)   :: ipiv(*), ipiv2(*), info
            real(ep),  intent(out)   :: work(*)
        end subroutine dsytrf_aa_2stage_quad

        subroutine dsytrs_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ltb, ldb
            real(ep),  intent(in)    :: A(lda,*), TB(*)
            integer,   intent(in)    :: ipiv(*), ipiv2(*)
            real(ep),  intent(inout) :: B(ldb,*)
            integer,   intent(out)   :: info
        end subroutine dsytrs_aa_2stage_quad

        subroutine dsysv_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, nrhs, lda, ltb, ldb, lwork
            real(ep),  intent(inout) :: A(lda,*), B(ldb,*), TB(*)
            integer,   intent(out)   :: ipiv(*), ipiv2(*), info
            real(ep),  intent(out)   :: work(*)
        end subroutine dsysv_aa_2stage_quad

        subroutine zhetrf_aa_quad(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zhetrf_aa_quad

        subroutine zhetrs_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
            complex(ep), intent(in)    :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhetrs_aa_quad

        subroutine zhesv_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zhesv_aa_quad

        subroutine zhetrf_aa_2stage_quad(uplo, n, A, lda, TB, ltb, ipiv, ipiv2, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, ltb, lwork
            complex(ep), intent(inout) :: A(lda,*), TB(*)
            integer,     intent(out)   :: ipiv(*), ipiv2(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zhetrf_aa_2stage_quad

        subroutine zhetrs_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ltb, ldb
            complex(ep), intent(in)    :: A(lda,*), TB(*)
            integer,     intent(in)    :: ipiv(*), ipiv2(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zhetrs_aa_2stage_quad

        subroutine zhesv_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ltb, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), TB(*)
            integer,     intent(out)   :: ipiv(*), ipiv2(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zhesv_aa_2stage_quad

        subroutine zsytrf_aa_quad(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zsytrf_aa_quad

        subroutine zsytrs_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
            complex(ep), intent(in)    :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsytrs_aa_quad

        subroutine zsysv_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zsysv_aa_quad

        subroutine zsytrf_aa_2stage_quad(uplo, n, A, lda, TB, ltb, ipiv, ipiv2, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, ltb, lwork
            complex(ep), intent(inout) :: A(lda,*), TB(*)
            integer,     intent(out)   :: ipiv(*), ipiv2(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zsytrf_aa_2stage_quad

        subroutine zsytrs_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ltb, ldb
            complex(ep), intent(in)    :: A(lda,*), TB(*)
            integer,     intent(in)    :: ipiv(*), ipiv2(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zsytrs_aa_2stage_quad

        subroutine zsysv_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ltb, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*), TB(*)
            integer,     intent(out)   :: ipiv(*), ipiv2(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zsysv_aa_2stage_quad

        ! ── Z-symmetric base + Z-sym _rook + inverse helpers ────────────
        subroutine zsytrf_quad(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zsytrf_quad

        subroutine zsytrs_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zsytrs_quad

        subroutine zsysv_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zsysv_quad

        subroutine zsycon_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zsycon_quad

        subroutine zsytri_quad(uplo, n, A, lda, ipiv, work, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsytri_quad

        subroutine zsytrf_rook_quad(uplo, n, A, lda, ipiv, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, lwork
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zsytrf_rook_quad

        subroutine zsytf2_rook_quad(uplo, n, A, lda, ipiv, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zsytf2_rook_quad

        subroutine zsytrs_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(in)    :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: B(ldb,*)
            integer,     intent(out)   :: info
        end subroutine zsytrs_rook_quad

        subroutine zsytri_rook_quad(uplo, n, A, lda, ipiv, work, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zsytri_rook_quad

        subroutine zsycon_rook_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
            import :: ep
            character,   intent(in)  :: uplo
            integer,     intent(in)  :: n, lda
            complex(ep), intent(in)  :: A(lda,*)
            integer,     intent(in)  :: ipiv(*)
            real(ep),    intent(in)  :: anorm
            real(ep),    intent(out) :: rcond
            complex(ep), intent(out) :: work(*)
            integer,     intent(out) :: info
        end subroutine zsycon_rook_quad

        subroutine zsysv_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
            complex(ep), intent(out)   :: work(*)
        end subroutine zsysv_rook_quad

        subroutine dsytri_quad(uplo, n, A, lda, ipiv, work, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(out)   :: work(*)
            integer,   intent(out)   :: info
        end subroutine dsytri_quad

        subroutine zhetri_quad(uplo, n, A, lda, ipiv, work, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhetri_quad

        subroutine dsytri2x_quad(uplo, n, A, lda, ipiv, work, nb, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, nb
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(out)   :: work(n+nb+1, nb+3)
            integer,   intent(out)   :: info
        end subroutine dsytri2x_quad

        subroutine zhetri2x_quad(uplo, n, A, lda, ipiv, work, nb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, nb
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(n+nb+1, nb+3)
            integer,     intent(out)   :: info
        end subroutine zhetri2x_quad

        subroutine zsytri2x_quad(uplo, n, A, lda, ipiv, work, nb, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, nb
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: work(n+nb+1, nb+3)
            integer,     intent(out)   :: info
        end subroutine zsytri2x_quad

        ! ── BLAS-style Z-symmetric / Z-packed fillers ───────────────────
        subroutine zsymv_quad(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, incx, incy
            complex(ep), intent(in)    :: alpha, beta, A(lda,*), x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine zsymv_quad

        subroutine zsyr_quad(uplo, n, alpha, x, incx, A, lda)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, incx
            complex(ep), intent(in)    :: alpha, x(*)
            complex(ep), intent(inout) :: A(lda,*)
        end subroutine zsyr_quad

        subroutine zspmv_quad(uplo, n, alpha, AP, x, incx, beta, y, incy)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, incx, incy
            complex(ep), intent(in)    :: alpha, beta, AP(*), x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine zspmv_quad

        subroutine zspr_quad(uplo, n, alpha, x, incx, AP)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, incx
            complex(ep), intent(in)    :: alpha, x(*)
            complex(ep), intent(inout) :: AP(*)
        end subroutine zspr_quad

        subroutine dsyswapr_quad(uplo, n, A, lda, i1, i2)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, i1, i2
            real(ep),  intent(inout) :: A(lda,*)
        end subroutine dsyswapr_quad

        subroutine zheswapr_quad(uplo, n, A, lda, i1, i2)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, i1, i2
            complex(ep), intent(inout) :: A(lda,*)
        end subroutine zheswapr_quad

        subroutine zsyswapr_quad(uplo, n, A, lda, i1, i2)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, i1, i2
            complex(ep), intent(inout) :: A(lda,*)
        end subroutine zsyswapr_quad

        subroutine dsyconv_quad(uplo, way, n, A, lda, ipiv, e, info)
            import :: ep
            character, intent(in)    :: uplo, way
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*)
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(out)   :: e(*)
            integer,   intent(out)   :: info
        end subroutine dsyconv_quad

        subroutine dsyconvf_quad(uplo, way, n, A, lda, e, ipiv, info)
            import :: ep
            character, intent(in)    :: uplo, way
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*), e(*)
            integer,   intent(inout) :: ipiv(*)
            integer,   intent(out)   :: info
        end subroutine dsyconvf_quad

        subroutine dsyconvf_rook_quad(uplo, way, n, A, lda, e, ipiv, info)
            import :: ep
            character, intent(in)    :: uplo, way
            integer,   intent(in)    :: n, lda
            real(ep),  intent(inout) :: A(lda,*), e(*)
            integer,   intent(in)    :: ipiv(*)
            integer,   intent(out)   :: info
        end subroutine dsyconvf_rook_quad

        subroutine zsyconv_quad(uplo, way, n, A, lda, ipiv, e, info)
            import :: ep
            character,   intent(in)    :: uplo, way
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(out)   :: e(*)
            integer,     intent(out)   :: info
        end subroutine zsyconv_quad

        subroutine zsyconvf_quad(uplo, way, n, A, lda, e, ipiv, info)
            import :: ep
            character,   intent(in)    :: uplo, way
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*), e(*)
            integer,     intent(inout) :: ipiv(*)
            integer,     intent(out)   :: info
        end subroutine zsyconvf_quad

        subroutine zsyconvf_rook_quad(uplo, way, n, A, lda, e, ipiv, info)
            import :: ep
            character,   intent(in)    :: uplo, way
            integer,     intent(in)    :: n, lda
            complex(ep), intent(inout) :: A(lda,*), e(*)
            integer,     intent(in)    :: ipiv(*)
            integer,     intent(out)   :: info
        end subroutine zsyconvf_rook_quad

        subroutine zhprfs_quad(uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            complex(ep), intent(in)    :: AP(*), AFP(*), B(ldb,*)
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhprfs_quad

        subroutine zhpsvx_quad(fact, uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, &
                          rcond, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, uplo
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            complex(ep), intent(in)    :: AP(*), B(ldb,*)
            complex(ep), intent(inout) :: AFP(*)
            integer,     intent(inout) :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zhpsvx_quad

        subroutine zspsvx_quad(fact, uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, &
                          rcond, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, uplo
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            complex(ep), intent(in)    :: AP(*), B(ldb,*)
            complex(ep), intent(inout) :: AFP(*)
            integer,     intent(inout) :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zspsvx_quad

        subroutine dspsvx_quad(fact, uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, &
                          rcond, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: fact, uplo
            integer,   intent(in)    :: n, nrhs, ldb, ldx
            real(ep),  intent(in)    :: AP(*), B(ldb,*)
            real(ep),  intent(inout) :: AFP(*)
            integer,   intent(inout) :: ipiv(*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: rcond, ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dspsvx_quad

        subroutine dgtsv_quad(n, nrhs, dl, d, du, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, ldb
            real(ep), intent(inout) :: dl(*), d(*), du(*), B(ldb,*)
            integer,  intent(out)   :: info
        end subroutine dgtsv_quad

        subroutine dptsv_quad(n, nrhs, d, e, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, ldb
            real(ep), intent(inout) :: d(*), e(*), B(ldb,*)
            integer,  intent(out)   :: info
        end subroutine dptsv_quad

        ! ── gtsvx / ptsvx / trsna ────────────────────────────────────
        subroutine dgtsvx_quad(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, &
                          B, ldb, X, ldx, rcond, ferr, berr, work, iwork, info)
            import :: ep
            character, intent(in)    :: fact, trans
            integer,   intent(in)    :: n, nrhs, ldb, ldx
            real(ep),  intent(in)    :: dl(*), d(*), du(*), B(ldb,*)
            real(ep),  intent(inout) :: dlf(*), df(*), duf(*), du2(*)
            integer,   intent(inout) :: ipiv(*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: rcond, ferr(*), berr(*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgtsvx_quad

        subroutine zgtsvx_quad(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, &
                          B, ldb, X, ldx, rcond, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, trans
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            complex(ep), intent(in)    :: dl(*), d(*), du(*), B(ldb,*)
            complex(ep), intent(inout) :: dlf(*), df(*), duf(*), du2(*)
            integer,     intent(inout) :: ipiv(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zgtsvx_quad

        subroutine dptsvx_quad(fact, n, nrhs, d, e, df, ef, B, ldb, X, ldx, &
                          rcond, ferr, berr, work, info)
            import :: ep
            character, intent(in)    :: fact
            integer,   intent(in)    :: n, nrhs, ldb, ldx
            real(ep),  intent(in)    :: d(*), e(*), B(ldb,*)
            real(ep),  intent(inout) :: df(*), ef(*)
            real(ep),  intent(inout) :: X(ldx,*)
            real(ep),  intent(out)   :: rcond, ferr(*), berr(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dptsvx_quad

        subroutine zptsvx_quad(fact, n, nrhs, d, e, df, ef, B, ldb, X, ldx, &
                          rcond, ferr, berr, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact
            integer,     intent(in)    :: n, nrhs, ldb, ldx
            real(ep),    intent(in)    :: d(*)
            complex(ep), intent(in)    :: e(*), B(ldb,*)
            real(ep),    intent(inout) :: df(*)
            complex(ep), intent(inout) :: ef(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
            complex(ep), intent(out)   :: work(*)
            integer,     intent(out)   :: info
        end subroutine zptsvx_quad

        subroutine dtrsna_quad(job, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, &
                          S, SEP, mm, m, work, ldwork, iwork, info)
            import :: ep
            character, intent(in)    :: job, howmny
            logical,   intent(in)    :: sel(*)
            integer,   intent(in)    :: n, ldt, ldvl, ldvr, mm, ldwork
            real(ep),  intent(in)    :: T(ldt,*), VL(ldvl,*), VR(ldvr,*)
            real(ep),  intent(out)   :: S(*), SEP(*), work(ldwork,*)
            integer,   intent(out)   :: m, iwork(*), info
        end subroutine dtrsna_quad

        subroutine ztrsna_quad(job, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, &
                          S, SEP, mm, m, work, ldwork, rwork, info)
            import :: ep
            character,   intent(in)    :: job, howmny
            logical,     intent(in)    :: sel(*)
            integer,     intent(in)    :: n, ldt, ldvl, ldvr, mm, ldwork
            complex(ep), intent(in)    :: T(ldt,*), VL(ldvl,*), VR(ldvr,*)
            real(ep),    intent(out)   :: S(*), SEP(*), rwork(*)
            complex(ep), intent(out)   :: work(ldwork,*)
            integer,     intent(out)   :: m, info
        end subroutine ztrsna_quad

        ! ── 2-stage tridiag reductions ──────────────────────────────
        subroutine dsytrd_2stage_quad(vect, uplo, n, A, lda, D, E, tau, hous2, lhous2, work, lwork, info)
            import :: ep
            character, intent(in)    :: vect, uplo
            integer,   intent(in)    :: n, lda, lhous2, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: D(*), E(*), tau(*), hous2(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsytrd_2stage_quad

        subroutine zhetrd_2stage_quad(vect, uplo, n, A, lda, D, E, tau, hous2, lhous2, work, lwork, info)
            import :: ep
            character,   intent(in)    :: vect, uplo
            integer,     intent(in)    :: n, lda, lhous2, lwork
            complex(ep), intent(inout) :: A(lda,*)
            real(ep),    intent(out)   :: D(*), E(*)
            complex(ep), intent(out)   :: tau(*), hous2(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhetrd_2stage_quad

        subroutine dsytrd_sy2sb_quad(uplo, n, kd, A, lda, AB, ldab, tau, work, lwork, info)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, kd, lda, ldab, lwork
            real(ep),  intent(inout) :: A(lda,*)
            real(ep),  intent(out)   :: AB(ldab,*), tau(*), work(*)
            integer,   intent(out)   :: info
        end subroutine dsytrd_sy2sb_quad

        subroutine zhetrd_he2hb_quad(uplo, n, kd, A, lda, AB, ldab, tau, work, lwork, info)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, kd, lda, ldab, lwork
            complex(ep), intent(inout) :: A(lda,*)
            complex(ep), intent(out)   :: AB(ldab,*), tau(*), work(*)
            integer,     intent(out)   :: info
        end subroutine zhetrd_he2hb_quad

        subroutine dgedmd_quad(jobs, jobz, jobr, jobf, whtsvd, m, n, X, ldx, Y, ldy, &
                          nrnk, tol, k, reig, imeig, Z, ldz, res, B, ldb,        &
                          W, ldw, S, lds, work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobs, jobz, jobr, jobf
            integer,   intent(in)    :: whtsvd, m, n, ldx, ldy, nrnk, ldz, ldb, &
                                         ldw, lds, lwork, liwork
            real(ep),  intent(in)    :: tol
            real(ep),  intent(inout) :: X(ldx,*), Y(ldy,*)
            integer,   intent(out)   :: k, info
            real(ep),  intent(out)   :: reig(*), imeig(*), Z(ldz,*), res(*),    &
                                         B(ldb,*), W(ldw,*), S(lds,*), work(*)
            integer,   intent(out)   :: iwork(*)
        end subroutine dgedmd_quad

        subroutine zgedmd_quad(jobs, jobz, jobr, jobf, whtsvd, m, n, X, ldx, Y, ldy, &
                          nrnk, tol, k, eigs, Z, ldz, res, B, ldb,               &
                          W, ldw, S, lds, zwork, lzwork, rwork, lrwork,          &
                          iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobs, jobz, jobr, jobf
            integer,     intent(in)    :: whtsvd, m, n, ldx, ldy, nrnk, ldz, ldb, &
                                           ldw, lds, lzwork, lrwork, liwork
            real(ep),    intent(in)    :: tol
            complex(ep), intent(inout) :: X(ldx,*), Y(ldy,*)
            integer,     intent(out)   :: k, info
            complex(ep), intent(out)   :: eigs(*), Z(ldz,*), B(ldb,*), W(ldw,*), &
                                           S(lds,*), zwork(*)
            real(ep),    intent(out)   :: res(*), rwork(*)
            integer,     intent(out)   :: iwork(*)
        end subroutine zgedmd_quad

        subroutine dgedmdq_quad(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n,    &
                            F, ldf, X, ldx, Y, ldy, nrnk, tol, k, reig, imeig,  &
                            Z, ldz, res, B, ldb, V, ldv, S, lds,                 &
                            work, lwork, iwork, liwork, info)
            import :: ep
            character, intent(in)    :: jobs, jobz, jobr, jobq, jobt, jobf
            integer,   intent(in)    :: whtsvd, m, n, ldf, ldx, ldy, nrnk, ldz,  &
                                         ldb, ldv, lds, lwork, liwork
            real(ep),  intent(in)    :: tol
            real(ep),  intent(inout) :: F(ldf,*)
            integer,   intent(out)   :: k, info
            real(ep),  intent(out)   :: X(ldx,*), Y(ldy,*), Z(ldz,*), B(ldb,*),  &
                                         V(ldv,*), S(lds,*), reig(*), imeig(*),  &
                                         res(*), work(*)
            integer,   intent(out)   :: iwork(*)
        end subroutine dgedmdq_quad

        subroutine zgedmdq_quad(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n,    &
                            F, ldf, X, ldx, Y, ldy, nrnk, tol, k, eigs,          &
                            Z, ldz, res, B, ldb, V, ldv, S, lds,                 &
                            zwork, lzwork, work, lwork, iwork, liwork, info)
            import :: ep
            character,   intent(in)    :: jobs, jobz, jobr, jobq, jobt, jobf
            integer,     intent(in)    :: whtsvd, m, n, ldf, ldx, ldy, nrnk,     &
                                           ldz, ldb, ldv, lds, lzwork, lwork,    &
                                           liwork
            real(ep),    intent(in)    :: tol
            complex(ep), intent(inout) :: F(ldf,*)
            integer,     intent(out)   :: k, info
            complex(ep), intent(out)   :: X(ldx,*), Y(ldy,*), Z(ldz,*), B(ldb,*), &
                                           V(ldv,*), S(lds,*), eigs(*), zwork(*)
            real(ep),    intent(out)   :: res(*), work(*)
            integer,     intent(out)   :: iwork(*)
        end subroutine zgedmdq_quad

        ! ── Phase P17xx — extra-precise iterative-refinement xx family ─
        !
        ! These xx-suffix drivers (FACT='E' equilibrate-then-factor and
        ! refine, or RFSX refine-only) perform residual refinement in
        ! higher precision than the working precision via an XBLAS
        ! (head, tail) accumulator pair. They return per-RHS BERR plus
        ! a 2-D ERR_BNDS_NORM/COMP array of dimension
        ! (NRHS, N_ERR_BNDS) where N_ERR_BNDS = 3:
        !   col 1 = TRUST  (1.0 if bound is reliable, 0.0 if not),
        !   col 2 = NORMWISE/COMPONENTWISE error bound,
        !   col 3 = reciprocal of the condition number used.
        ! NPARAMS / PARAMS are tuning knobs; passing NPARAMS=0 takes
        ! all defaults (refinement on, threshold 10.0, no
        ! componentwise refinement) and the PARAMS array is not
        ! referenced.

        subroutine dgesvxx_quad(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, &
                           R, C, B, ldb, X, ldx, rcond, rpvgrw, berr,           &
                           n_err_bnds, err_bnds_norm, err_bnds_comp,            &
                           nparams, params, work, iwork, info)
            import :: ep
            character, intent(in)    :: fact, trans
            character, intent(inout) :: equed
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            real(ep),  intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*), R(*), C(*), params(*)
            integer,   intent(inout) :: ipiv(*)
            real(ep),  intent(out)   :: X(ldx,*), rcond, rpvgrw, berr(*),       &
                                        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgesvxx_quad

        subroutine zgesvxx_quad(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, &
                           R, C, B, ldb, X, ldx, rcond, rpvgrw, berr,           &
                           n_err_bnds, err_bnds_norm, err_bnds_comp,            &
                           nparams, params, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, trans
            character,   intent(inout) :: equed
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),    intent(inout) :: R(*), C(*), params(*)
            integer,     intent(inout) :: ipiv(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, rpvgrw, berr(*),               &
                                          err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zgesvxx_quad

        subroutine dgerfsx_quad(trans, equed, n, nrhs, A, lda, AF, ldaf, ipiv,       &
                           R, C, B, ldb, X, ldx, rcond, berr, n_err_bnds,       &
                           err_bnds_norm, err_bnds_comp, nparams, params,       &
                           work, iwork, info)
            import :: ep
            character, intent(in)    :: trans, equed
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*), R(*), C(*)
            real(ep),  intent(inout) :: X(ldx,*), params(*)
            real(ep),  intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*),  &
                                        err_bnds_comp(nrhs,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgerfsx_quad

        subroutine zgerfsx_quad(trans, equed, n, nrhs, A, lda, AF, ldaf, ipiv,       &
                           R, C, B, ldb, X, ldx, rcond, berr, n_err_bnds,       &
                           err_bnds_norm, err_bnds_comp, nparams, params,       &
                           work, rwork, info)
            import :: ep
            character,   intent(in)    :: trans, equed
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),    intent(in)    :: R(*), C(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(inout) :: params(*)
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
                                          err_bnds_comp(nrhs,*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zgerfsx_quad

        subroutine dgbsvxx_quad(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb,  &
                           ipiv, equed, R, C, B, ldb, X, ldx,                    &
                           rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm,       &
                           err_bnds_comp, nparams, params, work, iwork, info)
            import :: ep
            character, intent(in)    :: fact, trans
            character, intent(inout) :: equed
            integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx, n_err_bnds, nparams
            real(ep),  intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*), R(*), C(*), params(*)
            integer,   intent(inout) :: ipiv(*)
            real(ep),  intent(out)   :: X(ldx,*), rcond, rpvgrw, berr(*),       &
                                        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgbsvxx_quad

        subroutine zgbsvxx_quad(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb,  &
                           ipiv, equed, R, C, B, ldb, X, ldx,                    &
                           rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm,       &
                           err_bnds_comp, nparams, params, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, trans
            character,   intent(inout) :: equed
            integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx, n_err_bnds, nparams
            complex(ep), intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            real(ep),    intent(inout) :: R(*), C(*), params(*)
            integer,     intent(inout) :: ipiv(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, rpvgrw, berr(*),               &
                                          err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zgbsvxx_quad

        subroutine dgbrfsx_quad(trans, equed, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, &
                           ipiv, R, C, B, ldb, X, ldx, rcond, berr,             &
                           n_err_bnds, err_bnds_norm, err_bnds_comp,            &
                           nparams, params, work, iwork, info)
            import :: ep
            character, intent(in)    :: trans, equed
            integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx, n_err_bnds, nparams
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*), R(*), C(*)
            real(ep),  intent(inout) :: X(ldx,*), params(*)
            real(ep),  intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
                                        err_bnds_comp(nrhs,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dgbrfsx_quad

        subroutine zgbrfsx_quad(trans, equed, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, &
                           ipiv, R, C, B, ldb, X, ldx, rcond, berr,             &
                           n_err_bnds, err_bnds_norm, err_bnds_comp,            &
                           nparams, params, work, rwork, info)
            import :: ep
            character,   intent(in)    :: trans, equed
            integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx, n_err_bnds, nparams
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
            real(ep),    intent(in)    :: R(*), C(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(inout) :: params(*)
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
                                          err_bnds_comp(nrhs,*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zgbrfsx_quad

        subroutine dposvxx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S,     &
                           B, ldb, X, ldx, rcond, rpvgrw, berr, n_err_bnds,     &
                           err_bnds_norm, err_bnds_comp, nparams, params,       &
                           work, iwork, info)
            import :: ep
            character, intent(in)    :: fact, uplo
            character, intent(inout) :: equed
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            real(ep),  intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*), S(*), params(*)
            real(ep),  intent(out)   :: X(ldx,*), rcond, rpvgrw, berr(*),       &
                                        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dposvxx_quad

        subroutine zposvxx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S,     &
                           B, ldb, X, ldx, rcond, rpvgrw, berr, n_err_bnds,     &
                           err_bnds_norm, err_bnds_comp, nparams, params,       &
                           work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, uplo
            character,   intent(inout) :: equed
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),    intent(inout) :: S(*), params(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, rpvgrw, berr(*),               &
                                          err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zposvxx_quad

        subroutine dporfsx_quad(uplo, equed, n, nrhs, A, lda, AF, ldaf, S,           &
                           B, ldb, X, ldx, rcond, berr, n_err_bnds,             &
                           err_bnds_norm, err_bnds_comp, nparams, params,       &
                           work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo, equed
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), S(*), B(ldb,*)
            real(ep),  intent(inout) :: X(ldx,*), params(*)
            real(ep),  intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*),  &
                                        err_bnds_comp(nrhs,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dporfsx_quad

        subroutine zporfsx_quad(uplo, equed, n, nrhs, A, lda, AF, ldaf, S,           &
                           B, ldb, X, ldx, rcond, berr, n_err_bnds,             &
                           err_bnds_norm, err_bnds_comp, nparams, params,       &
                           work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo, equed
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),    intent(in)    :: S(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(inout) :: params(*)
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
                                          err_bnds_comp(nrhs,*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zporfsx_quad

        subroutine dsysvxx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, equed,  &
                           S, B, ldb, X, ldx, rcond, rpvgrw, berr,              &
                           n_err_bnds, err_bnds_norm, err_bnds_comp,            &
                           nparams, params, work, iwork, info)
            import :: ep
            character, intent(in)    :: fact, uplo
            character, intent(inout) :: equed
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            real(ep),  intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*), S(*), params(*)
            integer,   intent(inout) :: ipiv(*)
            real(ep),  intent(out)   :: X(ldx,*), rcond, rpvgrw, berr(*),       &
                                        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsysvxx_quad

        subroutine zsysvxx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, equed,  &
                           S, B, ldb, X, ldx, rcond, rpvgrw, berr,              &
                           n_err_bnds, err_bnds_norm, err_bnds_comp,            &
                           nparams, params, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, uplo
            character,   intent(inout) :: equed
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),    intent(inout) :: S(*), params(*)
            integer,     intent(inout) :: ipiv(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, rpvgrw, berr(*),               &
                                          err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zsysvxx_quad

        subroutine dsyrfsx_quad(uplo, equed, n, nrhs, A, lda, AF, ldaf, ipiv, S,     &
                           B, ldb, X, ldx, rcond, berr, n_err_bnds,             &
                           err_bnds_norm, err_bnds_comp, nparams, params,       &
                           work, iwork, info)
            import :: ep
            character, intent(in)    :: uplo, equed
            integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            integer,   intent(in)    :: ipiv(*)
            real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), S(*), B(ldb,*)
            real(ep),  intent(inout) :: X(ldx,*), params(*)
            real(ep),  intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
                                        err_bnds_comp(nrhs,*), work(*)
            integer,   intent(out)   :: iwork(*), info
        end subroutine dsyrfsx_quad

        subroutine zsyrfsx_quad(uplo, equed, n, nrhs, A, lda, AF, ldaf, ipiv, S,     &
                           B, ldb, X, ldx, rcond, berr, n_err_bnds,             &
                           err_bnds_norm, err_bnds_comp, nparams, params,       &
                           work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo, equed
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),    intent(in)    :: S(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(inout) :: params(*)
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
                                          err_bnds_comp(nrhs,*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zsyrfsx_quad

        subroutine zhesvxx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, equed,  &
                           S, B, ldb, X, ldx, rcond, rpvgrw, berr,              &
                           n_err_bnds, err_bnds_norm, err_bnds_comp,            &
                           nparams, params, work, rwork, info)
            import :: ep
            character,   intent(in)    :: fact, uplo
            character,   intent(inout) :: equed
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),    intent(inout) :: S(*), params(*)
            integer,     intent(inout) :: ipiv(*)
            complex(ep), intent(out)   :: X(ldx,*), work(*)
            real(ep),    intent(out)   :: rcond, rpvgrw, berr(*),               &
                                          err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zhesvxx_quad

        subroutine zherfsx_quad(uplo, equed, n, nrhs, A, lda, AF, ldaf, ipiv, S,     &
                           B, ldb, X, ldx, rcond, berr, n_err_bnds,             &
                           err_bnds_norm, err_bnds_comp, nparams, params,       &
                           work, rwork, info)
            import :: ep
            character,   intent(in)    :: uplo, equed
            integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
            integer,     intent(in)    :: ipiv(*)
            complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
            real(ep),    intent(in)    :: S(*)
            complex(ep), intent(inout) :: X(ldx,*)
            real(ep),    intent(inout) :: params(*)
            complex(ep), intent(out)   :: work(*)
            real(ep),    intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
                                          err_bnds_comp(nrhs,*), rwork(*)
            integer,     intent(out)   :: info
        end subroutine zherfsx_quad
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

    subroutine zgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        call zgesv_quad(n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine zgesv

    subroutine zgetrf(m, n, A, lda, ipiv, info)
        integer,     intent(in)    :: m, n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ipiv(*), info
        call zgetrf_quad(m, n, A, lda, ipiv, info)
    end subroutine zgetrf

    subroutine zgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(in)    :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zgetrs_quad(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine zgetrs

    subroutine zpotrf(uplo, n, A, lda, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: info
        call zpotrf_quad(uplo, n, A, lda, info)
    end subroutine zpotrf

    subroutine zpotrs(uplo, n, nrhs, A, lda, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zpotrs_quad(uplo, n, nrhs, A, lda, B, ldb, info)
    end subroutine zpotrs

    subroutine zungqr(m, n, k, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, k, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: tau(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zungqr_quad(m, n, k, A, lda, tau, work, lwork, info)
    end subroutine zungqr

    subroutine zgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, rwork, info)
        character,   intent(in)    :: jobu, jobvt
        integer,     intent(in)    :: m, n, lda, ldu, ldvt, lwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: s(*), rwork(*)
        complex(ep), intent(out)   :: U(ldu,*), VT(ldvt,*), work(*)
        integer,     intent(out)   :: info
        call zgesvd_quad(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, rwork, info)
    end subroutine zgesvd

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

    subroutine zlacpy(uplo, m, n, A, lda, B, ldb)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: m, n, lda, ldb
        complex(ep), intent(in)  :: A(lda,*)
        complex(ep), intent(out) :: B(ldb,*)
        call zlacpy_quad(uplo, m, n, A, lda, B, ldb)
    end subroutine zlacpy

    subroutine zlaset(uplo, m, n, alpha, beta, A, lda)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: m, n, lda
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(inout) :: A(lda,*)
        call zlaset_quad(uplo, m, n, alpha, beta, A, lda)
    end subroutine zlaset

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

    subroutine zgeqrf(m, n, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zgeqrf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine zgeqrf

    subroutine dsyev(jobz, uplo, n, A, lda, w, work, lwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: w(*), work(*)
        integer,   intent(out)   :: info
        call dsyev_quad(jobz, uplo, n, A, lda, w, work, lwork, info)
    end subroutine dsyev

    subroutine zheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zheev_quad(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
    end subroutine zheev

    subroutine dgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, info)
        character, intent(in)    :: jobu, jobvt
        integer,   intent(in)    :: m, n, lda, ldu, ldvt, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*), work(*)
        integer,   intent(out)   :: info
        call dgesvd_quad(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, info)
    end subroutine dgesvd

    function dlange(norm, m, n, A, lda, work) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: m, n, lda
        real(ep),  intent(in) :: A(lda,*)
        real(ep) :: work(*)
        real(ep) :: r
        r = dlange_quad(norm, m, n, A, lda, work)
    end function dlange

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

    subroutine dposv(uplo, n, nrhs, A, lda, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        integer,   intent(out)   :: info
        call dposv_quad(uplo, n, nrhs, A, lda, B, ldb, info)
    end subroutine dposv

    subroutine dsysv(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        integer,   intent(out)   :: ipiv(*), info
        real(ep),  intent(out)   :: work(*)
        call dsysv_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine dsysv

    subroutine dsytrf(uplo, n, A, lda, ipiv, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: ipiv(*), info
        real(ep),  intent(out)   :: work(*)
        call dsytrf_quad(uplo, n, A, lda, ipiv, work, lwork, info)
    end subroutine dsytrf

    subroutine dsytrs(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dsytrs_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine dsytrs

    subroutine dgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, nrhs, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dgels_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
    end subroutine dgels

    subroutine dsyevd(jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, lda, lwork, liwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: w(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsyevd_quad(jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info)
    end subroutine dsyevd

    subroutine dgeev(jobvl, jobvr, n, A, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
        character, intent(in)    :: jobvl, jobvr
        integer,   intent(in)    :: n, lda, ldvl, ldvr, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: wr(*), wi(*), vl(ldvl,*), vr(ldvr,*), work(*)
        integer,   intent(out)   :: info
        call dgeev_quad(jobvl, jobvr, n, A, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
    end subroutine dgeev

    subroutine dgesdd(jobz, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, iwork, info)
        character, intent(in)    :: jobz
        integer,   intent(in)    :: m, n, lda, ldu, ldvt, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: s(*), U(ldu,*), VT(ldvt,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgesdd_quad(jobz, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, iwork, info)
    end subroutine dgesdd

    subroutine dpotri(uplo, n, A, lda, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: info
        call dpotri_quad(uplo, n, A, lda, info)
    end subroutine dpotri

    subroutine dgetri(n, A, lda, ipiv, work, lwork, info)
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dgetri_quad(n, A, lda, ipiv, work, lwork, info)
    end subroutine dgetri

    subroutine dtrtri(uplo, diag, n, A, lda, info)
        character, intent(in)    :: uplo, diag
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: info
        call dtrtri_quad(uplo, diag, n, A, lda, info)
    end subroutine dtrtri

    subroutine dtrtrs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dtrtrs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
    end subroutine dtrtrs

    subroutine ztrtri(uplo, diag, n, A, lda, info)
        character,   intent(in)    :: uplo, diag
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: info
        call ztrtri_quad(uplo, diag, n, A, lda, info)
    end subroutine ztrtri

    subroutine ztrtrs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
        character,   intent(in)    :: uplo, trans, diag
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call ztrtrs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, info)
    end subroutine ztrtrs

    subroutine zposv(uplo, n, nrhs, A, lda, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: info
        call zposv_quad(uplo, n, nrhs, A, lda, B, ldb, info)
    end subroutine zposv

    subroutine zhesv(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zhesv_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine zhesv

    subroutine zheevd(jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, lda, lwork, lrwork, liwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: iwork(*), info
        call zheevd_quad(jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
    end subroutine zheevd

    subroutine zgeev(jobvl, jobvr, n, A, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
        character,   intent(in)    :: jobvl, jobvr
        integer,     intent(in)    :: n, lda, ldvl, ldvr, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: w(*), vl(ldvl,*), vr(ldvr,*), work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: info
        call zgeev_quad(jobvl, jobvr, n, A, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    end subroutine zgeev

    subroutine zgesdd(jobz, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, rwork, iwork, info)
        character,   intent(in)    :: jobz
        integer,     intent(in)    :: m, n, lda, ldu, ldvt, lwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: s(*), rwork(*)
        complex(ep), intent(out)   :: U(ldu,*), VT(ldvt,*), work(*)
        integer,     intent(out)   :: iwork(*), info
        call zgesdd_quad(jobz, m, n, A, lda, s, U, ldu, VT, ldvt, work, lwork, rwork, iwork, info)
    end subroutine zgesdd

    subroutine zgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgels_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
    end subroutine zgels

    subroutine dsyevr(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, &
                      work, lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, range, uplo
        integer,   intent(in)    :: n, lda, il, iu, ldz, lwork, liwork
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: m, isuppz(*), iwork(*), info
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        call dsyevr_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, &
                         work, lwork, iwork, liwork, info)
    end subroutine dsyevr

    subroutine zheevr(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, &
                      work, lwork, rwork, lrwork, iwork, liwork, info)
        character,   intent(in)    :: jobz, range, uplo
        integer,     intent(in)    :: n, lda, il, iu, ldz, lwork, lrwork, liwork
        real(ep),    intent(in)    :: vl, vu, abstol
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: m, isuppz(*), iwork(*), info
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: z(ldz,*), work(*)
        call zheevr_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, &
                         work, lwork, rwork, lrwork, iwork, liwork, info)
    end subroutine zheevr

    subroutine dsygv(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: itype, n, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: w(*), work(*)
        integer,   intent(out)   :: info
        call dsygv_quad(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, info)
    end subroutine dsygv

    subroutine zhegv(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, rwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: itype, n, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zhegv_quad(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, rwork, info)
    end subroutine zhegv

    subroutine dgehrd(n, ilo, ihi, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: n, ilo, ihi, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dgehrd_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
    end subroutine dgehrd

    subroutine dorghr(n, ilo, ihi, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: n, ilo, ihi, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: tau(*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dorghr_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
    end subroutine dorghr

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

    subroutine zgetri(n, A, lda, ipiv, work, lwork, info)
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgetri_quad(n, A, lda, ipiv, work, lwork, info)
    end subroutine zgetri

    subroutine zpotri(uplo, n, A, lda, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: info
        call zpotri_quad(uplo, n, A, lda, info)
    end subroutine zpotri

    subroutine zhetrf(uplo, n, A, lda, ipiv, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zhetrf_quad(uplo, n, A, lda, ipiv, work, lwork, info)
    end subroutine zhetrf

    subroutine zhetrs(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(in)    :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zhetrs_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine zhetrs

    subroutine zgehrd(n, ilo, ihi, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: n, ilo, ihi, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zgehrd_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
    end subroutine zgehrd

    subroutine zunghr(n, ilo, ihi, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: n, ilo, ihi, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: tau(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunghr_quad(n, ilo, ihi, A, lda, tau, work, lwork, info)
    end subroutine zunghr

    subroutine dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        integer,  intent(inout) :: jpvt(*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dgeqp3_quad(m, n, A, lda, jpvt, tau, work, lwork, info)
    end subroutine dgeqp3

    subroutine zgeqp3(m, n, A, lda, jpvt, tau, work, lwork, rwork, info)
        integer,     intent(in)    :: m, n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(inout) :: jpvt(*)
        complex(ep), intent(out)   :: tau(*), work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: info
        call zgeqp3_quad(m, n, A, lda, jpvt, tau, work, lwork, rwork, info)
    end subroutine zgeqp3

    subroutine dormqr(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormqr_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormqr

    subroutine zunmqr(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmqr_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmqr

    subroutine dgbtrf(m, n, kl, ku, AB, ldab, ipiv, info)
        integer,  intent(in)    :: m, n, kl, ku, ldab
        real(ep), intent(inout) :: AB(ldab,*)
        integer,  intent(out)   :: ipiv(*), info
        call dgbtrf_quad(m, n, kl, ku, AB, ldab, ipiv, info)
    end subroutine dgbtrf

    subroutine dgbtrs(trans, n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldb
        real(ep),  intent(in)    :: AB(ldab,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dgbtrs_quad(trans, n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
    end subroutine dgbtrs

    subroutine dgbsv(n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
        integer,  intent(in)    :: n, kl, ku, nrhs, ldab, ldb
        real(ep), intent(inout) :: AB(ldab,*), B(ldb,*)
        integer,  intent(out)   :: ipiv(*), info
        call dgbsv_quad(n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
    end subroutine dgbsv

    subroutine dpbtrf(uplo, n, kd, AB, ldab, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, kd, ldab
        real(ep),  intent(inout) :: AB(ldab,*)
        integer,   intent(out)   :: info
        call dpbtrf_quad(uplo, n, kd, AB, ldab, info)
    end subroutine dpbtrf

    subroutine dpbtrs(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, kd, nrhs, ldab, ldb
        real(ep),  intent(in)    :: AB(ldab,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dpbtrs_quad(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
    end subroutine dpbtrs

    subroutine dpbsv(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, kd, nrhs, ldab, ldb
        real(ep),  intent(inout) :: AB(ldab,*), B(ldb,*)
        integer,   intent(out)   :: info
        call dpbsv_quad(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
    end subroutine dpbsv

    subroutine dgttrf(n, dl, d, du, du2, ipiv, info)
        integer,  intent(in)    :: n
        real(ep), intent(inout) :: dl(*), d(*), du(*)
        real(ep), intent(out)   :: du2(*)
        integer,  intent(out)   :: ipiv(*), info
        call dgttrf_quad(n, dl, d, du, du2, ipiv, info)
    end subroutine dgttrf

    subroutine dgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, B, ldb, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: n, nrhs, ldb
        real(ep),  intent(in)    :: dl(*), d(*), du(*), du2(*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dgttrs_quad(trans, n, nrhs, dl, d, du, du2, ipiv, B, ldb, info)
    end subroutine dgttrs

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

    subroutine dsyevx(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, &
                      lwork, iwork, ifail, info)
        character, intent(in)    :: jobz, range, uplo
        integer,   intent(in)    :: n, lda, il, iu, ldz, lwork
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: m, iwork(*), ifail(*), info
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        call dsyevx_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, &
                         lwork, iwork, ifail, info)
    end subroutine dsyevx

    subroutine zheevx(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, &
                      lwork, rwork, iwork, ifail, info)
        character,   intent(in)    :: jobz, range, uplo
        integer,     intent(in)    :: n, lda, il, iu, ldz, lwork
        real(ep),    intent(in)    :: vl, vu, abstol
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: m, iwork(*), ifail(*), info
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: z(ldz,*), work(*)
        call zheevx_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, &
                         lwork, rwork, iwork, ifail, info)
    end subroutine zheevx

    subroutine dstev(jobz, n, d, e, z, ldz, work, info)
        character, intent(in)    :: jobz
        integer,   intent(in)    :: n, ldz
        real(ep),  intent(inout) :: d(*), e(*)
        real(ep),  intent(out)   :: z(ldz,*), work(*)
        integer,   intent(out)   :: info
        call dstev_quad(jobz, n, d, e, z, ldz, work, info)
    end subroutine dstev

    subroutine dstevd(jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info)
        character, intent(in)    :: jobz
        integer,   intent(in)    :: n, ldz, lwork, liwork
        real(ep),  intent(inout) :: d(*), e(*)
        real(ep),  intent(out)   :: z(ldz,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dstevd_quad(jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info)
    end subroutine dstevd

    subroutine dstevx(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, &
                      info)
        character, intent(in)    :: jobz, range
        integer,   intent(in)    :: n, il, iu, ldz
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: d(*), e(*)
        integer,   intent(out)   :: m, iwork(*), ifail(*), info
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        call dstevx_quad(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, &
                         ifail, info)
    end subroutine dstevx

    subroutine dstevr(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, &
                      lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, range
        integer,   intent(in)    :: n, il, iu, ldz, lwork, liwork
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: d(*), e(*)
        integer,   intent(out)   :: m, isuppz(*), iwork(*), info
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        call dstevr_quad(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, &
                         lwork, iwork, liwork, info)
    end subroutine dstevr

    subroutine dsbev(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, kd, ldab, ldz
        real(ep),  intent(inout) :: AB(ldab,*)
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        integer,   intent(out)   :: info
        call dsbev_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, info)
    end subroutine dsbev

    subroutine dsbevd(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, kd, ldab, ldz, lwork, liwork
        real(ep),  intent(inout) :: AB(ldab,*)
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsbevd_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, iwork, liwork, info)
    end subroutine dsbevd

    subroutine dsbevx(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, vl, vu, il, iu, abstol, m, w, z, &
                      ldz, work, iwork, ifail, info)
        character, intent(in)    :: jobz, range, uplo
        integer,   intent(in)    :: n, kd, ldab, ldq, il, iu, ldz
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: AB(ldab,*)
        real(ep),  intent(out)   :: Q(ldq,*), w(*), z(ldz,*), work(*)
        integer,   intent(out)   :: m, iwork(*), ifail(*), info
        call dsbevx_quad(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, vl, vu, il, iu, abstol, m, w, z, &
                         ldz, work, iwork, ifail, info)
    end subroutine dsbevx

    subroutine zhbev(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, rwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, kd, ldab, ldz
        complex(ep), intent(inout) :: AB(ldab,*)
        complex(ep), intent(out)   :: z(ldz,*), work(*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        integer,     intent(out)   :: info
        call zhbev_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, rwork, info)
    end subroutine zhbev

    subroutine zhbevd(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, &
                      liwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, kd, ldab, ldz, lwork, lrwork, liwork
        complex(ep), intent(inout) :: AB(ldab,*)
        complex(ep), intent(out)   :: z(ldz,*), work(*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        integer,     intent(out)   :: iwork(*), info
        call zhbevd_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork, &
                         liwork, info)
    end subroutine zhbevd

    subroutine zhbevx(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, vl, vu, il, iu, abstol, m, w, z, &
                      ldz, work, rwork, iwork, ifail, info)
        character,   intent(in)    :: jobz, range, uplo
        integer,     intent(in)    :: n, kd, ldab, ldq, il, iu, ldz
        real(ep),    intent(in)    :: vl, vu, abstol
        complex(ep), intent(inout) :: AB(ldab,*)
        complex(ep), intent(out)   :: Q(ldq,*), z(ldz,*), work(*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        integer,     intent(out)   :: m, iwork(*), ifail(*), info
        call zhbevx_quad(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, vl, vu, il, iu, abstol, m, w, z, &
                         ldz, work, rwork, iwork, ifail, info)
    end subroutine zhbevx

    subroutine dspev(jobz, uplo, n, AP, w, z, ldz, work, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, ldz
        real(ep),  intent(inout) :: AP(*)
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        integer,   intent(out)   :: info
        call dspev_quad(jobz, uplo, n, AP, w, z, ldz, work, info)
    end subroutine dspev

    subroutine dspevd(jobz, uplo, n, AP, w, z, ldz, work, lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, ldz, lwork, liwork
        real(ep),  intent(inout) :: AP(*)
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dspevd_quad(jobz, uplo, n, AP, w, z, ldz, work, lwork, iwork, liwork, info)
    end subroutine dspevd

    subroutine dspevx(jobz, range, uplo, n, AP, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, &
                      ifail, info)
        character, intent(in)    :: jobz, range, uplo
        integer,   intent(in)    :: n, il, iu, ldz
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: AP(*)
        integer,   intent(out)   :: m, iwork(*), ifail(*), info
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        call dspevx_quad(jobz, range, uplo, n, AP, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, &
                         ifail, info)
    end subroutine dspevx

    subroutine zhpev(jobz, uplo, n, AP, w, z, ldz, work, rwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, ldz
        complex(ep), intent(inout) :: AP(*)
        complex(ep), intent(out)   :: z(ldz,*), work(*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        integer,     intent(out)   :: info
        call zhpev_quad(jobz, uplo, n, AP, w, z, ldz, work, rwork, info)
    end subroutine zhpev

    subroutine zhpevd(jobz, uplo, n, AP, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, ldz, lwork, lrwork, liwork
        complex(ep), intent(inout) :: AP(*)
        complex(ep), intent(out)   :: z(ldz,*), work(*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        integer,     intent(out)   :: iwork(*), info
        call zhpevd_quad(jobz, uplo, n, AP, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, &
                         info)
    end subroutine zhpevd

    subroutine zhpevx(jobz, range, uplo, n, AP, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, &
                      iwork, ifail, info)
        character,   intent(in)    :: jobz, range, uplo
        integer,     intent(in)    :: n, il, iu, ldz
        real(ep),    intent(in)    :: vl, vu, abstol
        complex(ep), intent(inout) :: AP(*)
        complex(ep), intent(out)   :: z(ldz,*), work(*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        integer,     intent(out)   :: m, iwork(*), ifail(*), info
        call zhpevx_quad(jobz, range, uplo, n, AP, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, &
                         iwork, ifail, info)
    end subroutine zhpevx

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

    subroutine dsptrf(uplo, n, AP, ipiv, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n
        real(ep),  intent(inout) :: AP(*)
        integer,   intent(out)   :: ipiv(*), info
        call dsptrf_quad(uplo, n, AP, ipiv, info)
    end subroutine dsptrf

    subroutine dsptrs(uplo, n, nrhs, AP, ipiv, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, ldb
        real(ep),  intent(in)    :: AP(*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dsptrs_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
    end subroutine dsptrs

    subroutine dspsv(uplo, n, nrhs, AP, ipiv, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, ldb
        real(ep),  intent(inout) :: AP(*), B(ldb,*)
        integer,   intent(out)   :: ipiv(*), info
        call dspsv_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
    end subroutine dspsv

    subroutine dsptri(uplo, n, AP, ipiv, work, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n
        real(ep),  intent(inout) :: AP(*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dsptri_quad(uplo, n, AP, ipiv, work, info)
    end subroutine dsptri

    subroutine dpptrf(uplo, n, AP, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n
        real(ep),  intent(inout) :: AP(*)
        integer,   intent(out)   :: info
        call dpptrf_quad(uplo, n, AP, info)
    end subroutine dpptrf

    subroutine dpptrs(uplo, n, nrhs, AP, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, ldb
        real(ep),  intent(in)    :: AP(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dpptrs_quad(uplo, n, nrhs, AP, B, ldb, info)
    end subroutine dpptrs

    subroutine dppsv(uplo, n, nrhs, AP, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, ldb
        real(ep),  intent(inout) :: AP(*), B(ldb,*)
        integer,   intent(out)   :: info
        call dppsv_quad(uplo, n, nrhs, AP, B, ldb, info)
    end subroutine dppsv

    subroutine dpptri(uplo, n, AP, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n
        real(ep),  intent(inout) :: AP(*)
        integer,   intent(out)   :: info
        call dpptri_quad(uplo, n, AP, info)
    end subroutine dpptri

    subroutine zhptrf(uplo, n, AP, ipiv, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: AP(*)
        integer,     intent(out)   :: ipiv(*), info
        call zhptrf_quad(uplo, n, AP, ipiv, info)
    end subroutine zhptrf

    subroutine zhptrs(uplo, n, nrhs, AP, ipiv, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, ldb
        complex(ep), intent(in)    :: AP(*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zhptrs_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
    end subroutine zhptrs

    subroutine zhpsv(uplo, n, nrhs, AP, ipiv, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, ldb
        complex(ep), intent(inout) :: AP(*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        call zhpsv_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
    end subroutine zhpsv

    subroutine zhptri(uplo, n, AP, ipiv, work, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: AP(*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zhptri_quad(uplo, n, AP, ipiv, work, info)
    end subroutine zhptri

    subroutine zsptrf(uplo, n, AP, ipiv, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: AP(*)
        integer,     intent(out)   :: ipiv(*), info
        call zsptrf_quad(uplo, n, AP, ipiv, info)
    end subroutine zsptrf

    subroutine zsptrs(uplo, n, nrhs, AP, ipiv, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, ldb
        complex(ep), intent(in)    :: AP(*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zsptrs_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
    end subroutine zsptrs

    subroutine zspsv(uplo, n, nrhs, AP, ipiv, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, ldb
        complex(ep), intent(inout) :: AP(*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        call zspsv_quad(uplo, n, nrhs, AP, ipiv, B, ldb, info)
    end subroutine zspsv

    subroutine zsptri(uplo, n, AP, ipiv, work, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: AP(*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zsptri_quad(uplo, n, AP, ipiv, work, info)
    end subroutine zsptri

    subroutine zpptrf(uplo, n, AP, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: AP(*)
        integer,     intent(out)   :: info
        call zpptrf_quad(uplo, n, AP, info)
    end subroutine zpptrf

    subroutine zpptrs(uplo, n, nrhs, AP, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, ldb
        complex(ep), intent(in)    :: AP(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zpptrs_quad(uplo, n, nrhs, AP, B, ldb, info)
    end subroutine zpptrs

    subroutine zppsv(uplo, n, nrhs, AP, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, ldb
        complex(ep), intent(inout) :: AP(*), B(ldb,*)
        integer,     intent(out)   :: info
        call zppsv_quad(uplo, n, nrhs, AP, B, ldb, info)
    end subroutine zppsv

    subroutine zpptri(uplo, n, AP, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: AP(*)
        integer,     intent(out)   :: info
        call zpptri_quad(uplo, n, AP, info)
    end subroutine zpptri

    subroutine zgbtrf(m, n, kl, ku, AB, ldab, ipiv, info)
        integer,     intent(in)    :: m, n, kl, ku, ldab
        complex(ep), intent(inout) :: AB(ldab,*)
        integer,     intent(out)   :: ipiv(*), info
        call zgbtrf_quad(m, n, kl, ku, AB, ldab, ipiv, info)
    end subroutine zgbtrf

    subroutine zgbtrs(trans, n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldb
        complex(ep), intent(in)    :: AB(ldab,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zgbtrs_quad(trans, n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
    end subroutine zgbtrs

    subroutine zgbsv(n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
        integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldb
        complex(ep), intent(inout) :: AB(ldab,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        call zgbsv_quad(n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, info)
    end subroutine zgbsv

    subroutine zpbtrf(uplo, n, kd, AB, ldab, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, kd, ldab
        complex(ep), intent(inout) :: AB(ldab,*)
        integer,     intent(out)   :: info
        call zpbtrf_quad(uplo, n, kd, AB, ldab, info)
    end subroutine zpbtrf

    subroutine zpbtrs(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, kd, nrhs, ldab, ldb
        complex(ep), intent(in)    :: AB(ldab,*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zpbtrs_quad(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
    end subroutine zpbtrs

    subroutine zpbsv(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, kd, nrhs, ldab, ldb
        complex(ep), intent(inout) :: AB(ldab,*), B(ldb,*)
        integer,     intent(out)   :: info
        call zpbsv_quad(uplo, n, kd, nrhs, AB, ldab, B, ldb, info)
    end subroutine zpbsv

    subroutine zgttrf(n, dl, d, du, du2, ipiv, info)
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: dl(*), d(*), du(*)
        complex(ep), intent(out)   :: du2(*)
        integer,     intent(out)   :: ipiv(*), info
        call zgttrf_quad(n, dl, d, du, du2, ipiv, info)
    end subroutine zgttrf

    subroutine zgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, B, ldb, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: n, nrhs, ldb
        complex(ep), intent(in)    :: dl(*), d(*), du(*), du2(*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zgttrs_quad(trans, n, nrhs, dl, d, du, du2, ipiv, B, ldb, info)
    end subroutine zgttrs

    subroutine zgtsv(n, nrhs, dl, d, du, B, ldb, info)
        integer,     intent(in)    :: n, nrhs, ldb
        complex(ep), intent(inout) :: dl(*), d(*), du(*), B(ldb,*)
        integer,     intent(out)   :: info
        call zgtsv_quad(n, nrhs, dl, d, du, B, ldb, info)
    end subroutine zgtsv

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

    subroutine zptsv(n, nrhs, d, e, B, ldb, info)
        integer,     intent(in)    :: n, nrhs, ldb
        real(ep),    intent(inout) :: d(*)
        complex(ep), intent(inout) :: e(*), B(ldb,*)
        integer,     intent(out)   :: info
        call zptsv_quad(n, nrhs, d, e, B, ldb, info)
    end subroutine zptsv

    subroutine dtbtrs(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, info)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, kd, nrhs, ldab, ldb
        real(ep),  intent(in)    :: AB(ldab,*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dtbtrs_quad(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, info)
    end subroutine dtbtrs

    subroutine ztbtrs(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, info)
        character,   intent(in)    :: uplo, trans, diag
        integer,     intent(in)    :: n, kd, nrhs, ldab, ldb
        complex(ep), intent(in)    :: AB(ldab,*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call ztbtrs_quad(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, info)
    end subroutine ztbtrs

    subroutine dtptrs(uplo, trans, diag, n, nrhs, AP, B, ldb, info)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, nrhs, ldb
        real(ep),  intent(in)    :: AP(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dtptrs_quad(uplo, trans, diag, n, nrhs, AP, B, ldb, info)
    end subroutine dtptrs

    subroutine ztptrs(uplo, trans, diag, n, nrhs, AP, B, ldb, info)
        character,   intent(in)    :: uplo, trans, diag
        integer,     intent(in)    :: n, nrhs, ldb
        complex(ep), intent(in)    :: AP(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call ztptrs_quad(uplo, trans, diag, n, nrhs, AP, B, ldb, info)
    end subroutine ztptrs

    subroutine dtptri(uplo, diag, n, AP, info)
        character, intent(in)    :: uplo, diag
        integer,   intent(in)    :: n
        real(ep),  intent(inout) :: AP(*)
        integer,   intent(out)   :: info
        call dtptri_quad(uplo, diag, n, AP, info)
    end subroutine dtptri

    subroutine ztptri(uplo, diag, n, AP, info)
        character,   intent(in)    :: uplo, diag
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: AP(*)
        integer,     intent(out)   :: info
        call ztptri_quad(uplo, diag, n, AP, info)
    end subroutine ztptri

    subroutine dgecon(norm, n, A, lda, anorm, rcond, work, iwork, info)
        character, intent(in)  :: norm
        integer,   intent(in)  :: n, lda
        real(ep),  intent(in)  :: A(lda,*), anorm
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dgecon_quad(norm, n, A, lda, anorm, rcond, work, iwork, info)
    end subroutine dgecon

    subroutine zgecon(norm, n, A, lda, anorm, rcond, work, rwork, info)
        character,   intent(in)  :: norm
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond, rwork(*)
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zgecon_quad(norm, n, A, lda, anorm, rcond, work, rwork, info)
    end subroutine zgecon

    subroutine dpocon(uplo, n, A, lda, anorm, rcond, work, iwork, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n, lda
        real(ep),  intent(in)  :: A(lda,*), anorm
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dpocon_quad(uplo, n, A, lda, anorm, rcond, work, iwork, info)
    end subroutine dpocon

    subroutine zpocon(uplo, n, A, lda, anorm, rcond, work, rwork, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond, rwork(*)
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zpocon_quad(uplo, n, A, lda, anorm, rcond, work, rwork, info)
    end subroutine zpocon

    subroutine dgbcon(norm, n, kl, ku, AB, ldab, ipiv, anorm, rcond, work, iwork, info)
        character, intent(in)  :: norm
        integer,   intent(in)  :: n, kl, ku, ldab
        real(ep),  intent(in)  :: AB(ldab,*), anorm
        integer,   intent(in)  :: ipiv(*)
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dgbcon_quad(norm, n, kl, ku, AB, ldab, ipiv, anorm, rcond, work, iwork, info)
    end subroutine dgbcon

    subroutine zgbcon(norm, n, kl, ku, AB, ldab, ipiv, anorm, rcond, work, rwork, info)
        character,   intent(in)  :: norm
        integer,     intent(in)  :: n, kl, ku, ldab
        complex(ep), intent(in)  :: AB(ldab,*)
        integer,     intent(in)  :: ipiv(*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond, rwork(*)
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zgbcon_quad(norm, n, kl, ku, AB, ldab, ipiv, anorm, rcond, work, rwork, info)
    end subroutine zgbcon

    subroutine dgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info)
        character, intent(in)  :: norm
        integer,   intent(in)  :: n
        real(ep),  intent(in)  :: dl(*), d(*), du(*), du2(*), anorm
        integer,   intent(in)  :: ipiv(*)
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dgtcon_quad(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info)
    end subroutine dgtcon

    subroutine zgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info)
        character,   intent(in)  :: norm
        integer,     intent(in)  :: n
        complex(ep), intent(in)  :: dl(*), d(*), du(*), du2(*)
        integer,     intent(in)  :: ipiv(*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zgtcon_quad(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info)
    end subroutine zgtcon

    subroutine dsycon(uplo, n, A, lda, ipiv, anorm, rcond, work, iwork, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n, lda
        real(ep),  intent(in)  :: A(lda,*), anorm
        integer,   intent(in)  :: ipiv(*)
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dsycon_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, iwork, info)
    end subroutine dsycon

    subroutine zhecon(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        integer,     intent(in)  :: ipiv(*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zhecon_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
    end subroutine zhecon

    subroutine dpbcon(uplo, n, kd, AB, ldab, anorm, rcond, work, iwork, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n, kd, ldab
        real(ep),  intent(in)  :: AB(ldab,*), anorm
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dpbcon_quad(uplo, n, kd, AB, ldab, anorm, rcond, work, iwork, info)
    end subroutine dpbcon

    subroutine zpbcon(uplo, n, kd, AB, ldab, anorm, rcond, work, rwork, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, kd, ldab
        complex(ep), intent(in)  :: AB(ldab,*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond, rwork(*)
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zpbcon_quad(uplo, n, kd, AB, ldab, anorm, rcond, work, rwork, info)
    end subroutine zpbcon

    subroutine dppcon(uplo, n, AP, anorm, rcond, work, iwork, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n
        real(ep),  intent(in)  :: AP(*), anorm
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dppcon_quad(uplo, n, AP, anorm, rcond, work, iwork, info)
    end subroutine dppcon

    subroutine zppcon(uplo, n, AP, anorm, rcond, work, rwork, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n
        complex(ep), intent(in)  :: AP(*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond, rwork(*)
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zppcon_quad(uplo, n, AP, anorm, rcond, work, rwork, info)
    end subroutine zppcon

    subroutine dptcon(n, d, e, anorm, rcond, work, info)
        integer,  intent(in)  :: n
        real(ep), intent(in)  :: d(*), e(*), anorm
        real(ep), intent(out) :: rcond, work(*)
        integer,  intent(out) :: info
        call dptcon_quad(n, d, e, anorm, rcond, work, info)
    end subroutine dptcon

    subroutine zptcon(n, d, e, anorm, rcond, rwork, info)
        integer,     intent(in)  :: n
        real(ep),    intent(in)  :: d(*), anorm
        complex(ep), intent(in)  :: e(*)
        real(ep),    intent(out) :: rcond, rwork(*)
        integer,     intent(out) :: info
        call zptcon_quad(n, d, e, anorm, rcond, rwork, info)
    end subroutine zptcon

    subroutine dspcon(uplo, n, AP, ipiv, anorm, rcond, work, iwork, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n
        real(ep),  intent(in)  :: AP(*), anorm
        integer,   intent(in)  :: ipiv(*)
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dspcon_quad(uplo, n, AP, ipiv, anorm, rcond, work, iwork, info)
    end subroutine dspcon

    subroutine zhpcon(uplo, n, AP, ipiv, anorm, rcond, work, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n
        complex(ep), intent(in)  :: AP(*)
        integer,     intent(in)  :: ipiv(*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zhpcon_quad(uplo, n, AP, ipiv, anorm, rcond, work, info)
    end subroutine zhpcon

    subroutine zspcon(uplo, n, AP, ipiv, anorm, rcond, work, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n
        complex(ep), intent(in)  :: AP(*)
        integer,     intent(in)  :: ipiv(*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zspcon_quad(uplo, n, AP, ipiv, anorm, rcond, work, info)
    end subroutine zspcon

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
        real(ep),    intent(out) :: rcond, rwork(*)
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call ztrcon_quad(norm, uplo, diag, n, A, lda, rcond, work, rwork, info)
    end subroutine ztrcon

    subroutine dtpcon(norm, uplo, diag, n, AP, rcond, work, iwork, info)
        character, intent(in)  :: norm, uplo, diag
        integer,   intent(in)  :: n
        real(ep),  intent(in)  :: AP(*)
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dtpcon_quad(norm, uplo, diag, n, AP, rcond, work, iwork, info)
    end subroutine dtpcon

    subroutine ztpcon(norm, uplo, diag, n, AP, rcond, work, rwork, info)
        character,   intent(in)  :: norm, uplo, diag
        integer,     intent(in)  :: n
        complex(ep), intent(in)  :: AP(*)
        real(ep),    intent(out) :: rcond, rwork(*)
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call ztpcon_quad(norm, uplo, diag, n, AP, rcond, work, rwork, info)
    end subroutine ztpcon

    subroutine dtbcon(norm, uplo, diag, n, kd, AB, ldab, rcond, work, iwork, info)
        character, intent(in)  :: norm, uplo, diag
        integer,   intent(in)  :: n, kd, ldab
        real(ep),  intent(in)  :: AB(ldab,*)
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dtbcon_quad(norm, uplo, diag, n, kd, AB, ldab, rcond, work, iwork, info)
    end subroutine dtbcon

    subroutine ztbcon(norm, uplo, diag, n, kd, AB, ldab, rcond, work, rwork, info)
        character,   intent(in)  :: norm, uplo, diag
        integer,     intent(in)  :: n, kd, ldab
        complex(ep), intent(in)  :: AB(ldab,*)
        real(ep),    intent(out) :: rcond, rwork(*)
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call ztbcon_quad(norm, uplo, diag, n, kd, AB, ldab, rcond, work, rwork, info)
    end subroutine ztbcon

    subroutine dgerfs(trans, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, ferr, berr, work, &
                      iwork, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
        real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: X(ldx,*)
        real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgerfs_quad(trans, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, ferr, berr, work, &
                         iwork, info)
    end subroutine dgerfs

    subroutine zgerfs(trans, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, ferr, berr, work, &
                      rwork, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
        complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgerfs_quad(trans, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, ferr, berr, work, &
                         rwork, info)
    end subroutine zgerfs

    subroutine dporfs(uplo, n, nrhs, A, lda, AF, ldaf, B, ldb, X, ldx, ferr, berr, work, iwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
        real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
        real(ep),  intent(inout) :: X(ldx,*)
        real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dporfs_quad(uplo, n, nrhs, A, lda, AF, ldaf, B, ldb, X, ldx, ferr, berr, work, iwork, &
                         info)
    end subroutine dporfs

    subroutine zporfs(uplo, n, nrhs, A, lda, AF, ldaf, B, ldb, X, ldx, ferr, berr, work, rwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
        complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zporfs_quad(uplo, n, nrhs, A, lda, AF, ldaf, B, ldb, X, ldx, ferr, berr, work, rwork, &
                         info)
    end subroutine zporfs

    subroutine dgbrfs(trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, B, ldb, X, ldx, ferr, berr, &
                      work, iwork, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx
        real(ep),  intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: X(ldx,*)
        real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgbrfs_quad(trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, B, ldb, X, ldx, ferr, &
                         berr, work, iwork, info)
    end subroutine dgbrfs

    subroutine zgbrfs(trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, B, ldb, X, ldx, ferr, berr, &
                      work, rwork, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx
        complex(ep), intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgbrfs_quad(trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, B, ldb, X, ldx, ferr, &
                         berr, work, rwork, info)
    end subroutine zgbrfs

    subroutine dpbrfs(uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, B, ldb, X, ldx, ferr, berr, work, &
                      iwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, kd, nrhs, ldab, ldafb, ldb, ldx
        real(ep),  intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
        real(ep),  intent(inout) :: X(ldx,*)
        real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dpbrfs_quad(uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, B, ldb, X, ldx, ferr, berr, work, &
                         iwork, info)
    end subroutine dpbrfs

    subroutine zpbrfs(uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, B, ldb, X, ldx, ferr, berr, work, &
                      rwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, kd, nrhs, ldab, ldafb, ldb, ldx
        complex(ep), intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zpbrfs_quad(uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, B, ldb, X, ldx, ferr, berr, work, &
                         rwork, info)
    end subroutine zpbrfs

    subroutine dgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, B, ldb, X, ldx, ferr, berr, &
                      work, iwork, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: n, nrhs, ldb, ldx
        real(ep),  intent(in)    :: dl(*), d(*), du(*), dlf(*), df(*), duf(*), du2(*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(in)    :: B(ldb,*)
        real(ep),  intent(inout) :: X(ldx,*)
        real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgtrfs_quad(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, B, ldb, X, ldx, ferr, &
                         berr, work, iwork, info)
    end subroutine dgtrfs

    subroutine zgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, B, ldb, X, ldx, ferr, berr, &
                      work, rwork, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: n, nrhs, ldb, ldx
        complex(ep), intent(in)    :: dl(*), d(*), du(*), dlf(*), df(*), duf(*), du2(*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(in)    :: B(ldb,*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgtrfs_quad(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, B, ldb, X, ldx, ferr, &
                         berr, work, rwork, info)
    end subroutine zgtrfs

    subroutine dptrfs(n, nrhs, d, e, df, ef, B, ldb, X, ldx, ferr, berr, work, info)
        integer,  intent(in)    :: n, nrhs, ldb, ldx
        real(ep), intent(in)    :: d(*), e(*), df(*), ef(*), B(ldb,*)
        real(ep), intent(inout) :: X(ldx,*)
        real(ep), intent(out)   :: ferr(*), berr(*), work(*)
        integer,  intent(out)   :: info
        call dptrfs_quad(n, nrhs, d, e, df, ef, B, ldb, X, ldx, ferr, berr, work, info)
    end subroutine dptrfs

    subroutine zptrfs(uplo, n, nrhs, d, e, df, ef, B, ldb, X, ldx, ferr, berr, work, rwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, ldb, ldx
        real(ep),    intent(in)    :: d(*), df(*)
        complex(ep), intent(in)    :: e(*), ef(*), B(ldb,*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zptrfs_quad(uplo, n, nrhs, d, e, df, ef, B, ldb, X, ldx, ferr, berr, work, rwork, info)
    end subroutine zptrfs

    subroutine dsyrfs(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, ferr, berr, work, iwork, &
                      info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
        real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: X(ldx,*)
        real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsyrfs_quad(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, ferr, berr, work, &
                         iwork, info)
    end subroutine dsyrfs

    subroutine zsyrfs(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, ferr, berr, work, rwork, &
                      info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
        complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zsyrfs_quad(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, ferr, berr, work, &
                         rwork, info)
    end subroutine zsyrfs

    subroutine zherfs(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, ferr, berr, work, rwork, &
                      info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
        complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zherfs_quad(uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, ferr, berr, work, &
                         rwork, info)
    end subroutine zherfs

    subroutine dsprfs(uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, ferr, berr, work, iwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, ldb, ldx
        real(ep),  intent(in)    :: AP(*), AFP(*), B(ldb,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: X(ldx,*)
        real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsprfs_quad(uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, ferr, berr, work, iwork, info)
    end subroutine dsprfs

    subroutine zsprfs(uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, ferr, berr, work, rwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, ldb, ldx
        complex(ep), intent(in)    :: AP(*), AFP(*), B(ldb,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zsprfs_quad(uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, ferr, berr, work, rwork, info)
    end subroutine zsprfs

    subroutine dpprfs(uplo, n, nrhs, AP, AFP, B, ldb, X, ldx, ferr, berr, work, iwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, ldb, ldx
        real(ep),  intent(in)    :: AP(*), AFP(*), B(ldb,*)
        real(ep),  intent(inout) :: X(ldx,*)
        real(ep),  intent(out)   :: ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dpprfs_quad(uplo, n, nrhs, AP, AFP, B, ldb, X, ldx, ferr, berr, work, iwork, info)
    end subroutine dpprfs

    subroutine zpprfs(uplo, n, nrhs, AP, AFP, B, ldb, X, ldx, ferr, berr, work, rwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, ldb, ldx
        complex(ep), intent(in)    :: AP(*), AFP(*), B(ldb,*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zpprfs_quad(uplo, n, nrhs, AP, AFP, B, ldb, X, ldx, ferr, berr, work, rwork, info)
    end subroutine zpprfs

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

    subroutine ztrrfs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, X, ldx, ferr, berr, work, rwork, &
                      info)
        character,   intent(in)  :: uplo, trans, diag
        integer,     intent(in)  :: n, nrhs, lda, ldb, ldx
        complex(ep), intent(in)  :: A(lda,*), B(ldb,*), X(ldx,*)
        real(ep),    intent(out) :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call ztrrfs_quad(uplo, trans, diag, n, nrhs, A, lda, B, ldb, X, ldx, ferr, berr, work, rwork, &
                         info)
    end subroutine ztrrfs

    subroutine dtprfs(uplo, trans, diag, n, nrhs, AP, B, ldb, X, ldx, ferr, berr, work, iwork, info)
        character, intent(in)  :: uplo, trans, diag
        integer,   intent(in)  :: n, nrhs, ldb, ldx
        real(ep),  intent(in)  :: AP(*), B(ldb,*), X(ldx,*)
        real(ep),  intent(out) :: ferr(*), berr(*), work(*)
        integer,   intent(out) :: iwork(*), info
        call dtprfs_quad(uplo, trans, diag, n, nrhs, AP, B, ldb, X, ldx, ferr, berr, work, iwork, &
                         info)
    end subroutine dtprfs

    subroutine ztprfs(uplo, trans, diag, n, nrhs, AP, B, ldb, X, ldx, ferr, berr, work, rwork, info)
        character,   intent(in)  :: uplo, trans, diag
        integer,     intent(in)  :: n, nrhs, ldb, ldx
        complex(ep), intent(in)  :: AP(*), B(ldb,*), X(ldx,*)
        real(ep),    intent(out) :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call ztprfs_quad(uplo, trans, diag, n, nrhs, AP, B, ldb, X, ldx, ferr, berr, work, rwork, &
                         info)
    end subroutine ztprfs

    subroutine dtbrfs(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, X, ldx, ferr, berr, work, &
                      iwork, info)
        character, intent(in)  :: uplo, trans, diag
        integer,   intent(in)  :: n, kd, nrhs, ldab, ldb, ldx
        real(ep),  intent(in)  :: AB(ldab,*), B(ldb,*), X(ldx,*)
        real(ep),  intent(out) :: ferr(*), berr(*), work(*)
        integer,   intent(out) :: iwork(*), info
        call dtbrfs_quad(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, X, ldx, ferr, berr, work, &
                         iwork, info)
    end subroutine dtbrfs

    subroutine ztbrfs(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, X, ldx, ferr, berr, work, &
                      rwork, info)
        character,   intent(in)  :: uplo, trans, diag
        integer,     intent(in)  :: n, kd, nrhs, ldab, ldb, ldx
        complex(ep), intent(in)  :: AB(ldab,*), B(ldb,*), X(ldx,*)
        real(ep),    intent(out) :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call ztbrfs_quad(uplo, trans, diag, n, kd, nrhs, AB, ldab, B, ldb, X, ldx, ferr, berr, work, &
                         rwork, info)
    end subroutine ztbrfs

    subroutine dgeequ(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
        integer,  intent(in)  :: m, n, lda
        real(ep), intent(in)  :: A(lda,*)
        real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
        integer,  intent(out) :: info
        call dgeequ_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
    end subroutine dgeequ

    subroutine dgbequ(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
        integer,  intent(in)  :: m, n, kl, ku, ldab
        real(ep), intent(in)  :: AB(ldab,*)
        real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
        integer,  intent(out) :: info
        call dgbequ_quad(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
    end subroutine dgbequ

    subroutine dpoequ(n, A, lda, S, scond, amax, info)
        integer,  intent(in)  :: n, lda
        real(ep), intent(in)  :: A(lda,*)
        real(ep), intent(out) :: S(*), scond, amax
        integer,  intent(out) :: info
        call dpoequ_quad(n, A, lda, S, scond, amax, info)
    end subroutine dpoequ

    subroutine dppequ(uplo, n, AP, S, scond, amax, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n
        real(ep),  intent(in)  :: AP(*)
        real(ep),  intent(out) :: S(*), scond, amax
        integer,   intent(out) :: info
        call dppequ_quad(uplo, n, AP, S, scond, amax, info)
    end subroutine dppequ

    subroutine dpbequ(uplo, n, kd, AB, ldab, S, scond, amax, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n, kd, ldab
        real(ep),  intent(in)  :: AB(ldab,*)
        real(ep),  intent(out) :: S(*), scond, amax
        integer,   intent(out) :: info
        call dpbequ_quad(uplo, n, kd, AB, ldab, S, scond, amax, info)
    end subroutine dpbequ

    subroutine zgeequ(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
        integer,     intent(in)  :: m, n, lda
        complex(ep), intent(in)  :: A(lda,*)
        real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
        integer,     intent(out) :: info
        call zgeequ_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
    end subroutine zgeequ

    subroutine zgbequ(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
        integer,     intent(in)  :: m, n, kl, ku, ldab
        complex(ep), intent(in)  :: AB(ldab,*)
        real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
        integer,     intent(out) :: info
        call zgbequ_quad(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
    end subroutine zgbequ

    subroutine zpoequ(n, A, lda, S, scond, amax, info)
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        real(ep),    intent(out) :: S(*), scond, amax
        integer,     intent(out) :: info
        call zpoequ_quad(n, A, lda, S, scond, amax, info)
    end subroutine zpoequ

    subroutine zppequ(uplo, n, AP, S, scond, amax, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n
        complex(ep), intent(in)  :: AP(*)
        real(ep),    intent(out) :: S(*), scond, amax
        integer,     intent(out) :: info
        call zppequ_quad(uplo, n, AP, S, scond, amax, info)
    end subroutine zppequ

    subroutine zpbequ(uplo, n, kd, AB, ldab, S, scond, amax, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, kd, ldab
        complex(ep), intent(in)  :: AB(ldab,*)
        real(ep),    intent(out) :: S(*), scond, amax
        integer,     intent(out) :: info
        call zpbequ_quad(uplo, n, kd, AB, ldab, S, scond, amax, info)
    end subroutine zpbequ

    subroutine dgeequb(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
        integer,  intent(in)  :: m, n, lda
        real(ep), intent(in)  :: A(lda,*)
        real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
        integer,  intent(out) :: info
        call dgeequb_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
    end subroutine dgeequb

    subroutine dgbequb(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
        integer,  intent(in)  :: m, n, kl, ku, ldab
        real(ep), intent(in)  :: AB(ldab,*)
        real(ep), intent(out) :: R(*), C(*), rowcnd, colcnd, amax
        integer,  intent(out) :: info
        call dgbequb_quad(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
    end subroutine dgbequb

    subroutine dpoequb(n, A, lda, S, scond, amax, info)
        integer,  intent(in)  :: n, lda
        real(ep), intent(in)  :: A(lda,*)
        real(ep), intent(out) :: S(*), scond, amax
        integer,  intent(out) :: info
        call dpoequb_quad(n, A, lda, S, scond, amax, info)
    end subroutine dpoequb

    subroutine dsyequb(uplo, n, A, lda, S, scond, amax, work, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n, lda
        real(ep),  intent(in)  :: A(lda,*)
        real(ep),  intent(out) :: S(*), scond, amax, work(*)
        integer,   intent(out) :: info
        call dsyequb_quad(uplo, n, A, lda, S, scond, amax, work, info)
    end subroutine dsyequb

    subroutine zgeequb(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
        integer,     intent(in)  :: m, n, lda
        complex(ep), intent(in)  :: A(lda,*)
        real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
        integer,     intent(out) :: info
        call zgeequb_quad(m, n, A, lda, R, C, rowcnd, colcnd, amax, info)
    end subroutine zgeequb

    subroutine zgbequb(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
        integer,     intent(in)  :: m, n, kl, ku, ldab
        complex(ep), intent(in)  :: AB(ldab,*)
        real(ep),    intent(out) :: R(*), C(*), rowcnd, colcnd, amax
        integer,     intent(out) :: info
        call zgbequb_quad(m, n, kl, ku, AB, ldab, R, C, rowcnd, colcnd, amax, info)
    end subroutine zgbequb

    subroutine zpoequb(n, A, lda, S, scond, amax, info)
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        real(ep),    intent(out) :: S(*), scond, amax
        integer,     intent(out) :: info
        call zpoequb_quad(n, A, lda, S, scond, amax, info)
    end subroutine zpoequb

    subroutine zsyequb(uplo, n, A, lda, S, scond, amax, work, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        real(ep),    intent(out) :: S(*), scond, amax
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zsyequb_quad(uplo, n, A, lda, S, scond, amax, work, info)
    end subroutine zsyequb

    subroutine zheequb(uplo, n, A, lda, S, scond, amax, work, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        real(ep),    intent(out) :: S(*), scond, amax
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zheequb_quad(uplo, n, A, lda, S, scond, amax, work, info)
    end subroutine zheequb

    subroutine dsytrd(uplo, n, A, lda, D, E, tau, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: D(*), E(*), tau(*), work(*)
        integer,   intent(out)   :: info
        call dsytrd_quad(uplo, n, A, lda, D, E, tau, work, lwork, info)
    end subroutine dsytrd

    subroutine zhetrd(uplo, n, A, lda, D, E, tau, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: D(*), E(*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zhetrd_quad(uplo, n, A, lda, D, E, tau, work, lwork, info)
    end subroutine zhetrd

    subroutine dsptrd(uplo, n, AP, D, E, tau, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n
        real(ep),  intent(inout) :: AP(*)
        real(ep),  intent(out)   :: D(*), E(*), tau(*)
        integer,   intent(out)   :: info
        call dsptrd_quad(uplo, n, AP, D, E, tau, info)
    end subroutine dsptrd

    subroutine zhptrd(uplo, n, AP, D, E, tau, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: AP(*)
        real(ep),    intent(out)   :: D(*), E(*)
        complex(ep), intent(out)   :: tau(*)
        integer,     intent(out)   :: info
        call zhptrd_quad(uplo, n, AP, D, E, tau, info)
    end subroutine zhptrd

    subroutine dsbtrd(vect, uplo, n, kd, AB, ldab, D, E, Q, ldq, work, info)
        character, intent(in)    :: vect, uplo
        integer,   intent(in)    :: n, kd, ldab, ldq
        real(ep),  intent(inout) :: AB(ldab,*), Q(ldq,*)
        real(ep),  intent(out)   :: D(*), E(*), work(*)
        integer,   intent(out)   :: info
        call dsbtrd_quad(vect, uplo, n, kd, AB, ldab, D, E, Q, ldq, work, info)
    end subroutine dsbtrd

    subroutine zhbtrd(vect, uplo, n, kd, AB, ldab, D, E, Q, ldq, work, info)
        character,   intent(in)    :: vect, uplo
        integer,     intent(in)    :: n, kd, ldab, ldq
        complex(ep), intent(inout) :: AB(ldab,*), Q(ldq,*)
        real(ep),    intent(out)   :: D(*), E(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zhbtrd_quad(vect, uplo, n, kd, AB, ldab, D, E, Q, ldq, work, info)
    end subroutine zhbtrd

    subroutine dorgtr(uplo, n, A, lda, tau, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(in)    :: tau(*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dorgtr_quad(uplo, n, A, lda, tau, work, lwork, info)
    end subroutine dorgtr

    subroutine zungtr(uplo, n, A, lda, tau, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: tau(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zungtr_quad(uplo, n, A, lda, tau, work, lwork, info)
    end subroutine zungtr

    subroutine dormtr(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, uplo, trans
        integer,   intent(in)    :: m, n, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormtr_quad(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormtr

    subroutine zunmtr(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, uplo, trans
        integer,     intent(in)    :: m, n, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmtr_quad(side, uplo, trans, m, n, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmtr

    subroutine dopgtr(uplo, n, AP, tau, Q, ldq, work, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n, ldq
        real(ep),  intent(in)  :: AP(*), tau(*)
        real(ep),  intent(out) :: Q(ldq,*), work(*)
        integer,   intent(out) :: info
        call dopgtr_quad(uplo, n, AP, tau, Q, ldq, work, info)
    end subroutine dopgtr

    subroutine zupgtr(uplo, n, AP, tau, Q, ldq, work, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, ldq
        complex(ep), intent(in)  :: AP(*), tau(*)
        complex(ep), intent(out) :: Q(ldq,*), work(*)
        integer,     intent(out) :: info
        call zupgtr_quad(uplo, n, AP, tau, Q, ldq, work, info)
    end subroutine zupgtr

    subroutine dopmtr(side, uplo, trans, m, n, AP, tau, C, ldc, work, info)
        character, intent(in)    :: side, uplo, trans
        integer,   intent(in)    :: m, n, ldc
        real(ep),  intent(in)    :: AP(*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dopmtr_quad(side, uplo, trans, m, n, AP, tau, C, ldc, work, info)
    end subroutine dopmtr

    subroutine zupmtr(side, uplo, trans, m, n, AP, tau, C, ldc, work, info)
        character,   intent(in)    :: side, uplo, trans
        integer,     intent(in)    :: m, n, ldc
        complex(ep), intent(in)    :: AP(*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zupmtr_quad(side, uplo, trans, m, n, AP, tau, C, ldc, work, info)
    end subroutine zupmtr

    subroutine dgebrd(m, n, A, lda, D, E, tauq, taup, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: D(*), E(*), tauq(*), taup(*), work(*)
        integer,  intent(out)   :: info
        call dgebrd_quad(m, n, A, lda, D, E, tauq, taup, work, lwork, info)
    end subroutine dgebrd

    subroutine zgebrd(m, n, A, lda, D, E, tauq, taup, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: D(*), E(*)
        complex(ep), intent(out)   :: tauq(*), taup(*), work(*)
        integer,     intent(out)   :: info
        call zgebrd_quad(m, n, A, lda, D, E, tauq, taup, work, lwork, info)
    end subroutine zgebrd

    subroutine dorgbr(vect, m, n, k, A, lda, tau, work, lwork, info)
        character, intent(in)    :: vect
        integer,   intent(in)    :: m, n, k, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(in)    :: tau(*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dorgbr_quad(vect, m, n, k, A, lda, tau, work, lwork, info)
    end subroutine dorgbr

    subroutine zungbr(vect, m, n, k, A, lda, tau, work, lwork, info)
        character,   intent(in)    :: vect
        integer,     intent(in)    :: m, n, k, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: tau(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zungbr_quad(vect, m, n, k, A, lda, tau, work, lwork, info)
    end subroutine zungbr

    subroutine dormbr(vect, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: vect, side, trans
        integer,   intent(in)    :: m, n, k, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormbr_quad(vect, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormbr

    subroutine zunmbr(vect, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: vect, side, trans
        integer,     intent(in)    :: m, n, k, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmbr_quad(vect, side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmbr

    subroutine dspgv(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: itype, n, ldz
        real(ep),  intent(inout) :: AP(*), BP(*)
        real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
        integer,   intent(out)   :: info
        call dspgv_quad(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, info)
    end subroutine dspgv

    subroutine zhpgv(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, rwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: itype, n, ldz
        complex(ep), intent(inout) :: AP(*), BP(*)
        real(ep),    intent(out)   :: W(*), rwork(*)
        complex(ep), intent(out)   :: Z(ldz,*), work(*)
        integer,     intent(out)   :: info
        call zhpgv_quad(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, rwork, info)
    end subroutine zhpgv

    subroutine dsbgv(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, work, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, ka, kb, ldab, ldbb, ldz
        real(ep),  intent(inout) :: AB(ldab,*), BB(ldbb,*)
        real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
        integer,   intent(out)   :: info
        call dsbgv_quad(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, work, info)
    end subroutine dsbgv

    subroutine zhbgv(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, work, rwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, ka, kb, ldab, ldbb, ldz
        complex(ep), intent(inout) :: AB(ldab,*), BB(ldbb,*)
        real(ep),    intent(out)   :: W(*), rwork(*)
        complex(ep), intent(out)   :: Z(ldz,*), work(*)
        integer,     intent(out)   :: info
        call zhbgv_quad(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, work, rwork, info)
    end subroutine zhbgv

    subroutine dsygvd(itype, jobz, uplo, n, A, lda, B, ldb, W, work, lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: itype, n, lda, ldb, lwork, liwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: W(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsygvd_quad(itype, jobz, uplo, n, A, lda, B, ldb, W, work, lwork, iwork, liwork, info)
    end subroutine dsygvd

    subroutine zhegvd(itype, jobz, uplo, n, A, lda, B, ldb, W, work, lwork, rwork, lrwork, iwork, &
                      liwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: itype, n, lda, ldb, lwork, lrwork, liwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),    intent(out)   :: W(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: iwork(*), info
        call zhegvd_quad(itype, jobz, uplo, n, A, lda, B, ldb, W, work, lwork, rwork, lrwork, iwork, &
                         liwork, info)
    end subroutine zhegvd

    subroutine dspgvd(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: itype, n, ldz, lwork, liwork
        real(ep),  intent(inout) :: AP(*), BP(*)
        real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dspgvd_quad(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, lwork, iwork, liwork, info)
    end subroutine dspgvd

    subroutine zhpgvd(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, lwork, rwork, lrwork, iwork, &
                      liwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: itype, n, ldz, lwork, lrwork, liwork
        complex(ep), intent(inout) :: AP(*), BP(*)
        real(ep),    intent(out)   :: W(*), rwork(*)
        complex(ep), intent(out)   :: Z(ldz,*), work(*)
        integer,     intent(out)   :: iwork(*), info
        call zhpgvd_quad(itype, jobz, uplo, n, AP, BP, W, Z, ldz, work, lwork, rwork, lrwork, iwork, &
                         liwork, info)
    end subroutine zhpgvd

    subroutine dsbgvd(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, work, lwork, iwork, &
                      liwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, ka, kb, ldab, ldbb, ldz, lwork, liwork
        real(ep),  intent(inout) :: AB(ldab,*), BB(ldbb,*)
        real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsbgvd_quad(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, work, lwork, iwork, &
                         liwork, info)
    end subroutine dsbgvd

    subroutine zhbgvd(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, work, lwork, rwork, &
                      lrwork, iwork, liwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, ka, kb, ldab, ldbb, ldz, lwork, &
        lrwork, liwork
        complex(ep), intent(inout) :: AB(ldab,*), BB(ldbb,*)
        real(ep),    intent(out)   :: W(*), rwork(*)
        complex(ep), intent(out)   :: Z(ldz,*), work(*)
        integer,     intent(out)   :: iwork(*), info
        call zhbgvd_quad(jobz, uplo, n, ka, kb, AB, ldab, BB, ldbb, W, Z, ldz, work, lwork, rwork, &
                         lrwork, iwork, liwork, info)
    end subroutine zhbgvd

    subroutine dsterf(n, D, E, info)
        integer,  intent(in)    :: n
        real(ep), intent(inout) :: D(*), E(*)
        integer,  intent(out)   :: info
        call dsterf_quad(n, D, E, info)
    end subroutine dsterf

    subroutine dsteqr(compz, n, D, E, Z, ldz, work, info)
        character, intent(in)    :: compz
        integer,   intent(in)    :: n, ldz
        real(ep),  intent(inout) :: D(*), E(*), Z(ldz,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dsteqr_quad(compz, n, D, E, Z, ldz, work, info)
    end subroutine dsteqr

    subroutine zsteqr(compz, n, D, E, Z, ldz, work, info)
        character,   intent(in)    :: compz
        integer,     intent(in)    :: n, ldz
        real(ep),    intent(inout) :: D(*), E(*)
        complex(ep), intent(inout) :: Z(ldz,*)
        real(ep),    intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zsteqr_quad(compz, n, D, E, Z, ldz, work, info)
    end subroutine zsteqr

    subroutine dstedc(compz, n, D, E, Z, ldz, work, lwork, iwork, liwork, info)
        character, intent(in)    :: compz
        integer,   intent(in)    :: n, ldz, lwork, liwork
        real(ep),  intent(inout) :: D(*), E(*), Z(ldz,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: iwork(*), info
        call dstedc_quad(compz, n, D, E, Z, ldz, work, lwork, iwork, liwork, info)
    end subroutine dstedc

    subroutine zstedc(compz, n, D, E, Z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info)
        character,   intent(in)    :: compz
        integer,     intent(in)    :: n, ldz, lwork, lrwork, liwork
        real(ep),    intent(inout) :: D(*), E(*)
        complex(ep), intent(inout) :: Z(ldz,*)
        real(ep),    intent(out)   :: rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: iwork(*), info
        call zstedc_quad(compz, n, D, E, Z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info)
    end subroutine zstedc

    subroutine dstebz(range, order, n, vl, vu, il, iu, abstol, D, E, m, nsplit, W, iblock, isplit, &
                      work, iwork, info)
        character, intent(in)  :: range, order
        integer,   intent(in)  :: n, il, iu
        real(ep),  intent(in)  :: vl, vu, abstol, D(*), E(*)
        integer,   intent(out) :: m, nsplit
        real(ep),  intent(out) :: W(*), work(*)
        integer,   intent(out) :: iblock(*), isplit(*), iwork(*), info
        call dstebz_quad(range, order, n, vl, vu, il, iu, abstol, D, E, m, nsplit, W, iblock, isplit, &
                         work, iwork, info)
    end subroutine dstebz

    subroutine dbdsqr(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc, work, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, ncvt, nru, ncc, ldvt, ldu, ldc
        real(ep),  intent(inout) :: D(*), E(*), VT(ldvt,*), U(ldu,*), C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dbdsqr_quad(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc, work, info)
    end subroutine dbdsqr

    subroutine zbdsqr(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc, rwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, ncvt, nru, ncc, ldvt, ldu, ldc
        real(ep),    intent(inout) :: D(*), E(*)
        complex(ep), intent(inout) :: VT(ldvt,*), U(ldu,*), C(ldc,*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: info
        call zbdsqr_quad(uplo, n, ncvt, nru, ncc, D, E, VT, ldvt, U, ldu, C, ldc, rwork, info)
    end subroutine zbdsqr

    subroutine dbdsdc(uplo, compq, n, D, E, U, ldu, VT, ldvt, Q, IQ, work, iwork, info)
        character, intent(in)    :: uplo, compq
        integer,   intent(in)    :: n, ldu, ldvt
        real(ep),  intent(inout) :: D(*), E(*)
        real(ep),  intent(out)   :: U(ldu,*), VT(ldvt,*), Q(*), work(*)
        integer,   intent(out)   :: IQ(*), iwork(*), info
        call dbdsdc_quad(uplo, compq, n, D, E, U, ldu, VT, ldvt, Q, IQ, work, iwork, info)
    end subroutine dbdsdc

    subroutine ddisna(job, m, n, D, sep, info)
        character, intent(in)  :: job
        integer,   intent(in)  :: m, n
        real(ep),  intent(in)  :: D(*)
        real(ep),  intent(out) :: sep(*)
        integer,   intent(out) :: info
        call ddisna_quad(job, m, n, D, sep, info)
    end subroutine ddisna

    subroutine dgebal(job, n, A, lda, ilo, ihi, scale, info)
        character, intent(in)    :: job
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: ilo, ihi
        real(ep),  intent(out)   :: scale(*)
        integer,   intent(out)   :: info
        call dgebal_quad(job, n, A, lda, ilo, ihi, scale, info)
    end subroutine dgebal

    subroutine zgebal(job, n, A, lda, ilo, ihi, scale, info)
        character,   intent(in)    :: job
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ilo, ihi
        real(ep),    intent(out)   :: scale(*)
        integer,     intent(out)   :: info
        call zgebal_quad(job, n, A, lda, ilo, ihi, scale, info)
    end subroutine zgebal

    subroutine dgebak(job, side, n, ilo, ihi, scale, m, V, ldv, info)
        character, intent(in)    :: job, side
        integer,   intent(in)    :: n, ilo, ihi, m, ldv
        real(ep),  intent(in)    :: scale(*)
        real(ep),  intent(inout) :: V(ldv,*)
        integer,   intent(out)   :: info
        call dgebak_quad(job, side, n, ilo, ihi, scale, m, V, ldv, info)
    end subroutine dgebak

    subroutine zgebak(job, side, n, ilo, ihi, scale, m, V, ldv, info)
        character,   intent(in)    :: job, side
        integer,     intent(in)    :: n, ilo, ihi, m, ldv
        real(ep),    intent(in)    :: scale(*)
        complex(ep), intent(inout) :: V(ldv,*)
        integer,     intent(out)   :: info
        call zgebak_quad(job, side, n, ilo, ihi, scale, m, V, ldv, info)
    end subroutine zgebak

    subroutine dggbal(job, n, A, lda, B, ldb, ilo, ihi, lscale, rscale, work, info)
        character, intent(in)    :: job
        integer,   intent(in)    :: n, lda, ldb
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        integer,   intent(out)   :: ilo, ihi
        real(ep),  intent(out)   :: lscale(*), rscale(*), work(*)
        integer,   intent(out)   :: info
        call dggbal_quad(job, n, A, lda, B, ldb, ilo, ihi, lscale, rscale, work, info)
    end subroutine dggbal

    subroutine zggbal(job, n, A, lda, B, ldb, ilo, ihi, lscale, rscale, work, info)
        character,   intent(in)    :: job
        integer,     intent(in)    :: n, lda, ldb
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ilo, ihi
        real(ep),    intent(out)   :: lscale(*), rscale(*), work(*)
        integer,     intent(out)   :: info
        call zggbal_quad(job, n, A, lda, B, ldb, ilo, ihi, lscale, rscale, work, info)
    end subroutine zggbal

    subroutine dggbak(job, side, n, ilo, ihi, lscale, rscale, m, V, ldv, info)
        character, intent(in)    :: job, side
        integer,   intent(in)    :: n, ilo, ihi, m, ldv
        real(ep),  intent(in)    :: lscale(*), rscale(*)
        real(ep),  intent(inout) :: V(ldv,*)
        integer,   intent(out)   :: info
        call dggbak_quad(job, side, n, ilo, ihi, lscale, rscale, m, V, ldv, info)
    end subroutine dggbak

    subroutine zggbak(job, side, n, ilo, ihi, lscale, rscale, m, V, ldv, info)
        character,   intent(in)    :: job, side
        integer,     intent(in)    :: n, ilo, ihi, m, ldv
        real(ep),    intent(in)    :: lscale(*), rscale(*)
        complex(ep), intent(inout) :: V(ldv,*)
        integer,     intent(out)   :: info
        call zggbak_quad(job, side, n, ilo, ihi, lscale, rscale, m, V, ldv, info)
    end subroutine zggbak

    subroutine dgghrd(compq, compz, n, ilo, ihi, A, lda, B, ldb, Q, ldq, Z, ldz, info)
        character, intent(in)    :: compq, compz
        integer,   intent(in)    :: n, ilo, ihi, lda, ldb, ldq, ldz
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
        integer,   intent(out)   :: info
        call dgghrd_quad(compq, compz, n, ilo, ihi, A, lda, B, ldb, Q, ldq, Z, ldz, info)
    end subroutine dgghrd

    subroutine zgghrd(compq, compz, n, ilo, ihi, A, lda, B, ldb, Q, ldq, Z, ldz, info)
        character,   intent(in)    :: compq, compz
        integer,     intent(in)    :: n, ilo, ihi, lda, ldb, ldq, ldz
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
        integer,     intent(out)   :: info
        call zgghrd_quad(compq, compz, n, ilo, ihi, A, lda, B, ldb, Q, ldq, Z, ldz, info)
    end subroutine zgghrd

    subroutine dggqrf(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
        integer,  intent(in)    :: n, m, p, lda, ldb, lwork
        real(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep), intent(out)   :: taua(*), taub(*), work(*)
        integer,  intent(out)   :: info
        call dggqrf_quad(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
    end subroutine dggqrf

    subroutine zggqrf(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
        integer,     intent(in)    :: n, m, p, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: taua(*), taub(*), work(*)
        integer,     intent(out)   :: info
        call zggqrf_quad(n, m, p, A, lda, taua, B, ldb, taub, work, lwork, info)
    end subroutine zggqrf

    subroutine dggrqf(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
        integer,  intent(in)    :: m, p, n, lda, ldb, lwork
        real(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep), intent(out)   :: taua(*), taub(*), work(*)
        integer,  intent(out)   :: info
        call dggrqf_quad(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
    end subroutine dggrqf

    subroutine zggrqf(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
        integer,     intent(in)    :: m, p, n, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: taua(*), taub(*), work(*)
        integer,     intent(out)   :: info
        call zggrqf_quad(m, p, n, A, lda, taua, B, ldb, taub, work, lwork, info)
    end subroutine zggrqf

    subroutine dtrsyl(trana, tranb, isgn, m, n, A, lda, B, ldb, C, ldc, scale, info)
        character, intent(in)    :: trana, tranb
        integer,   intent(in)    :: isgn, m, n, lda, ldb, ldc
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: scale
        integer,   intent(out)   :: info
        call dtrsyl_quad(trana, tranb, isgn, m, n, A, lda, B, ldb, C, ldc, scale, info)
    end subroutine dtrsyl

    subroutine ztrsyl(trana, tranb, isgn, m, n, A, lda, B, ldb, C, ldc, scale, info)
        character,   intent(in)    :: trana, tranb
        integer,     intent(in)    :: isgn, m, n, lda, ldb, ldc
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
        complex(ep), intent(inout) :: C(ldc,*)
        real(ep),    intent(out)   :: scale
        integer,     intent(out)   :: info
        call ztrsyl_quad(trana, tranb, isgn, m, n, A, lda, B, ldb, C, ldc, scale, info)
    end subroutine ztrsyl

    subroutine dhseqr(job, compz, n, ilo, ihi, H, ldh, WR, WI, Z, ldz, work, lwork, info)
        character, intent(in)    :: job, compz
        integer,   intent(in)    :: n, ilo, ihi, ldh, ldz, lwork
        real(ep),  intent(inout) :: H(ldh,*), Z(ldz,*)
        real(ep),  intent(out)   :: WR(*), WI(*), work(*)
        integer,   intent(out)   :: info
        call dhseqr_quad(job, compz, n, ilo, ihi, H, ldh, WR, WI, Z, ldz, work, lwork, info)
    end subroutine dhseqr

    subroutine zhseqr(job, compz, n, ilo, ihi, H, ldh, W, Z, ldz, work, lwork, info)
        character,   intent(in)    :: job, compz
        integer,     intent(in)    :: n, ilo, ihi, ldh, ldz, lwork
        complex(ep), intent(inout) :: H(ldh,*), Z(ldz,*)
        complex(ep), intent(out)   :: W(*), work(*)
        integer,     intent(out)   :: info
        call zhseqr_quad(job, compz, n, ilo, ihi, H, ldh, W, Z, ldz, work, lwork, info)
    end subroutine zhseqr

    subroutine dtrevc(side, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, mm, m, work, info)
        character, intent(in)    :: side, howmny
        integer,   intent(in)    :: n, ldt, ldvl, ldvr, mm
        logical,   intent(inout) :: sel(*)
        real(ep),  intent(in)    :: T(ldt,*)
        real(ep),  intent(inout) :: VL(ldvl,*), VR(ldvr,*)
        integer,   intent(out)   :: m
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dtrevc_quad(side, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, mm, m, work, info)
    end subroutine dtrevc

    subroutine ztrevc(side, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, mm, m, work, rwork, info)
        character,   intent(in)    :: side, howmny
        integer,     intent(in)    :: n, ldt, ldvl, ldvr, mm
        logical,     intent(in)    :: sel(*)
        complex(ep), intent(inout) :: T(ldt,*), VL(ldvl,*), VR(ldvr,*)
        integer,     intent(out)   :: m
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: info
        call ztrevc_quad(side, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, mm, m, work, rwork, info)
    end subroutine ztrevc

    subroutine dggev(jobvl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, VL, ldvl, VR, ldvr, work, &
                     lwork, info)
        character, intent(in)    :: jobvl, jobvr
        integer,   intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
        real(ep),  intent(out)   :: VL(ldvl,*), VR(ldvr,*), work(*)
        integer,   intent(out)   :: info
        call dggev_quad(jobvl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, VL, ldvl, VR, ldvr, &
                        work, lwork, info)
    end subroutine dggev

    subroutine zggev(jobvl, jobvr, n, A, lda, B, ldb, alpha, beta, VL, ldvl, VR, ldvr, work, lwork, &
                     rwork, info)
        character,   intent(in)    :: jobvl, jobvr
        integer,     intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: alpha(*), beta(*)
        complex(ep), intent(out)   :: VL(ldvl,*), VR(ldvr,*), work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: info
        call zggev_quad(jobvl, jobvr, n, A, lda, B, ldb, alpha, beta, VL, ldvl, VR, ldvr, work, &
                        lwork, rwork, info)
    end subroutine zggev

    subroutine dgeqr(m, n, A, lda, T, tsize, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, tsize, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: T(*), work(*)
        integer,  intent(out)   :: info
        call dgeqr_quad(m, n, A, lda, T, tsize, work, lwork, info)
    end subroutine dgeqr

    subroutine zgeqr(m, n, A, lda, T, tsize, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, tsize, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: T(*), work(*)
        integer,     intent(out)   :: info
        call zgeqr_quad(m, n, A, lda, T, tsize, work, lwork, info)
    end subroutine zgeqr

    subroutine dgelq(m, n, A, lda, T, tsize, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, tsize, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: T(*), work(*)
        integer,  intent(out)   :: info
        call dgelq_quad(m, n, A, lda, T, tsize, work, lwork, info)
    end subroutine dgelq

    subroutine zgelq(m, n, A, lda, T, tsize, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, tsize, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: T(*), work(*)
        integer,     intent(out)   :: info
        call zgelq_quad(m, n, A, lda, T, tsize, work, lwork, info)
    end subroutine zgelq

    subroutine dgemqr(side, trans, m, n, k, A, lda, T, tsize, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, lda, tsize, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), T(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dgemqr_quad(side, trans, m, n, k, A, lda, T, tsize, C, ldc, work, lwork, info)
    end subroutine dgemqr

    subroutine zgemqr(side, trans, m, n, k, A, lda, T, tsize, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, lda, tsize, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), T(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgemqr_quad(side, trans, m, n, k, A, lda, T, tsize, C, ldc, work, lwork, info)
    end subroutine zgemqr

    subroutine dgemlq(side, trans, m, n, k, A, lda, T, tsize, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, lda, tsize, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), T(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dgemlq_quad(side, trans, m, n, k, A, lda, T, tsize, C, ldc, work, lwork, info)
    end subroutine dgemlq

    subroutine zgemlq(side, trans, m, n, k, A, lda, T, tsize, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, lda, tsize, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), T(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgemlq_quad(side, trans, m, n, k, A, lda, T, tsize, C, ldc, work, lwork, info)
    end subroutine zgemlq

    subroutine dgeqrt(m, n, nb, A, lda, T, ldt, work, info)
        integer,  intent(in)    :: m, n, nb, lda, ldt
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: T(ldt,*), work(*)
        integer,  intent(out)   :: info
        call dgeqrt_quad(m, n, nb, A, lda, T, ldt, work, info)
    end subroutine dgeqrt

    subroutine zgeqrt(m, n, nb, A, lda, T, ldt, work, info)
        integer,     intent(in)    :: m, n, nb, lda, ldt
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: T(ldt,*), work(*)
        integer,     intent(out)   :: info
        call zgeqrt_quad(m, n, nb, A, lda, T, ldt, work, info)
    end subroutine zgeqrt

    subroutine dgelqt(m, n, mb, A, lda, T, ldt, work, info)
        integer,  intent(in)    :: m, n, mb, lda, ldt
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: T(ldt,*), work(*)
        integer,  intent(out)   :: info
        call dgelqt_quad(m, n, mb, A, lda, T, ldt, work, info)
    end subroutine dgelqt

    subroutine zgelqt(m, n, mb, A, lda, T, ldt, work, info)
        integer,     intent(in)    :: m, n, mb, lda, ldt
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: T(ldt,*), work(*)
        integer,     intent(out)   :: info
        call zgelqt_quad(m, n, mb, A, lda, T, ldt, work, info)
    end subroutine zgelqt

    subroutine dgemqrt(side, trans, m, n, k, nb, V, ldv, T, ldt, C, ldc, work, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, nb, ldv, ldt, ldc
        real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dgemqrt_quad(side, trans, m, n, k, nb, V, ldv, T, ldt, C, ldc, work, info)
    end subroutine dgemqrt

    subroutine zgemqrt(side, trans, m, n, k, nb, V, ldv, T, ldt, C, ldc, work, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, nb, ldv, ldt, ldc
        complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgemqrt_quad(side, trans, m, n, k, nb, V, ldv, T, ldt, C, ldc, work, info)
    end subroutine zgemqrt

    subroutine dgemlqt(side, trans, m, n, k, mb, V, ldv, T, ldt, C, ldc, work, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, mb, ldv, ldt, ldc
        real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dgemlqt_quad(side, trans, m, n, k, mb, V, ldv, T, ldt, C, ldc, work, info)
    end subroutine dgemlqt

    subroutine zgemlqt(side, trans, m, n, k, mb, V, ldv, T, ldt, C, ldc, work, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, mb, ldv, ldt, ldc
        complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgemlqt_quad(side, trans, m, n, k, mb, V, ldv, T, ldt, C, ldc, work, info)
    end subroutine zgemlqt

    subroutine dgetsls(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, nrhs, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dgetsls_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
    end subroutine dgetsls

    subroutine zgetsls(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgetsls_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
    end subroutine zgetsls

    subroutine dgesvx(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, R, C, B, ldb, X, ldx, &
                      rcond, ferr, berr, work, iwork, info)
        character, intent(in)    :: fact, trans
        character, intent(inout) :: equed
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
        real(ep),  intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*), R(*), C(*)
        integer,   intent(inout) :: ipiv(*)
        real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgesvx_quad(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, R, C, B, ldb, X, ldx, &
                         rcond, ferr, berr, work, iwork, info)
    end subroutine dgesvx

    subroutine zgesvx(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, R, C, B, ldb, X, ldx, &
                      rcond, ferr, berr, work, rwork, info)
        character,   intent(in)    :: fact, trans
        character,   intent(inout) :: equed
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
        complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
        real(ep),    intent(inout) :: R(*), C(*)
        integer,     intent(inout) :: ipiv(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
        integer,     intent(out)   :: info
        call zgesvx_quad(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, R, C, B, ldb, X, ldx, &
                         rcond, ferr, berr, work, rwork, info)
    end subroutine zgesvx

    subroutine dposvx(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S, B, ldb, X, ldx, rcond, ferr, &
                      berr, work, iwork, info)
        character, intent(in)    :: fact, uplo
        character, intent(inout) :: equed
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
        real(ep),  intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*), S(*)
        real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dposvx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S, B, ldb, X, ldx, rcond, &
                         ferr, berr, work, iwork, info)
    end subroutine dposvx

    subroutine zposvx(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S, B, ldb, X, ldx, rcond, ferr, &
                      berr, work, rwork, info)
        character,   intent(in)    :: fact, uplo
        character,   intent(inout) :: equed
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx
        complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
        real(ep),    intent(inout) :: S(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
        integer,     intent(out)   :: info
        call zposvx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S, B, ldb, X, ldx, rcond, &
                         ferr, berr, work, rwork, info)
    end subroutine zposvx

    subroutine dpbsvx(fact, uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, equed, S, B, ldb, X, ldx, rcond, &
                      ferr, berr, work, iwork, info)
        character, intent(in)    :: fact, uplo
        character, intent(inout) :: equed
        integer,   intent(in)    :: n, kd, nrhs, ldab, ldafb, ldb, ldx
        real(ep),  intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*), S(*)
        real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dpbsvx_quad(fact, uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, equed, S, B, ldb, X, ldx, &
                         rcond, ferr, berr, work, iwork, info)
    end subroutine dpbsvx

    subroutine zpbsvx(fact, uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, equed, S, B, ldb, X, ldx, rcond, &
                      ferr, berr, work, rwork, info)
        character,   intent(in)    :: fact, uplo
        character,   intent(inout) :: equed
        integer,     intent(in)    :: n, kd, nrhs, ldab, ldafb, ldb, ldx
        complex(ep), intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
        real(ep),    intent(inout) :: S(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
        integer,     intent(out)   :: info
        call zpbsvx_quad(fact, uplo, n, kd, nrhs, AB, ldab, AFB, ldafb, equed, S, B, ldb, X, ldx, &
                         rcond, ferr, berr, work, rwork, info)
    end subroutine zpbsvx

    subroutine dppsvx(fact, uplo, n, nrhs, AP, AFP, equed, S, B, ldb, X, ldx, rcond, ferr, berr, &
                      work, iwork, info)
        character, intent(in)    :: fact, uplo
        character, intent(inout) :: equed
        integer,   intent(in)    :: n, nrhs, ldb, ldx
        real(ep),  intent(inout) :: AP(*), AFP(*), B(ldb,*), S(*)
        real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dppsvx_quad(fact, uplo, n, nrhs, AP, AFP, equed, S, B, ldb, X, ldx, rcond, ferr, berr, &
                         work, iwork, info)
    end subroutine dppsvx

    subroutine zppsvx(fact, uplo, n, nrhs, AP, AFP, equed, S, B, ldb, X, ldx, rcond, ferr, berr, &
                      work, rwork, info)
        character,   intent(in)    :: fact, uplo
        character,   intent(inout) :: equed
        integer,     intent(in)    :: n, nrhs, ldb, ldx
        complex(ep), intent(inout) :: AP(*), AFP(*), B(ldb,*)
        real(ep),    intent(inout) :: S(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
        integer,     intent(out)   :: info
        call zppsvx_quad(fact, uplo, n, nrhs, AP, AFP, equed, S, B, ldb, X, ldx, rcond, ferr, berr, &
                         work, rwork, info)
    end subroutine zppsvx

    subroutine dgbsvx(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, equed, R, C, B, ldb, &
                      X, ldx, rcond, ferr, berr, work, iwork, info)
        character, intent(in)    :: fact, trans
        character, intent(inout) :: equed
        integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx
        real(ep),  intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*), R(*), C(*)
        integer,   intent(inout) :: ipiv(*)
        real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgbsvx_quad(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, equed, R, C, B, &
                         ldb, X, ldx, rcond, ferr, berr, work, iwork, info)
    end subroutine dgbsvx

    subroutine zgbsvx(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, equed, R, C, B, ldb, &
                      X, ldx, rcond, ferr, berr, work, rwork, info)
        character,   intent(in)    :: fact, trans
        character,   intent(inout) :: equed
        integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx
        complex(ep), intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
        real(ep),    intent(inout) :: R(*), C(*)
        integer,     intent(inout) :: ipiv(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
        integer,     intent(out)   :: info
        call zgbsvx_quad(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, equed, R, C, B, &
                         ldb, X, ldx, rcond, ferr, berr, work, rwork, info)
    end subroutine zgbsvx

    subroutine dsysvx(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, rcond, ferr, berr, &
                      work, lwork, iwork, info)
        character, intent(in)    :: fact, uplo
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, lwork
        real(ep),  intent(inout) :: AF(ldaf,*)
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        integer,   intent(inout) :: ipiv(*)
        real(ep),  intent(out)   :: X(ldx,*), rcond, ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsysvx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, rcond, ferr, &
                         berr, work, lwork, iwork, info)
    end subroutine dsysvx

    subroutine zsysvx(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, rcond, ferr, berr, &
                      work, lwork, rwork, info)
        character,   intent(in)    :: fact, uplo
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, lwork
        complex(ep), intent(inout) :: AF(ldaf,*)
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
        integer,     intent(inout) :: ipiv(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
        integer,     intent(out)   :: info
        call zsysvx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, rcond, ferr, &
                         berr, work, lwork, rwork, info)
    end subroutine zsysvx

    subroutine zhesvx(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, rcond, ferr, berr, &
                      work, lwork, rwork, info)
        character,   intent(in)    :: fact, uplo
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, lwork
        complex(ep), intent(inout) :: AF(ldaf,*)
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
        integer,     intent(inout) :: ipiv(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
        integer,     intent(out)   :: info
        call zhesvx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, B, ldb, X, ldx, rcond, ferr, &
                         berr, work, lwork, rwork, info)
    end subroutine zhesvx

    subroutine dsyev_2stage(jobz, uplo, n, A, lda, w, work, lwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: w(*), work(*)
        integer,   intent(out)   :: info
        call dsyev_2stage_quad(jobz, uplo, n, A, lda, w, work, lwork, info)
    end subroutine dsyev_2stage

    subroutine zheev_2stage(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zheev_2stage_quad(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
    end subroutine zheev_2stage

    subroutine dsyevd_2stage(jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, lda, lwork, liwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: w(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsyevd_2stage_quad(jobz, uplo, n, A, lda, w, work, lwork, iwork, liwork, info)
    end subroutine dsyevd_2stage

    subroutine zheevd_2stage(jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, &
                             info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, lda, lwork, lrwork, liwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: iwork(*), info
        call zheevd_2stage_quad(jobz, uplo, n, A, lda, w, work, lwork, rwork, lrwork, iwork, liwork, &
                                info)
    end subroutine zheevd_2stage

    subroutine dsyevr_2stage(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, &
                             isuppz, work, lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, range, uplo
        integer,   intent(in)    :: n, lda, il, iu, ldz, lwork, liwork
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        integer,   intent(out)   :: m, isuppz(*), iwork(*), info
        call dsyevr_2stage_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, &
                                isuppz, work, lwork, iwork, liwork, info)
    end subroutine dsyevr_2stage

    subroutine zheevr_2stage(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, &
                             isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
        character,   intent(in)    :: jobz, range, uplo
        integer,     intent(in)    :: n, lda, il, iu, ldz, lwork, lrwork, liwork
        real(ep),    intent(in)    :: vl, vu, abstol
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: z(ldz,*), work(*)
        integer,     intent(out)   :: m, isuppz(*), iwork(*), info
        call zheevr_2stage_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, &
                                isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
    end subroutine zheevr_2stage

    subroutine dsyevx_2stage(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, &
                             work, lwork, iwork, ifail, info)
        character, intent(in)    :: jobz, range, uplo
        integer,   intent(in)    :: n, lda, il, iu, ldz, lwork
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        integer,   intent(out)   :: m, iwork(*), ifail(*), info
        call dsyevx_2stage_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, &
                                work, lwork, iwork, ifail, info)
    end subroutine dsyevx_2stage

    subroutine zheevx_2stage(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, &
                             work, lwork, rwork, iwork, ifail, info)
        character,   intent(in)    :: jobz, range, uplo
        integer,     intent(in)    :: n, lda, il, iu, ldz, lwork
        real(ep),    intent(in)    :: vl, vu, abstol
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: z(ldz,*), work(*)
        integer,     intent(out)   :: m, iwork(*), ifail(*), info
        call zheevx_2stage_quad(jobz, range, uplo, n, A, lda, vl, vu, il, iu, abstol, m, w, z, ldz, &
                                work, lwork, rwork, iwork, ifail, info)
    end subroutine zheevx_2stage

    subroutine dsbev_2stage(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, kd, ldab, ldz, lwork
        real(ep),  intent(inout) :: AB(ldab,*)
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        integer,   intent(out)   :: info
        call dsbev_2stage_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, info)
    end subroutine dsbev_2stage

    subroutine zhbev_2stage(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, rwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, kd, ldab, ldz, lwork
        complex(ep), intent(inout) :: AB(ldab,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: z(ldz,*), work(*)
        integer,     intent(out)   :: info
        call zhbev_2stage_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, rwork, info)
    end subroutine zhbev_2stage

    subroutine dsbevd_2stage(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, iwork, liwork, &
                             info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: n, kd, ldab, ldz, lwork, liwork
        real(ep),  intent(inout) :: AB(ldab,*)
        real(ep),  intent(out)   :: w(*), z(ldz,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsbevd_2stage_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, iwork, liwork, &
                                info)
    end subroutine dsbevd_2stage

    subroutine zhbevd_2stage(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, rwork, lrwork, &
                             iwork, liwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: n, kd, ldab, ldz, lwork, lrwork, liwork
        complex(ep), intent(inout) :: AB(ldab,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: z(ldz,*), work(*)
        integer,     intent(out)   :: iwork(*), info
        call zhbevd_2stage_quad(jobz, uplo, n, kd, AB, ldab, w, z, ldz, work, lwork, rwork, lrwork, &
                                iwork, liwork, info)
    end subroutine zhbevd_2stage

    subroutine dsbevx_2stage(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, vl, vu, il, iu, abstol, m, &
                             w, z, ldz, work, lwork, iwork, ifail, info)
        character, intent(in)    :: jobz, range, uplo
        integer,   intent(in)    :: n, kd, ldab, ldq, il, iu, ldz, lwork
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: AB(ldab,*)
        real(ep),  intent(out)   :: Q(ldq,*), w(*), z(ldz,*), work(*)
        integer,   intent(out)   :: m, iwork(*), ifail(*), info
        call dsbevx_2stage_quad(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, vl, vu, il, iu, abstol, &
                                m, w, z, ldz, work, lwork, iwork, ifail, info)
    end subroutine dsbevx_2stage

    subroutine zhbevx_2stage(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, vl, vu, il, iu, abstol, m, &
                             w, z, ldz, work, lwork, rwork, iwork, ifail, info)
        character,   intent(in)    :: jobz, range, uplo
        integer,     intent(in)    :: n, kd, ldab, ldq, il, iu, ldz, lwork
        real(ep),    intent(in)    :: vl, vu, abstol
        complex(ep), intent(inout) :: AB(ldab,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: Q(ldq,*), z(ldz,*), work(*)
        integer,     intent(out)   :: m, iwork(*), ifail(*), info
        call zhbevx_2stage_quad(jobz, range, uplo, n, kd, AB, ldab, Q, ldq, vl, vu, il, iu, abstol, &
                                m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
    end subroutine zhbevx_2stage

    subroutine dsygv_2stage(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, info)
        character, intent(in)    :: jobz, uplo
        integer,   intent(in)    :: itype, n, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: w(*), work(*)
        integer,   intent(out)   :: info
        call dsygv_2stage_quad(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, info)
    end subroutine dsygv_2stage

    subroutine zhegv_2stage(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, rwork, info)
        character,   intent(in)    :: jobz, uplo
        integer,     intent(in)    :: itype, n, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zhegv_2stage_quad(itype, jobz, uplo, n, A, lda, B, ldb, w, work, lwork, rwork, info)
    end subroutine zhegv_2stage

    subroutine dgees(jobvs, sort, select, n, A, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, &
                     info)
        character, intent(in)    :: jobvs, sort
        interface
        logical function select(re, im)
        import :: ep
        real(ep), intent(in) :: re, im
        end function select
        end interface
        integer,   intent(in)    :: n, lda, ldvs, lwork
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: sdim, info
        real(ep),  intent(out)   :: wr(*), wi(*), vs(ldvs,*), work(*)
        logical,   intent(out)   :: bwork(*)
        call dgees_quad(jobvs, sort, select, n, A, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, &
                        info)
    end subroutine dgees

    subroutine zgees(jobvs, sort, select, n, A, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, &
                     info)
        character,   intent(in)    :: jobvs, sort
        interface
        logical function select(z)
        import :: ep
        complex(ep), intent(in) :: z
        end function select
        end interface
        integer,     intent(in)    :: n, lda, ldvs, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: sdim, info
        complex(ep), intent(out)   :: w(*), vs(ldvs,*), work(*)
        real(ep),    intent(out)   :: rwork(*)
        logical,     intent(out)   :: bwork(*)
        call zgees_quad(jobvs, sort, select, n, A, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, &
                        info)
    end subroutine zgees

    subroutine dgeesx(jobvs, sort, select, sense, n, A, lda, sdim, wr, wi, vs, ldvs, rconde, rcondv, &
                      work, lwork, iwork, liwork, bwork, info)
        character, intent(in)    :: jobvs, sort, sense
        interface
        logical function select(re, im)
        import :: ep
        real(ep), intent(in) :: re, im
        end function select
        end interface
        integer,   intent(in)    :: n, lda, ldvs, lwork, liwork
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: sdim, iwork(*), info
        real(ep),  intent(out)   :: wr(*), wi(*), vs(ldvs,*), work(*)
        real(ep),  intent(out)   :: rconde, rcondv
        logical,   intent(out)   :: bwork(*)
        call dgeesx_quad(jobvs, sort, select, sense, n, A, lda, sdim, wr, wi, vs, ldvs, rconde, &
                         rcondv, work, lwork, iwork, liwork, bwork, info)
    end subroutine dgeesx

    subroutine zgeesx(jobvs, sort, select, sense, n, A, lda, sdim, w, vs, ldvs, rconde, rcondv, work, &
                      lwork, rwork, bwork, info)
        character,   intent(in)    :: jobvs, sort, sense
        interface
        logical function select(z)
        import :: ep
        complex(ep), intent(in) :: z
        end function select
        end interface
        integer,     intent(in)    :: n, lda, ldvs, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: sdim, info
        complex(ep), intent(out)   :: w(*), vs(ldvs,*), work(*)
        real(ep),    intent(out)   :: rconde, rcondv, rwork(*)
        logical,     intent(out)   :: bwork(*)
        call zgeesx_quad(jobvs, sort, select, sense, n, A, lda, sdim, w, vs, ldvs, rconde, rcondv, &
                         work, lwork, rwork, bwork, info)
    end subroutine zgeesx

    subroutine dgeevx(balanc, jobvl, jobvr, sense, n, A, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, ihi, &
                      scale, abnrm, rconde, rcondv, work, lwork, iwork, info)
        character, intent(in)    :: balanc, jobvl, jobvr, sense
        integer,   intent(in)    :: n, lda, ldvl, ldvr, lwork
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: ilo, ihi, iwork(*), info
        real(ep),  intent(out)   :: wr(*), wi(*), vl(ldvl,*), vr(ldvr,*)
        real(ep),  intent(out)   :: scale(*), abnrm
        real(ep),  intent(out)   :: rconde(*), rcondv(*), work(*)
        call dgeevx_quad(balanc, jobvl, jobvr, sense, n, A, lda, wr, wi, vl, ldvl, vr, ldvr, ilo, &
                         ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info)
    end subroutine dgeevx

    subroutine zgeevx(balanc, jobvl, jobvr, sense, n, A, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, &
                      abnrm, rconde, rcondv, work, lwork, rwork, info)
        character,   intent(in)    :: balanc, jobvl, jobvr, sense
        integer,     intent(in)    :: n, lda, ldvl, ldvr, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ilo, ihi, info
        complex(ep), intent(out)   :: w(*), vl(ldvl,*), vr(ldvr,*), work(*)
        real(ep),    intent(out)   :: scale(*), abnrm
        real(ep),    intent(out)   :: rconde(*), rcondv(*), rwork(*)
        call zgeevx_quad(balanc, jobvl, jobvr, sense, n, A, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, &
                         scale, abnrm, rconde, rcondv, work, lwork, rwork, info)
    end subroutine zgeevx

    subroutine dhsein(side, eigsrc, initv, select, n, H, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m, &
                      work, ifaill, ifailr, info)
        character, intent(in)    :: side, eigsrc, initv
        integer,   intent(in)    :: n, ldh, ldvl, ldvr, mm
        logical,   intent(inout) :: select(*)
        real(ep),  intent(in)    :: H(ldh,*)
        real(ep),  intent(inout) :: wr(*)
        real(ep),  intent(in)    :: wi(*)
        real(ep),  intent(inout) :: vl(ldvl,*), vr(ldvr,*)
        integer,   intent(out)   :: m, ifaill(*), ifailr(*), info
        real(ep),  intent(out)   :: work(*)
        call dhsein_quad(side, eigsrc, initv, select, n, H, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m, &
                         work, ifaill, ifailr, info)
    end subroutine dhsein

    subroutine zhsein(side, eigsrc, initv, select, n, H, ldh, w, vl, ldvl, vr, ldvr, mm, m, work, &
                      rwork, ifaill, ifailr, info)
        character,   intent(in)    :: side, eigsrc, initv
        integer,     intent(in)    :: n, ldh, ldvl, ldvr, mm
        logical,     intent(in)    :: select(*)
        complex(ep), intent(in)    :: H(ldh,*)
        complex(ep), intent(inout) :: w(*), vl(ldvl,*), vr(ldvr,*)
        integer,     intent(out)   :: m, ifaill(*), ifailr(*), info
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rwork(*)
        call zhsein_quad(side, eigsrc, initv, select, n, H, ldh, w, vl, ldvl, vr, ldvr, mm, m, work, &
                         rwork, ifaill, ifailr, info)
    end subroutine zhsein

    subroutine dtrevc3(side, howmny, select, n, T, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, info)
        character, intent(in)    :: side, howmny
        integer,   intent(in)    :: n, ldt, ldvl, ldvr, mm, lwork
        logical,   intent(inout) :: select(*)
        real(ep),  intent(in)    :: T(ldt,*)
        real(ep),  intent(inout) :: vl(ldvl,*), vr(ldvr,*)
        integer,   intent(out)   :: m, info
        real(ep),  intent(out)   :: work(*)
        call dtrevc3_quad(side, howmny, select, n, T, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, &
                          info)
    end subroutine dtrevc3

    subroutine ztrevc3(side, howmny, select, n, T, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, &
                       rwork, lrwork, info)
        character,   intent(in)    :: side, howmny
        integer,     intent(in)    :: n, ldt, ldvl, ldvr, mm, lwork, lrwork
        logical,     intent(in)    :: select(*)
        complex(ep), intent(inout) :: T(ldt,*), vl(ldvl,*), vr(ldvr,*)
        integer,     intent(out)   :: m, info
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rwork(*)
        call ztrevc3_quad(side, howmny, select, n, T, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, &
                          rwork, lrwork, info)
    end subroutine ztrevc3

    subroutine dtrexc(compq, n, T, ldt, Q, ldq, ifst, ilst, work, info)
        character, intent(in)    :: compq
        integer,   intent(in)    :: n, ldt, ldq
        integer,   intent(inout) :: ifst, ilst
        real(ep),  intent(inout) :: T(ldt,*), Q(ldq,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dtrexc_quad(compq, n, T, ldt, Q, ldq, ifst, ilst, work, info)
    end subroutine dtrexc

    subroutine ztrexc(compq, n, T, ldt, Q, ldq, ifst, ilst, info)
        character,   intent(in)    :: compq
        integer,     intent(in)    :: n, ldt, ldq, ifst, ilst
        complex(ep), intent(inout) :: T(ldt,*), Q(ldq,*)
        integer,     intent(out)   :: info
        call ztrexc_quad(compq, n, T, ldt, Q, ldq, ifst, ilst, info)
    end subroutine ztrexc

    subroutine dtrsen(job, compq, select, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, iwork, &
                      liwork, info)
        character, intent(in)    :: job, compq
        integer,   intent(in)    :: n, ldt, ldq, lwork, liwork
        logical,   intent(in)    :: select(*)
        real(ep),  intent(inout) :: T(ldt,*), Q(ldq,*)
        real(ep),  intent(out)   :: wr(*), wi(*), s, sep, work(*)
        integer,   intent(out)   :: m, iwork(*), info
        call dtrsen_quad(job, compq, select, n, T, ldt, Q, ldq, wr, wi, m, s, sep, work, lwork, &
                         iwork, liwork, info)
    end subroutine dtrsen

    subroutine ztrsen(job, compq, select, n, T, ldt, Q, ldq, w, m, s, sep, work, lwork, info)
        character,   intent(in)    :: job, compq
        integer,     intent(in)    :: n, ldt, ldq, lwork
        logical,     intent(in)    :: select(*)
        complex(ep), intent(inout) :: T(ldt,*), Q(ldq,*)
        complex(ep), intent(out)   :: w(*), work(*)
        real(ep),    intent(out)   :: s, sep
        integer,     intent(out)   :: m, info
        call ztrsen_quad(job, compq, select, n, T, ldt, Q, ldq, w, m, s, sep, work, lwork, info)
    end subroutine ztrsen

    subroutine dgges(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, sdim, alphar, alphai, beta, &
                     vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info)
        character, intent(in) :: jobvsl, jobvsr, sort
        interface
        logical function selctg(re, im, b)
        import :: ep
        real(ep), intent(in) :: re, im, b
        end function selctg
        end interface
        integer,   intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        integer,   intent(out)   :: sdim, info
        real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
        real(ep),  intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
        logical,   intent(out)   :: bwork(*)
        call dgges_quad(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, sdim, alphar, alphai, beta, &
                        vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info)
    end subroutine dgges

    subroutine zgges(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, sdim, alpha, beta, vsl, ldvsl, &
                     vsr, ldvsr, work, lwork, rwork, bwork, info)
        character,   intent(in) :: jobvsl, jobvsr, sort
        interface
        logical function selctg(a, b)
        import :: ep
        complex(ep), intent(in) :: a, b
        end function selctg
        end interface
        integer,     intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: sdim, info
        complex(ep), intent(out)   :: alpha(*), beta(*)
        complex(ep), intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
        real(ep),    intent(out)   :: rwork(*)
        logical,     intent(out)   :: bwork(*)
        call zgges_quad(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, sdim, alpha, beta, vsl, &
                        ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info)
    end subroutine zgges

    subroutine dgges3(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, sdim, alphar, alphai, beta, &
                      vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info)
        character, intent(in) :: jobvsl, jobvsr, sort
        interface
        logical function selctg(re, im, b)
        import :: ep
        real(ep), intent(in) :: re, im, b
        end function selctg
        end interface
        integer,   intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        integer,   intent(out)   :: sdim, info
        real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
        real(ep),  intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
        logical,   intent(out)   :: bwork(*)
        call dgges3_quad(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, sdim, alphar, alphai, beta, &
                         vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info)
    end subroutine dgges3

    subroutine zgges3(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, sdim, alpha, beta, vsl, ldvsl, &
                      vsr, ldvsr, work, lwork, rwork, bwork, info)
        character,   intent(in) :: jobvsl, jobvsr, sort
        interface
        logical function selctg(a, b)
        import :: ep
        complex(ep), intent(in) :: a, b
        end function selctg
        end interface
        integer,     intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: sdim, info
        complex(ep), intent(out)   :: alpha(*), beta(*)
        complex(ep), intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
        real(ep),    intent(out)   :: rwork(*)
        logical,     intent(out)   :: bwork(*)
        call zgges3_quad(jobvsl, jobvsr, sort, selctg, n, A, lda, B, ldb, sdim, alpha, beta, vsl, &
                         ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info)
    end subroutine zgges3

    subroutine dggesx(jobvsl, jobvsr, sort, selctg, sense, n, A, lda, B, ldb, sdim, alphar, alphai, &
                      beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork, liwork, &
                      bwork, info)
        character, intent(in) :: jobvsl, jobvsr, sort, sense
        interface
        logical function selctg(re, im, b)
        import :: ep
        real(ep), intent(in) :: re, im, b
        end function selctg
        end interface
        integer,   intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork, liwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        integer,   intent(out)   :: sdim, iwork(*), info
        real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
        real(ep),  intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
        real(ep),  intent(out)   :: rconde(2), rcondv(2)
        logical,   intent(out)   :: bwork(*)
        call dggesx_quad(jobvsl, jobvsr, sort, selctg, sense, n, A, lda, B, ldb, sdim, alphar, &
                         alphai, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork, &
                         liwork, bwork, info)
    end subroutine dggesx

    subroutine zggesx(jobvsl, jobvsr, sort, selctg, sense, n, A, lda, B, ldb, sdim, alpha, beta, vsl, &
                      ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork, bwork, &
                      info)
        character,   intent(in) :: jobvsl, jobvsr, sort, sense
        interface
        logical function selctg(a, b)
        import :: ep
        complex(ep), intent(in) :: a, b
        end function selctg
        end interface
        integer,     intent(in)    :: n, lda, ldb, ldvsl, ldvsr, lwork, liwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: sdim, iwork(*), info
        complex(ep), intent(out)   :: alpha(*), beta(*)
        complex(ep), intent(out)   :: vsl(ldvsl,*), vsr(ldvsr,*), work(*)
        real(ep),    intent(out)   :: rconde(2), rcondv(2), rwork(*)
        logical,     intent(out)   :: bwork(*)
        call zggesx_quad(jobvsl, jobvsr, sort, selctg, sense, n, A, lda, B, ldb, sdim, alpha, beta, &
                         vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork, &
                         bwork, info)
    end subroutine zggesx

    subroutine dggev3(jobvl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, &
                      work, lwork, info)
        character, intent(in)    :: jobvl, jobvr
        integer,   intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        integer,   intent(out)   :: info
        real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
        real(ep),  intent(out)   :: vl(ldvl,*), vr(ldvr,*), work(*)
        call dggev3_quad(jobvl, jobvr, n, A, lda, B, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, &
                         work, lwork, info)
    end subroutine dggev3

    subroutine zggev3(jobvl, jobvr, n, A, lda, B, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, &
                      rwork, info)
        character,   intent(in)    :: jobvl, jobvr
        integer,     intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: info
        complex(ep), intent(out)   :: alpha(*), beta(*)
        complex(ep), intent(out)   :: vl(ldvl,*), vr(ldvr,*), work(*)
        real(ep),    intent(out)   :: rwork(*)
        call zggev3_quad(jobvl, jobvr, n, A, lda, B, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, &
                         lwork, rwork, info)
    end subroutine zggev3

    subroutine dggevx(balanc, jobvl, jobvr, sense, n, A, lda, B, ldb, alphar, alphai, beta, vl, ldvl, &
                      vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, &
                      iwork, bwork, info)
        character, intent(in)    :: balanc, jobvl, jobvr, sense
        integer,   intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        integer,   intent(out)   :: ilo, ihi, iwork(*), info
        real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*)
        real(ep),  intent(out)   :: vl(ldvl,*), vr(ldvr,*)
        real(ep),  intent(out)   :: lscale(*), rscale(*), abnrm, bbnrm
        real(ep),  intent(out)   :: rconde(*), rcondv(*), work(*)
        logical,   intent(out)   :: bwork(*)
        call dggevx_quad(balanc, jobvl, jobvr, sense, n, A, lda, B, ldb, alphar, alphai, beta, vl, &
                         ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, &
                         work, lwork, iwork, bwork, info)
    end subroutine dggevx

    subroutine zggevx(balanc, jobvl, jobvr, sense, n, A, lda, B, ldb, alpha, beta, vl, ldvl, vr, &
                      ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, &
                      rwork, iwork, bwork, info)
        character,   intent(in)    :: balanc, jobvl, jobvr, sense
        integer,     intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ilo, ihi, iwork(*), info
        complex(ep), intent(out)   :: alpha(*), beta(*)
        complex(ep), intent(out)   :: vl(ldvl,*), vr(ldvr,*), work(*)
        real(ep),    intent(out)   :: lscale(*), rscale(*), abnrm, bbnrm
        real(ep),    intent(out)   :: rconde(*), rcondv(*), rwork(*)
        logical,     intent(out)   :: bwork(*)
        call zggevx_quad(balanc, jobvl, jobvr, sense, n, A, lda, B, ldb, alpha, beta, vl, ldvl, vr, &
                         ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, &
                         rwork, iwork, bwork, info)
    end subroutine zggevx

    subroutine dhgeqz(job, compq, compz, n, ilo, ihi, H, ldh, T, ldt, alphar, alphai, beta, Q, ldq, &
                      Z, ldz, work, lwork, info)
        character, intent(in)    :: job, compq, compz
        integer,   intent(in)    :: n, ilo, ihi, ldh, ldt, ldq, ldz, lwork
        real(ep),  intent(inout) :: H(ldh,*), T(ldt,*), Q(ldq,*), Z(ldz,*)
        real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*), work(*)
        integer,   intent(out)   :: info
        call dhgeqz_quad(job, compq, compz, n, ilo, ihi, H, ldh, T, ldt, alphar, alphai, beta, Q, &
                         ldq, Z, ldz, work, lwork, info)
    end subroutine dhgeqz

    subroutine zhgeqz(job, compq, compz, n, ilo, ihi, H, ldh, T, ldt, alpha, beta, Q, ldq, Z, ldz, &
                      work, lwork, rwork, info)
        character,   intent(in)    :: job, compq, compz
        integer,     intent(in)    :: n, ilo, ihi, ldh, ldt, ldq, ldz, lwork
        complex(ep), intent(inout) :: H(ldh,*), T(ldt,*), Q(ldq,*), Z(ldz,*)
        complex(ep), intent(out)   :: alpha(*), beta(*), work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: info
        call zhgeqz_quad(job, compq, compz, n, ilo, ihi, H, ldh, T, ldt, alpha, beta, Q, ldq, Z, ldz, &
                         work, lwork, rwork, info)
    end subroutine zhgeqz

    subroutine dtgevc(side, howmny, select, n, S, lds, P, ldp, vl, ldvl, vr, ldvr, mm, m, work, info)
        character, intent(in)    :: side, howmny
        integer,   intent(in)    :: n, lds, ldp, ldvl, ldvr, mm
        logical,   intent(in)    :: select(*)
        real(ep),  intent(in)    :: S(lds,*), P(ldp,*)
        real(ep),  intent(inout) :: vl(ldvl,*), vr(ldvr,*)
        integer,   intent(out)   :: m, info
        real(ep),  intent(out)   :: work(*)
        call dtgevc_quad(side, howmny, select, n, S, lds, P, ldp, vl, ldvl, vr, ldvr, mm, m, work, &
                         info)
    end subroutine dtgevc

    subroutine ztgevc(side, howmny, select, n, S, lds, P, ldp, vl, ldvl, vr, ldvr, mm, m, work, &
                      rwork, info)
        character,   intent(in)    :: side, howmny
        integer,     intent(in)    :: n, lds, ldp, ldvl, ldvr, mm
        logical,     intent(in)    :: select(*)
        complex(ep), intent(in)    :: S(lds,*), P(ldp,*)
        complex(ep), intent(inout) :: vl(ldvl,*), vr(ldvr,*)
        integer,     intent(out)   :: m, info
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rwork(*)
        call ztgevc_quad(side, howmny, select, n, S, lds, P, ldp, vl, ldvl, vr, ldvr, mm, m, work, &
                         rwork, info)
    end subroutine ztgevc

    subroutine dtgexc(wantq, wantz, n, A, lda, B, ldb, Q, ldq, Z, ldz, ifst, ilst, work, lwork, info)
        logical,   intent(in)    :: wantq, wantz
        integer,   intent(in)    :: n, lda, ldb, ldq, ldz, lwork
        integer,   intent(inout) :: ifst, ilst
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dtgexc_quad(wantq, wantz, n, A, lda, B, ldb, Q, ldq, Z, ldz, ifst, ilst, work, lwork, &
                         info)
    end subroutine dtgexc

    subroutine ztgexc(wantq, wantz, n, A, lda, B, ldb, Q, ldq, Z, ldz, ifst, ilst, info)
        logical,     intent(in)    :: wantq, wantz
        integer,     intent(in)    :: n, lda, ldb, ldq, ldz, ifst, ilst
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
        integer,     intent(out)   :: info
        call ztgexc_quad(wantq, wantz, n, A, lda, B, ldb, Q, ldq, Z, ldz, ifst, ilst, info)
    end subroutine ztgexc

    subroutine dtgsen(ijob, wantq, wantz, select, n, A, lda, B, ldb, alphar, alphai, beta, Q, ldq, Z, &
                      ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info)
        integer,   intent(in)    :: ijob, n, lda, ldb, ldq, ldz, lwork, liwork
        logical,   intent(in)    :: wantq, wantz, select(*)
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
        real(ep),  intent(out)   :: alphar(*), alphai(*), beta(*), pl, pr
        real(ep),  intent(out)   :: dif(2), work(*)
        integer,   intent(out)   :: m, iwork(*), info
        call dtgsen_quad(ijob, wantq, wantz, select, n, A, lda, B, ldb, alphar, alphai, beta, Q, ldq, &
                         Z, ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info)
    end subroutine dtgsen

    subroutine ztgsen(ijob, wantq, wantz, select, n, A, lda, B, ldb, alpha, beta, Q, ldq, Z, ldz, m, &
                      pl, pr, dif, work, lwork, iwork, liwork, info)
        integer,     intent(in)    :: ijob, n, lda, ldb, ldq, ldz, lwork, liwork
        logical,     intent(in)    :: wantq, wantz, select(*)
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*), Q(ldq,*), Z(ldz,*)
        complex(ep), intent(out)   :: alpha(*), beta(*), work(*)
        real(ep),    intent(out)   :: pl, pr, dif(2)
        integer,     intent(out)   :: m, iwork(*), info
        call ztgsen_quad(ijob, wantq, wantz, select, n, A, lda, B, ldb, alpha, beta, Q, ldq, Z, ldz, &
                         m, pl, pr, dif, work, lwork, iwork, liwork, info)
    end subroutine ztgsen

    subroutine dtgsna(job, howmny, select, n, A, lda, B, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, &
                      work, lwork, iwork, info)
        character, intent(in)    :: job, howmny
        integer,   intent(in)    :: n, lda, ldb, ldvl, ldvr, mm, lwork
        logical,   intent(in)    :: select(*)
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*), vl(ldvl,*), vr(ldvr,*)
        real(ep),  intent(out)   :: s(*), dif(*), work(*)
        integer,   intent(out)   :: m, iwork(*), info
        call dtgsna_quad(job, howmny, select, n, A, lda, B, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, &
                         work, lwork, iwork, info)
    end subroutine dtgsna

    subroutine ztgsna(job, howmny, select, n, A, lda, B, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, &
                      work, lwork, iwork, info)
        character,   intent(in)    :: job, howmny
        integer,     intent(in)    :: n, lda, ldb, ldvl, ldvr, mm, lwork
        logical,     intent(in)    :: select(*)
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*), vl(ldvl,*), vr(ldvr,*)
        real(ep),    intent(out)   :: s(*), dif(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: m, iwork(*), info
        call ztgsna_quad(job, howmny, select, n, A, lda, B, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, &
                         work, lwork, iwork, info)
    end subroutine ztgsna

    subroutine dtgsja(jobu, jobv, jobq, m, p, n, k, l, A, lda, B, ldb, tola, tolb, alpha, beta, U, &
                      ldu, V, ldv, Q, ldq, work, ncycle, info)
        character, intent(in)    :: jobu, jobv, jobq
        integer,   intent(in)    :: m, p, n, k, l, lda, ldb, ldu, ldv, ldq
        real(ep),  intent(in)    :: tola, tolb
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: U(ldu,*), V(ldv,*), Q(ldq,*)
        real(ep),  intent(out)   :: alpha(*), beta(*), work(*)
        integer,   intent(out)   :: ncycle, info
        call dtgsja_quad(jobu, jobv, jobq, m, p, n, k, l, A, lda, B, ldb, tola, tolb, alpha, beta, U, &
                         ldu, V, ldv, Q, ldq, work, ncycle, info)
    end subroutine dtgsja

    subroutine ztgsja(jobu, jobv, jobq, m, p, n, k, l, A, lda, B, ldb, tola, tolb, alpha, beta, U, &
                      ldu, V, ldv, Q, ldq, work, ncycle, info)
        character,   intent(in)    :: jobu, jobv, jobq
        integer,     intent(in)    :: m, p, n, k, l, lda, ldb, ldu, ldv, ldq
        real(ep),    intent(in)    :: tola, tolb
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(inout) :: U(ldu,*), V(ldv,*), Q(ldq,*)
        real(ep),    intent(out)   :: alpha(*), beta(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: ncycle, info
        call ztgsja_quad(jobu, jobv, jobq, m, p, n, k, l, A, lda, B, ldb, tola, tolb, alpha, beta, U, &
                         ldu, V, ldv, Q, ldq, work, ncycle, info)
    end subroutine ztgsja

    subroutine dtgsyl(trans, ijob, m, n, A, lda, B, ldb, C, ldc, D, ldd, E, lde, F, ldf, scale, dif, &
                      work, lwork, iwork, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: ijob, m, n, lda, ldb, ldc, ldd, lde, ldf, lwork
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*), D(ldd,*), E(lde,*)
        real(ep),  intent(inout) :: C(ldc,*), F(ldf,*)
        real(ep),  intent(out)   :: scale, dif, work(*)
        integer,   intent(out)   :: iwork(*), info
        call dtgsyl_quad(trans, ijob, m, n, A, lda, B, ldb, C, ldc, D, ldd, E, lde, F, ldf, scale, &
                         dif, work, lwork, iwork, info)
    end subroutine dtgsyl

    subroutine ztgsyl(trans, ijob, m, n, A, lda, B, ldb, C, ldc, D, ldd, E, lde, F, ldf, scale, dif, &
                      work, lwork, iwork, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: ijob, m, n, lda, ldb, ldc, ldd, lde, ldf, lwork
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*), D(ldd,*), E(lde,*)
        complex(ep), intent(inout) :: C(ldc,*), F(ldf,*)
        real(ep),    intent(out)   :: scale, dif
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: iwork(*), info
        call ztgsyl_quad(trans, ijob, m, n, A, lda, B, ldb, C, ldc, D, ldd, E, lde, F, ldf, scale, &
                         dif, work, lwork, iwork, info)
    end subroutine ztgsyl

    subroutine dggglm(n, m, p, A, lda, B, ldb, d, x, y, work, lwork, info)
        integer,   intent(in)    :: n, m, p, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*), d(*)
        real(ep),  intent(out)   :: x(*), y(*), work(*)
        integer,   intent(out)   :: info
        call dggglm_quad(n, m, p, A, lda, B, ldb, d, x, y, work, lwork, info)
    end subroutine dggglm

    subroutine zggglm(n, m, p, A, lda, B, ldb, d, x, y, work, lwork, info)
        integer,     intent(in)    :: n, m, p, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*), d(*)
        complex(ep), intent(out)   :: x(*), y(*), work(*)
        integer,     intent(out)   :: info
        call zggglm_quad(n, m, p, A, lda, B, ldb, d, x, y, work, lwork, info)
    end subroutine zggglm

    subroutine dgglse(m, n, p, A, lda, B, ldb, c, d, x, work, lwork, info)
        integer,   intent(in)    :: m, n, p, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*), c(*), d(*)
        real(ep),  intent(out)   :: x(*), work(*)
        integer,   intent(out)   :: info
        call dgglse_quad(m, n, p, A, lda, B, ldb, c, d, x, work, lwork, info)
    end subroutine dgglse

    subroutine zgglse(m, n, p, A, lda, B, ldb, c, d, x, work, lwork, info)
        integer,     intent(in)    :: m, n, p, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*), c(*), d(*)
        complex(ep), intent(out)   :: x(*), work(*)
        integer,     intent(out)   :: info
        call zgglse_quad(m, n, p, A, lda, B, ldb, c, d, x, work, lwork, info)
    end subroutine zgglse

    subroutine dgelsd(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, iwork, info)
        integer,  intent(in)    :: m, n, nrhs, lda, ldb, lwork
        real(ep), intent(in)    :: rcond
        real(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep), intent(out)   :: S(*), work(*)
        integer,  intent(out)   :: rank, iwork(*), info
        call dgelsd_quad(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, iwork, info)
    end subroutine dgelsd

    subroutine zgelsd(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, rwork, iwork, info)
        integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
        real(ep),    intent(in)    :: rcond
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),    intent(out)   :: S(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: rank, iwork(*), info
        call zgelsd_quad(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, rwork, iwork, info)
    end subroutine zgelsd

    subroutine dgelss(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, info)
        integer,  intent(in)    :: m, n, nrhs, lda, ldb, lwork
        real(ep), intent(in)    :: rcond
        real(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep), intent(out)   :: S(*), work(*)
        integer,  intent(out)   :: rank, info
        call dgelss_quad(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, info)
    end subroutine dgelss

    subroutine zgelss(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, rwork, info)
        integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
        real(ep),    intent(in)    :: rcond
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),    intent(out)   :: S(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: rank, info
        call zgelss_quad(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, rwork, info)
    end subroutine zgelss

    subroutine dgelsy(m, n, nrhs, A, lda, B, ldb, jpvt, rcond, rank, work, lwork, info)
        integer,  intent(in)    :: m, n, nrhs, lda, ldb, lwork
        real(ep), intent(in)    :: rcond
        real(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(inout) :: jpvt(*)
        integer,  intent(out)   :: rank, info
        call dgelsy_quad(m, n, nrhs, A, lda, B, ldb, jpvt, rcond, rank, work, lwork, info)
    end subroutine dgelsy

    subroutine zgelsy(m, n, nrhs, A, lda, B, ldb, jpvt, rcond, rank, work, lwork, rwork, info)
        integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
        real(ep),    intent(in)    :: rcond
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(inout) :: jpvt(*)
        integer,     intent(out)   :: rank, info
        call zgelsy_quad(m, n, nrhs, A, lda, B, ldb, jpvt, rcond, rank, work, lwork, rwork, info)
    end subroutine zgelsy

    subroutine dgelst(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, nrhs, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dgelst_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
    end subroutine dgelst

    subroutine zgelst(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: m, n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgelst_quad(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info)
    end subroutine zgelst

    subroutine dgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, sva, U, ldu, V, ldv, work, &
                      lwork, iwork, info)
        character, intent(in)    :: joba, jobu, jobv, jobr, jobt, jobp
        integer,   intent(in)    :: m, n, lda, ldu, ldv, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: sva(*), U(ldu,*), V(ldv,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgejsv_quad(joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, sva, U, ldu, V, ldv, work, &
                         lwork, iwork, info)
    end subroutine dgejsv

    subroutine zgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, sva, U, ldu, V, ldv, cwork, &
                      lwork, rwork, lrwork, iwork, info)
        character,   intent(in)    :: joba, jobu, jobv, jobr, jobt, jobp
        integer,     intent(in)    :: m, n, lda, ldu, ldv, lwork, lrwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: sva(*), rwork(*)
        complex(ep), intent(out)   :: U(ldu,*), V(ldv,*), cwork(*)
        integer,     intent(out)   :: iwork(*), info
        call zgejsv_quad(joba, jobu, jobv, jobr, jobt, jobp, m, n, A, lda, sva, U, ldu, V, ldv, &
                         cwork, lwork, rwork, lrwork, iwork, info)
    end subroutine zgejsv

    subroutine dgesvj(joba, jobu, jobv, m, n, A, lda, sva, mv, V, ldv, work, lwork, info)
        character, intent(in)    :: joba, jobu, jobv
        integer,   intent(in)    :: m, n, lda, mv, ldv, lwork
        real(ep),  intent(inout) :: A(lda,*), V(ldv,*), work(*)
        real(ep),  intent(out)   :: sva(*)
        integer,   intent(out)   :: info
        call dgesvj_quad(joba, jobu, jobv, m, n, A, lda, sva, mv, V, ldv, work, lwork, info)
    end subroutine dgesvj

    subroutine zgesvj(joba, jobu, jobv, m, n, A, lda, sva, mv, V, ldv, cwork, lwork, rwork, lrwork, &
                      info)
        character,   intent(in)    :: joba, jobu, jobv
        integer,     intent(in)    :: m, n, lda, mv, ldv, lwork, lrwork
        complex(ep), intent(inout) :: A(lda,*), V(ldv,*), cwork(*)
        real(ep),    intent(inout) :: rwork(*)
        real(ep),    intent(out)   :: sva(*)
        integer,     intent(out)   :: info
        call zgesvj_quad(joba, jobu, jobv, m, n, A, lda, sva, mv, V, ldv, cwork, lwork, rwork, &
                         lrwork, info)
    end subroutine zgesvj

    function dlangb(norm, n, kl, ku, AB, ldab, work) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: n, kl, ku, ldab
        real(ep),  intent(in) :: AB(ldab,*)
        real(ep),  intent(out) :: work(*)
        real(ep) :: r
        r = dlangb_quad(norm, n, kl, ku, AB, ldab, work)
    end function dlangb

    function zlangb(norm, n, kl, ku, AB, ldab, work) result(r)
        character,   intent(in) :: norm
        integer,     intent(in) :: n, kl, ku, ldab
        complex(ep), intent(in) :: AB(ldab,*)
        real(ep),    intent(out) :: work(*)
        real(ep) :: r
        r = zlangb_quad(norm, n, kl, ku, AB, ldab, work)
    end function zlangb

    function dlangt(norm, n, dl, d, du) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: n
        real(ep),  intent(in) :: dl(*), d(*), du(*)
        real(ep) :: r
        r = dlangt_quad(norm, n, dl, d, du)
    end function dlangt

    function zlangt(norm, n, dl, d, du) result(r)
        character,   intent(in) :: norm
        integer,     intent(in) :: n
        complex(ep), intent(in) :: dl(*), d(*), du(*)
        real(ep) :: r
        r = zlangt_quad(norm, n, dl, d, du)
    end function zlangt

    function dlanhs(norm, n, A, lda, work) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: n, lda
        real(ep),  intent(in) :: A(lda,*)
        real(ep),  intent(out) :: work(*)
        real(ep) :: r
        r = dlanhs_quad(norm, n, A, lda, work)
    end function dlanhs

    function zlanhs(norm, n, A, lda, work) result(r)
        character,   intent(in) :: norm
        integer,     intent(in) :: n, lda
        complex(ep), intent(in) :: A(lda,*)
        real(ep),    intent(out) :: work(*)
        real(ep) :: r
        r = zlanhs_quad(norm, n, A, lda, work)
    end function zlanhs

    function dlansb(norm, uplo, n, k, AB, ldab, work) result(r)
        character, intent(in) :: norm, uplo
        integer,   intent(in) :: n, k, ldab
        real(ep),  intent(in) :: AB(ldab,*)
        real(ep),  intent(out) :: work(*)
        real(ep) :: r
        r = dlansb_quad(norm, uplo, n, k, AB, ldab, work)
    end function dlansb

    function zlanhb(norm, uplo, n, k, AB, ldab, work) result(r)
        character,   intent(in) :: norm, uplo
        integer,     intent(in) :: n, k, ldab
        complex(ep), intent(in) :: AB(ldab,*)
        real(ep),    intent(out) :: work(*)
        real(ep) :: r
        r = zlanhb_quad(norm, uplo, n, k, AB, ldab, work)
    end function zlanhb

    function dlansp(norm, uplo, n, AP, work) result(r)
        character, intent(in) :: norm, uplo
        integer,   intent(in) :: n
        real(ep),  intent(in) :: AP(*)
        real(ep),  intent(out) :: work(*)
        real(ep) :: r
        r = dlansp_quad(norm, uplo, n, AP, work)
    end function dlansp

    function zlanhp(norm, uplo, n, AP, work) result(r)
        character,   intent(in) :: norm, uplo
        integer,     intent(in) :: n
        complex(ep), intent(in) :: AP(*)
        real(ep),    intent(out) :: work(*)
        real(ep) :: r
        r = zlanhp_quad(norm, uplo, n, AP, work)
    end function zlanhp

    function dlanst(norm, n, d, e) result(r)
        character, intent(in) :: norm
        integer,   intent(in) :: n
        real(ep),  intent(in) :: d(*), e(*)
        real(ep) :: r
        r = dlanst_quad(norm, n, d, e)
    end function dlanst

    function zlanht(norm, n, d, e) result(r)
        character,   intent(in) :: norm
        integer,     intent(in) :: n
        real(ep),    intent(in) :: d(*)
        complex(ep), intent(in) :: e(*)
        real(ep) :: r
        r = zlanht_quad(norm, n, d, e)
    end function zlanht

    function dlansy(norm, uplo, n, A, lda, work) result(r)
        character, intent(in) :: norm, uplo
        integer,   intent(in) :: n, lda
        real(ep),  intent(in) :: A(lda,*)
        real(ep),  intent(out) :: work(*)
        real(ep) :: r
        r = dlansy_quad(norm, uplo, n, A, lda, work)
    end function dlansy

    function zlansy(norm, uplo, n, A, lda, work) result(r)
        character,   intent(in) :: norm, uplo
        integer,     intent(in) :: n, lda
        complex(ep), intent(in) :: A(lda,*)
        real(ep),    intent(out) :: work(*)
        real(ep) :: r
        r = zlansy_quad(norm, uplo, n, A, lda, work)
    end function zlansy

    function zlansb(norm, uplo, n, k, AB, ldab, work) result(r)
        character,   intent(in) :: norm, uplo
        integer,     intent(in) :: n, k, ldab
        complex(ep), intent(in) :: AB(ldab,*)
        real(ep),    intent(out) :: work(*)
        real(ep) :: r
        r = zlansb_quad(norm, uplo, n, k, AB, ldab, work)
    end function zlansb

    function zlansp(norm, uplo, n, AP, work) result(r)
        character,   intent(in) :: norm, uplo
        integer,     intent(in) :: n
        complex(ep), intent(in) :: AP(*)
        real(ep),    intent(out) :: work(*)
        real(ep) :: r
        r = zlansp_quad(norm, uplo, n, AP, work)
    end function zlansp

    function dlansf(norm, transr, uplo, n, A, work) result(r)
        character, intent(in) :: norm, transr, uplo
        integer,   intent(in) :: n
        real(ep),  intent(in) :: A(0:*)
        real(ep),  intent(out) :: work(0:*)
        real(ep) :: r
        r = dlansf_quad(norm, transr, uplo, n, A, work)
    end function dlansf

    function dlantb(norm, uplo, diag, n, k, AB, ldab, work) result(r)
        character, intent(in) :: norm, uplo, diag
        integer,   intent(in) :: n, k, ldab
        real(ep),  intent(in) :: AB(ldab,*)
        real(ep),  intent(out) :: work(*)
        real(ep) :: r
        r = dlantb_quad(norm, uplo, diag, n, k, AB, ldab, work)
    end function dlantb

    function zlantb(norm, uplo, diag, n, k, AB, ldab, work) result(r)
        character,   intent(in) :: norm, uplo, diag
        integer,     intent(in) :: n, k, ldab
        complex(ep), intent(in) :: AB(ldab,*)
        real(ep),    intent(out) :: work(*)
        real(ep) :: r
        r = zlantb_quad(norm, uplo, diag, n, k, AB, ldab, work)
    end function zlantb

    function dlantp(norm, uplo, diag, n, AP, work) result(r)
        character, intent(in) :: norm, uplo, diag
        integer,   intent(in) :: n
        real(ep),  intent(in) :: AP(*)
        real(ep),  intent(out) :: work(*)
        real(ep) :: r
        r = dlantp_quad(norm, uplo, diag, n, AP, work)
    end function dlantp

    function zlantp(norm, uplo, diag, n, AP, work) result(r)
        character,   intent(in) :: norm, uplo, diag
        integer,     intent(in) :: n
        complex(ep), intent(in) :: AP(*)
        real(ep),    intent(out) :: work(*)
        real(ep) :: r
        r = zlantp_quad(norm, uplo, diag, n, AP, work)
    end function zlantp

    function dlantr(norm, uplo, diag, m, n, A, lda, work) result(r)
        character, intent(in) :: norm, uplo, diag
        integer,   intent(in) :: m, n, lda
        real(ep),  intent(in) :: A(lda,*)
        real(ep),  intent(out) :: work(*)
        real(ep) :: r
        r = dlantr_quad(norm, uplo, diag, m, n, A, lda, work)
    end function dlantr

    function zlantr(norm, uplo, diag, m, n, A, lda, work) result(r)
        character,   intent(in) :: norm, uplo, diag
        integer,     intent(in) :: m, n, lda
        complex(ep), intent(in) :: A(lda,*)
        real(ep),    intent(out) :: work(*)
        real(ep) :: r
        r = zlantr_quad(norm, uplo, diag, m, n, A, lda, work)
    end function zlantr

    subroutine dlapmr(forwrd, m, n, X, ldx, K)
        logical,  intent(in)    :: forwrd
        integer,  intent(in)    :: m, n, ldx
        real(ep), intent(inout) :: X(ldx,*)
        integer,  intent(inout) :: K(*)
        call dlapmr_quad(forwrd, m, n, X, ldx, K)
    end subroutine dlapmr

    subroutine zlapmr(forwrd, m, n, X, ldx, K)
        logical,     intent(in)    :: forwrd
        integer,     intent(in)    :: m, n, ldx
        complex(ep), intent(inout) :: X(ldx,*)
        integer,     intent(inout) :: K(*)
        call zlapmr_quad(forwrd, m, n, X, ldx, K)
    end subroutine zlapmr

    subroutine dlapmt(forwrd, m, n, X, ldx, K)
        logical,  intent(in)    :: forwrd
        integer,  intent(in)    :: m, n, ldx
        real(ep), intent(inout) :: X(ldx,*)
        integer,  intent(inout) :: K(*)
        call dlapmt_quad(forwrd, m, n, X, ldx, K)
    end subroutine dlapmt

    subroutine zlapmt(forwrd, m, n, X, ldx, K)
        logical,     intent(in)    :: forwrd
        integer,     intent(in)    :: m, n, ldx
        complex(ep), intent(inout) :: X(ldx,*)
        integer,     intent(inout) :: K(*)
        call zlapmt_quad(forwrd, m, n, X, ldx, K)
    end subroutine zlapmt

    subroutine dlapll(n, X, incx, Y, incy, ssmin)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: X(*), Y(*)
        real(ep), intent(out)   :: ssmin
        call dlapll_quad(n, X, incx, Y, incy, ssmin)
    end subroutine dlapll

    subroutine zlapll(n, X, incx, Y, incy, ssmin)
        integer,     intent(in)    :: n, incx, incy
        complex(ep), intent(inout) :: X(*), Y(*)
        real(ep),    intent(out)   :: ssmin
        call zlapll_quad(n, X, incx, Y, incy, ssmin)
    end subroutine zlapll

    subroutine dlacn2(n, V, X, isgn, est, kase, isave)
        integer,  intent(in)    :: n
        real(ep), intent(out)   :: V(*)
        real(ep), intent(inout) :: X(*), est
        integer,  intent(out)   :: isgn(*)
        integer,  intent(inout) :: kase, isave(3)
        call dlacn2_quad(n, V, X, isgn, est, kase, isave)
    end subroutine dlacn2

    subroutine zlacn2(n, V, X, est, kase, isave)
        integer,     intent(in)    :: n
        complex(ep), intent(out)   :: V(*)
        complex(ep), intent(inout) :: X(*)
        real(ep),    intent(inout) :: est
        integer,     intent(inout) :: kase, isave(3)
        call zlacn2_quad(n, V, X, est, kase, isave)
    end subroutine zlacn2

    subroutine dlacon(n, V, X, isgn, est, kase)
        integer,  intent(in)    :: n
        real(ep), intent(out)   :: V(*)
        real(ep), intent(inout) :: X(*), est
        integer,  intent(out)   :: isgn(*)
        integer,  intent(inout) :: kase
        call dlacon_quad(n, V, X, isgn, est, kase)
    end subroutine dlacon

    subroutine zlacon(n, V, X, est, kase)
        integer,     intent(in)    :: n
        complex(ep), intent(out)   :: V(*)
        complex(ep), intent(inout) :: X(*)
        real(ep),    intent(inout) :: est
        integer,     intent(inout) :: kase
        call zlacon_quad(n, V, X, est, kase)
    end subroutine zlacon

    subroutine dlartg(f, g, c, s, r)
        real(ep), intent(in)  :: f, g
        real(ep), intent(out) :: c, s, r
        call dlartg_quad(f, g, c, s, r)
    end subroutine dlartg

    subroutine zlartg(f, g, c, s, r)
        complex(ep), intent(in)  :: f, g
        real(ep),    intent(out) :: c
        complex(ep), intent(out) :: s, r
        call zlartg_quad(f, g, c, s, r)
    end subroutine zlartg

    subroutine dlartgp(f, g, c, s, r)
        real(ep), intent(in)  :: f, g
        real(ep), intent(out) :: c, s, r
        call dlartgp_quad(f, g, c, s, r)
    end subroutine dlartgp

    subroutine dlartgs(x, y, sigma, c, s)
        real(ep), intent(in)  :: x, y, sigma
        real(ep), intent(out) :: c, s
        call dlartgs_quad(x, y, sigma, c, s)
    end subroutine dlartgs

    function dlapy2(x, y) result(r)
        real(ep), intent(in) :: x, y
        real(ep) :: r
        r = dlapy2_quad(x, y)
    end function dlapy2

    function dlapy3(x, y, z) result(r)
        real(ep), intent(in) :: x, y, z
        real(ep) :: r
        r = dlapy3_quad(x, y, z)
    end function dlapy3

    subroutine dladiv(a, b, c, d, p, q)
        real(ep), intent(in)  :: a, b, c, d
        real(ep), intent(out) :: p, q
        call dladiv_quad(a, b, c, d, p, q)
    end subroutine dladiv

    function zladiv(x, y) result(r)
        complex(ep), intent(in) :: x, y
        complex(ep) :: r
        r = zladiv_quad(x, y)
    end function zladiv

    subroutine dlamrg(n1, n2, A, dtrd1, dtrd2, index)
        integer,  intent(in)  :: n1, n2, dtrd1, dtrd2
        real(ep), intent(in)  :: A(*)
        integer,  intent(out) :: index(*)
        call dlamrg_quad(n1, n2, A, dtrd1, dtrd2, index)
    end subroutine dlamrg

    subroutine dlarnv(idist, iseed, n, X)
        integer,  intent(in)    :: idist, n
        integer,  intent(inout) :: iseed(4)
        real(ep), intent(out)   :: X(*)
        call dlarnv_quad(idist, iseed, n, X)
    end subroutine dlarnv

    subroutine zlarnv(idist, iseed, n, X)
        integer,     intent(in)    :: idist, n
        integer,     intent(inout) :: iseed(4)
        complex(ep), intent(out)   :: X(*)
        call zlarnv_quad(idist, iseed, n, X)
    end subroutine zlarnv

    subroutine dsbgst(vect, uplo, n, ka, kb, AB, ldab, BB, ldbb, X, ldx, work, info)
        character, intent(in)    :: vect, uplo
        integer,   intent(in)    :: n, ka, kb, ldab, ldbb, ldx
        real(ep),  intent(inout) :: AB(ldab,*)
        real(ep),  intent(in)    :: BB(ldbb,*)
        real(ep),  intent(out)   :: X(ldx,*), work(*)
        integer,   intent(out)   :: info
        call dsbgst_quad(vect, uplo, n, ka, kb, AB, ldab, BB, ldbb, X, ldx, work, info)
    end subroutine dsbgst

    subroutine zhbgst(vect, uplo, n, ka, kb, AB, ldab, BB, ldbb, X, ldx, work, rwork, info)
        character,   intent(in)    :: vect, uplo
        integer,     intent(in)    :: n, ka, kb, ldab, ldbb, ldx
        complex(ep), intent(inout) :: AB(ldab,*)
        complex(ep), intent(in)    :: BB(ldbb,*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: info
        call zhbgst_quad(vect, uplo, n, ka, kb, AB, ldab, BB, ldbb, X, ldx, work, rwork, info)
    end subroutine zhbgst

    subroutine dspgst(itype, uplo, n, AP, BP, info)
        integer,   intent(in)    :: itype, n
        character, intent(in)    :: uplo
        real(ep),  intent(inout) :: AP(*)
        real(ep),  intent(in)    :: BP(*)
        integer,   intent(out)   :: info
        call dspgst_quad(itype, uplo, n, AP, BP, info)
    end subroutine dspgst

    subroutine zhpgst(itype, uplo, n, AP, BP, info)
        integer,     intent(in)    :: itype, n
        character,   intent(in)    :: uplo
        complex(ep), intent(inout) :: AP(*)
        complex(ep), intent(in)    :: BP(*)
        integer,     intent(out)   :: info
        call zhpgst_quad(itype, uplo, n, AP, BP, info)
    end subroutine zhpgst

    subroutine dsygst(itype, uplo, n, A, lda, B, ldb, info)
        integer,   intent(in)    :: itype, n, lda, ldb
        character, intent(in)    :: uplo
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(in)    :: B(ldb,*)
        integer,   intent(out)   :: info
        call dsygst_quad(itype, uplo, n, A, lda, B, ldb, info)
    end subroutine dsygst

    subroutine zhegst(itype, uplo, n, A, lda, B, ldb, info)
        integer,     intent(in)    :: itype, n, lda, ldb
        character,   intent(in)    :: uplo
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: B(ldb,*)
        integer,     intent(out)   :: info
        call zhegst_quad(itype, uplo, n, A, lda, B, ldb, info)
    end subroutine zhegst

    subroutine dsbgvx(jobz, range, uplo, n, ka, kb, AB, ldab, BB, ldbb, Q, ldq, vl, vu, il, iu, &
                      abstol, m, w, Z, ldz, work, iwork, ifail, info)
        character, intent(in)    :: jobz, range, uplo
        integer,   intent(in)    :: n, ka, kb, ldab, ldbb, ldq, il, iu, ldz
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: AB(ldab,*), BB(ldbb,*)
        real(ep),  intent(out)   :: Q(ldq,*), w(*), Z(ldz,*), work(*)
        integer,   intent(out)   :: m, iwork(*), ifail(*), info
        call dsbgvx_quad(jobz, range, uplo, n, ka, kb, AB, ldab, BB, ldbb, Q, ldq, vl, vu, il, iu, &
                         abstol, m, w, Z, ldz, work, iwork, ifail, info)
    end subroutine dsbgvx

    subroutine zhbgvx(jobz, range, uplo, n, ka, kb, AB, ldab, BB, ldbb, Q, ldq, vl, vu, il, iu, &
                      abstol, m, w, Z, ldz, work, rwork, iwork, ifail, info)
        character,   intent(in)    :: jobz, range, uplo
        integer,     intent(in)    :: n, ka, kb, ldab, ldbb, ldq, il, iu, ldz
        real(ep),    intent(in)    :: vl, vu, abstol
        complex(ep), intent(inout) :: AB(ldab,*), BB(ldbb,*)
        complex(ep), intent(out)   :: Q(ldq,*), Z(ldz,*), work(*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        integer,     intent(out)   :: m, iwork(*), ifail(*), info
        call zhbgvx_quad(jobz, range, uplo, n, ka, kb, AB, ldab, BB, ldbb, Q, ldq, vl, vu, il, iu, &
                         abstol, m, w, Z, ldz, work, rwork, iwork, ifail, info)
    end subroutine zhbgvx

    subroutine dspgvx(itype, jobz, range, uplo, n, AP, BP, vl, vu, il, iu, abstol, m, w, Z, ldz, &
                      work, iwork, ifail, info)
        integer,   intent(in)    :: itype, n, il, iu, ldz
        character, intent(in)    :: jobz, range, uplo
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: AP(*), BP(*)
        real(ep),  intent(out)   :: w(*), Z(ldz,*), work(*)
        integer,   intent(out)   :: m, iwork(*), ifail(*), info
        call dspgvx_quad(itype, jobz, range, uplo, n, AP, BP, vl, vu, il, iu, abstol, m, w, Z, ldz, &
                         work, iwork, ifail, info)
    end subroutine dspgvx

    subroutine zhpgvx(itype, jobz, range, uplo, n, AP, BP, vl, vu, il, iu, abstol, m, w, Z, ldz, &
                      work, rwork, iwork, ifail, info)
        integer,     intent(in)    :: itype, n, il, iu, ldz
        character,   intent(in)    :: jobz, range, uplo
        real(ep),    intent(in)    :: vl, vu, abstol
        complex(ep), intent(inout) :: AP(*), BP(*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: Z(ldz,*), work(*)
        integer,     intent(out)   :: m, iwork(*), ifail(*), info
        call zhpgvx_quad(itype, jobz, range, uplo, n, AP, BP, vl, vu, il, iu, abstol, m, w, Z, ldz, &
                         work, rwork, iwork, ifail, info)
    end subroutine zhpgvx

    subroutine dsygvx(itype, jobz, range, uplo, n, A, lda, B, ldb, vl, vu, il, iu, abstol, m, w, Z, &
                      ldz, work, lwork, iwork, ifail, info)
        integer,   intent(in)    :: itype, n, lda, ldb, il, iu, ldz, lwork
        character, intent(in)    :: jobz, range, uplo
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: w(*), Z(ldz,*), work(*)
        integer,   intent(out)   :: m, iwork(*), ifail(*), info
        call dsygvx_quad(itype, jobz, range, uplo, n, A, lda, B, ldb, vl, vu, il, iu, abstol, m, w, &
                         Z, ldz, work, lwork, iwork, ifail, info)
    end subroutine dsygvx

    subroutine zhegvx(itype, jobz, range, uplo, n, A, lda, B, ldb, vl, vu, il, iu, abstol, m, w, Z, &
                      ldz, work, lwork, rwork, iwork, ifail, info)
        integer,     intent(in)    :: itype, n, lda, ldb, il, iu, ldz, lwork
        character,   intent(in)    :: jobz, range, uplo
        real(ep),    intent(in)    :: vl, vu, abstol
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),    intent(out)   :: w(*), rwork(*)
        complex(ep), intent(out)   :: Z(ldz,*), work(*)
        integer,     intent(out)   :: m, iwork(*), ifail(*), info
        call zhegvx_quad(itype, jobz, range, uplo, n, A, lda, B, ldb, vl, vu, il, iu, abstol, m, w, &
                         Z, ldz, work, lwork, rwork, iwork, ifail, info)
    end subroutine zhegvx

    subroutine dtfttp(transr, uplo, n, ARF, AP, info)
        character, intent(in)  :: transr, uplo
        integer,   intent(in)  :: n
        real(ep),  intent(in)  :: ARF(0:*)
        real(ep),  intent(out) :: AP(0:*)
        integer,   intent(out) :: info
        call dtfttp_quad(transr, uplo, n, ARF, AP, info)
    end subroutine dtfttp

    subroutine ztfttp(transr, uplo, n, ARF, AP, info)
        character,   intent(in)  :: transr, uplo
        integer,     intent(in)  :: n
        complex(ep), intent(in)  :: ARF(0:*)
        complex(ep), intent(out) :: AP(0:*)
        integer,     intent(out) :: info
        call ztfttp_quad(transr, uplo, n, ARF, AP, info)
    end subroutine ztfttp

    subroutine dtfttr(transr, uplo, n, ARF, A, lda, info)
        character, intent(in)    :: transr, uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(in)    :: ARF(0:*)
        real(ep),  intent(inout) :: A(0:lda-1, 0:*)
        integer,   intent(out)   :: info
        call dtfttr_quad(transr, uplo, n, ARF, A, lda, info)
    end subroutine dtfttr

    subroutine ztfttr(transr, uplo, n, ARF, A, lda, info)
        character,   intent(in)    :: transr, uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(in)    :: ARF(0:*)
        complex(ep), intent(inout) :: A(0:lda-1, 0:*)
        integer,     intent(out)   :: info
        call ztfttr_quad(transr, uplo, n, ARF, A, lda, info)
    end subroutine ztfttr

    subroutine dtpttf(transr, uplo, n, AP, ARF, info)
        character, intent(in)  :: transr, uplo
        integer,   intent(in)  :: n
        real(ep),  intent(in)  :: AP(0:*)
        real(ep),  intent(out) :: ARF(0:*)
        integer,   intent(out) :: info
        call dtpttf_quad(transr, uplo, n, AP, ARF, info)
    end subroutine dtpttf

    subroutine ztpttf(transr, uplo, n, AP, ARF, info)
        character,   intent(in)  :: transr, uplo
        integer,     intent(in)  :: n
        complex(ep), intent(in)  :: AP(0:*)
        complex(ep), intent(out) :: ARF(0:*)
        integer,     intent(out) :: info
        call ztpttf_quad(transr, uplo, n, AP, ARF, info)
    end subroutine ztpttf

    subroutine dtpttr(uplo, n, AP, A, lda, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(in)    :: AP(*)
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: info
        call dtpttr_quad(uplo, n, AP, A, lda, info)
    end subroutine dtpttr

    subroutine ztpttr(uplo, n, AP, A, lda, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(in)    :: AP(*)
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: info
        call ztpttr_quad(uplo, n, AP, A, lda, info)
    end subroutine ztpttr

    subroutine dtrttf(transr, uplo, n, A, lda, ARF, info)
        character, intent(in)  :: transr, uplo
        integer,   intent(in)  :: n, lda
        real(ep),  intent(in)  :: A(0:lda-1, 0:*)
        real(ep),  intent(out) :: ARF(0:*)
        integer,   intent(out) :: info
        call dtrttf_quad(transr, uplo, n, A, lda, ARF, info)
    end subroutine dtrttf

    subroutine ztrttf(transr, uplo, n, A, lda, ARF, info)
        character,   intent(in)  :: transr, uplo
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(0:lda-1, 0:*)
        complex(ep), intent(out) :: ARF(0:*)
        integer,     intent(out) :: info
        call ztrttf_quad(transr, uplo, n, A, lda, ARF, info)
    end subroutine ztrttf

    subroutine dtrttp(uplo, n, A, lda, AP, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n, lda
        real(ep),  intent(in)  :: A(lda,*)
        real(ep),  intent(out) :: AP(*)
        integer,   intent(out) :: info
        call dtrttp_quad(uplo, n, A, lda, AP, info)
    end subroutine dtrttp

    subroutine ztrttp(uplo, n, A, lda, AP, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        complex(ep), intent(out) :: AP(*)
        integer,     intent(out) :: info
        call ztrttp_quad(uplo, n, A, lda, AP, info)
    end subroutine ztrttp

    logical function disnan(din)
        real(ep), intent(in) :: din
        disnan = disnan_quad(din)
    end function disnan

    subroutine drscl(n, sa, sx, incx)
        integer,  intent(in)    :: n, incx
        real(ep), intent(in)    :: sa
        real(ep), intent(inout) :: sx(*)
        call drscl_quad(n, sa, sx, incx)
    end subroutine drscl

    subroutine zdrscl(n, sa, sx, incx)
        integer,     intent(in)    :: n, incx
        real(ep),    intent(in)    :: sa
        complex(ep), intent(inout) :: sx(*)
        call zdrscl_quad(n, sa, sx, incx)
    end subroutine zdrscl

    subroutine zrscl(n, a, x, incx)
        integer,     intent(in)    :: n, incx
        complex(ep), intent(in)    :: a
        complex(ep), intent(inout) :: x(*)
        call zrscl_quad(n, a, x, incx)
    end subroutine zrscl

    subroutine zrot(n, cx, incx, cy, incy, c, s)
        integer,     intent(in)    :: n, incx, incy
        real(ep),    intent(in)    :: c
        complex(ep), intent(in)    :: s
        complex(ep), intent(inout) :: cx(*), cy(*)
        call zrot_quad(n, cx, incx, cy, incy, c, s)
    end subroutine zrot

    subroutine dorg2l(m, n, k, A, lda, tau, work, info)
        integer,  intent(in)    :: m, n, k, lda
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: tau(*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dorg2l_quad(m, n, k, A, lda, tau, work, info)
    end subroutine dorg2l

    subroutine zung2l(m, n, k, A, lda, tau, work, info)
        integer,     intent(in)    :: m, n, k, lda
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: tau(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zung2l_quad(m, n, k, A, lda, tau, work, info)
    end subroutine zung2l

    subroutine dorg2r(m, n, k, A, lda, tau, work, info)
        integer,  intent(in)    :: m, n, k, lda
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: tau(*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dorg2r_quad(m, n, k, A, lda, tau, work, info)
    end subroutine dorg2r

    subroutine zung2r(m, n, k, A, lda, tau, work, info)
        integer,     intent(in)    :: m, n, k, lda
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: tau(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zung2r_quad(m, n, k, A, lda, tau, work, info)
    end subroutine zung2r

    subroutine dorm2l(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, lda, ldc
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dorm2l_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
    end subroutine dorm2l

    subroutine zunm2l(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, lda, ldc
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunm2l_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
    end subroutine zunm2l

    subroutine dorm2r(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, lda, ldc
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dorm2r_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
    end subroutine dorm2r

    subroutine zunm2r(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, lda, ldc
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunm2r_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, info)
    end subroutine zunm2r

    subroutine dormlq(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormlq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormlq

    subroutine zunmlq(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmlq_quad(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmlq

    subroutine dormhr(side, trans, m, n, ilo, ihi, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, ilo, ihi, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormhr_quad(side, trans, m, n, ilo, ihi, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormhr

    subroutine zunmhr(side, trans, m, n, ilo, ihi, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, ilo, ihi, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmhr_quad(side, trans, m, n, ilo, ihi, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmhr

    subroutine dtzrzf(m, n, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dtzrzf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine dtzrzf

    subroutine ztzrzf(m, n, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call ztzrzf_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine ztzrzf

    subroutine dormrz(side, trans, m, n, k, l, A, lda, tau, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, l, lda, ldc, lwork
        real(ep),  intent(in)    :: A(lda,*), tau(*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dormrz_quad(side, trans, m, n, k, l, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine dormrz

    subroutine zunmrz(side, trans, m, n, k, l, A, lda, tau, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, l, lda, ldc, lwork
        complex(ep), intent(in)    :: A(lda,*), tau(*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunmrz_quad(side, trans, m, n, k, l, A, lda, tau, C, ldc, work, lwork, info)
    end subroutine zunmrz

    subroutine dpftrf(transr, uplo, n, A, info)
        character, intent(in)    :: transr, uplo
        integer,   intent(in)    :: n
        real(ep),  intent(inout) :: A(0:*)
        integer,   intent(out)   :: info
        call dpftrf_quad(transr, uplo, n, A, info)
    end subroutine dpftrf

    subroutine zpftrf(transr, uplo, n, A, info)
        character,   intent(in)    :: transr, uplo
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: A(0:*)
        integer,     intent(out)   :: info
        call zpftrf_quad(transr, uplo, n, A, info)
    end subroutine zpftrf

    subroutine dpftri(transr, uplo, n, A, info)
        character, intent(in)    :: transr, uplo
        integer,   intent(in)    :: n
        real(ep),  intent(inout) :: A(0:*)
        integer,   intent(out)   :: info
        call dpftri_quad(transr, uplo, n, A, info)
    end subroutine dpftri

    subroutine zpftri(transr, uplo, n, A, info)
        character,   intent(in)    :: transr, uplo
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: A(0:*)
        integer,     intent(out)   :: info
        call zpftri_quad(transr, uplo, n, A, info)
    end subroutine zpftri

    subroutine dpftrs(transr, uplo, n, nrhs, A, B, ldb, info)
        character, intent(in)    :: transr, uplo
        integer,   intent(in)    :: n, nrhs, ldb
        real(ep),  intent(in)    :: A(0:*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dpftrs_quad(transr, uplo, n, nrhs, A, B, ldb, info)
    end subroutine dpftrs

    subroutine zpftrs(transr, uplo, n, nrhs, A, B, ldb, info)
        character,   intent(in)    :: transr, uplo
        integer,     intent(in)    :: n, nrhs, ldb
        complex(ep), intent(in)    :: A(0:*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zpftrs_quad(transr, uplo, n, nrhs, A, B, ldb, info)
    end subroutine zpftrs

    subroutine dpstrf(uplo, n, A, lda, piv, rank, tol, work, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: piv(*), rank, info
        real(ep),  intent(in)    :: tol
        real(ep),  intent(out)   :: work(*)
        call dpstrf_quad(uplo, n, A, lda, piv, rank, tol, work, info)
    end subroutine dpstrf

    subroutine zpstrf(uplo, n, A, lda, piv, rank, tol, work, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: piv(*), rank, info
        real(ep),    intent(in)    :: tol
        real(ep),    intent(out)   :: work(*)
        call zpstrf_quad(uplo, n, A, lda, piv, rank, tol, work, info)
    end subroutine zpstrf

    subroutine dpteqr(compz, n, D, E, Z, ldz, work, info)
        character, intent(in)    :: compz
        integer,   intent(in)    :: n, ldz
        real(ep),  intent(inout) :: D(*), E(*), Z(ldz,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dpteqr_quad(compz, n, D, E, Z, ldz, work, info)
    end subroutine dpteqr

    subroutine zpteqr(compz, n, D, E, Z, ldz, work, info)
        character,   intent(in)    :: compz
        integer,     intent(in)    :: n, ldz
        real(ep),    intent(inout) :: D(*), E(*)
        complex(ep), intent(inout) :: Z(ldz,*)
        real(ep),    intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zpteqr_quad(compz, n, D, E, Z, ldz, work, info)
    end subroutine zpteqr

    subroutine dpbstf(uplo, n, kd, AB, ldab, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, kd, ldab
        real(ep),  intent(inout) :: AB(ldab,*)
        integer,   intent(out)   :: info
        call dpbstf_quad(uplo, n, kd, AB, ldab, info)
    end subroutine dpbstf

    subroutine zpbstf(uplo, n, kd, AB, ldab, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, kd, ldab
        complex(ep), intent(inout) :: AB(ldab,*)
        integer,     intent(out)   :: info
        call zpbstf_quad(uplo, n, kd, AB, ldab, info)
    end subroutine zpbstf

    subroutine dsfrk(transr, uplo, trans, n, k, alpha, A, lda, beta, C)
        character, intent(in)    :: transr, uplo, trans
        integer,   intent(in)    :: n, k, lda
        real(ep),  intent(in)    :: alpha, beta, A(lda,*)
        real(ep),  intent(inout) :: C(0:*)
        call dsfrk_quad(transr, uplo, trans, n, k, alpha, A, lda, beta, C)
    end subroutine dsfrk

    subroutine zhfrk(transr, uplo, trans, n, k, alpha, A, lda, beta, C)
        character,   intent(in)    :: transr, uplo, trans
        integer,     intent(in)    :: n, k, lda
        real(ep),    intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: C(0:*)
        call zhfrk_quad(transr, uplo, trans, n, k, alpha, A, lda, beta, C)
    end subroutine zhfrk

    subroutine dtfsm(transr, side, uplo, trans, diag, m, n, alpha, A, B, ldb)
        character, intent(in)    :: transr, side, uplo, trans, diag
        integer,   intent(in)    :: m, n, ldb
        real(ep),  intent(in)    :: alpha, A(0:*)
        real(ep),  intent(inout) :: B(0:ldb-1, 0:*)
        call dtfsm_quad(transr, side, uplo, trans, diag, m, n, alpha, A, B, ldb)
    end subroutine dtfsm

    subroutine ztfsm(transr, side, uplo, trans, diag, m, n, alpha, A, B, ldb)
        character,   intent(in)    :: transr, side, uplo, trans, diag
        integer,     intent(in)    :: m, n, ldb
        complex(ep), intent(in)    :: alpha, A(0:*)
        complex(ep), intent(inout) :: B(0:ldb-1, 0:*)
        call ztfsm_quad(transr, side, uplo, trans, diag, m, n, alpha, A, B, ldb)
    end subroutine ztfsm

    subroutine dtftri(transr, uplo, diag, n, A, info)
        character, intent(in)    :: transr, uplo, diag
        integer,   intent(in)    :: n
        real(ep),  intent(inout) :: A(0:*)
        integer,   intent(out)   :: info
        call dtftri_quad(transr, uplo, diag, n, A, info)
    end subroutine dtftri

    subroutine ztftri(transr, uplo, diag, n, A, info)
        character,   intent(in)    :: transr, uplo, diag
        integer,     intent(in)    :: n
        complex(ep), intent(inout) :: A(0:*)
        integer,     intent(out)   :: info
        call ztftri_quad(transr, uplo, diag, n, A, info)
    end subroutine ztftri

    subroutine dtpqrt(m, n, l, nb, A, lda, B, ldb, T, ldt, work, info)
        integer,  intent(in)    :: m, n, l, nb, lda, ldb, ldt
        real(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dtpqrt_quad(m, n, l, nb, A, lda, B, ldb, T, ldt, work, info)
    end subroutine dtpqrt

    subroutine ztpqrt(m, n, l, nb, A, lda, B, ldb, T, ldt, work, info)
        integer,     intent(in)    :: m, n, l, nb, lda, ldb, ldt
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call ztpqrt_quad(m, n, l, nb, A, lda, B, ldb, T, ldt, work, info)
    end subroutine ztpqrt

    subroutine dtpqrt2(m, n, l, A, lda, B, ldb, T, ldt, info)
        integer,  intent(in)    :: m, n, l, lda, ldb, ldt
        real(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
        integer,  intent(out)   :: info
        call dtpqrt2_quad(m, n, l, A, lda, B, ldb, T, ldt, info)
    end subroutine dtpqrt2

    subroutine ztpqrt2(m, n, l, A, lda, B, ldb, T, ldt, info)
        integer,     intent(in)    :: m, n, l, lda, ldb, ldt
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
        integer,     intent(out)   :: info
        call ztpqrt2_quad(m, n, l, A, lda, B, ldb, T, ldt, info)
    end subroutine ztpqrt2

    subroutine dtpmqrt(side, trans, m, n, k, l, nb, V, ldv, T, ldt, A, lda, B, ldb, work, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, l, nb, ldv, ldt, lda, ldb
        real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dtpmqrt_quad(side, trans, m, n, k, l, nb, V, ldv, T, ldt, A, lda, B, ldb, work, info)
    end subroutine dtpmqrt

    subroutine ztpmqrt(side, trans, m, n, k, l, nb, V, ldv, T, ldt, A, lda, B, ldb, work, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, l, nb, ldv, ldt, lda, ldb
        complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call ztpmqrt_quad(side, trans, m, n, k, l, nb, V, ldv, T, ldt, A, lda, B, ldb, work, info)
    end subroutine ztpmqrt

    subroutine dtplqt(m, n, l, mb, A, lda, B, ldb, T, ldt, work, info)
        integer,  intent(in)    :: m, n, l, mb, lda, ldb, ldt
        real(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dtplqt_quad(m, n, l, mb, A, lda, B, ldb, T, ldt, work, info)
    end subroutine dtplqt

    subroutine ztplqt(m, n, l, mb, A, lda, B, ldb, T, ldt, work, info)
        integer,     intent(in)    :: m, n, l, mb, lda, ldb, ldt
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call ztplqt_quad(m, n, l, mb, A, lda, B, ldb, T, ldt, work, info)
    end subroutine ztplqt

    subroutine dtplqt2(m, n, l, A, lda, B, ldb, T, ldt, info)
        integer,  intent(in)    :: m, n, l, lda, ldb, ldt
        real(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
        integer,  intent(out)   :: info
        call dtplqt2_quad(m, n, l, A, lda, B, ldb, T, ldt, info)
    end subroutine dtplqt2

    subroutine ztplqt2(m, n, l, A, lda, B, ldb, T, ldt, info)
        integer,     intent(in)    :: m, n, l, lda, ldb, ldt
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*), T(ldt,*)
        integer,     intent(out)   :: info
        call ztplqt2_quad(m, n, l, A, lda, B, ldb, T, ldt, info)
    end subroutine ztplqt2

    subroutine dtpmlqt(side, trans, m, n, k, l, mb, V, ldv, T, ldt, A, lda, B, ldb, work, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, k, l, mb, ldv, ldt, lda, ldb
        real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dtpmlqt_quad(side, trans, m, n, k, l, mb, V, ldv, T, ldt, A, lda, B, ldb, work, info)
    end subroutine dtpmlqt

    subroutine ztpmlqt(side, trans, m, n, k, l, mb, V, ldv, T, ldt, A, lda, B, ldb, work, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, k, l, mb, ldv, ldt, lda, ldb
        complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call ztpmlqt_quad(side, trans, m, n, k, l, mb, V, ldv, T, ldt, A, lda, B, ldb, work, info)
    end subroutine ztpmlqt

    subroutine dtprfb(side, trans, direct, storev, m, n, k, l, V, ldv, T, ldt, A, lda, B, ldb, work, &
                      ldwork)
        character, intent(in)    :: side, trans, direct, storev
        integer,   intent(in)    :: m, n, k, l, ldv, ldt, lda, ldb, ldwork
        real(ep),  intent(in)    :: V(ldv,*), T(ldt,*)
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: work(ldwork,*)
        call dtprfb_quad(side, trans, direct, storev, m, n, k, l, V, ldv, T, ldt, A, lda, B, ldb, &
                         work, ldwork)
    end subroutine dtprfb

    subroutine ztprfb(side, trans, direct, storev, m, n, k, l, V, ldv, T, ldt, A, lda, B, ldb, work, &
                      ldwork)
        character,   intent(in)    :: side, trans, direct, storev
        integer,     intent(in)    :: m, n, k, l, ldv, ldt, lda, ldb, ldwork
        complex(ep), intent(in)    :: V(ldv,*), T(ldt,*)
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: work(ldwork,*)
        call ztprfb_quad(side, trans, direct, storev, m, n, k, l, V, ldv, T, ldt, A, lda, B, ldb, &
                         work, ldwork)
    end subroutine ztprfb

    subroutine dgeqr2p(m, n, A, lda, tau, work, info)
        integer,  intent(in)    :: m, n, lda
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dgeqr2p_quad(m, n, A, lda, tau, work, info)
    end subroutine dgeqr2p

    subroutine zgeqr2p(m, n, A, lda, tau, work, info)
        integer,     intent(in)    :: m, n, lda
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zgeqr2p_quad(m, n, A, lda, tau, work, info)
    end subroutine zgeqr2p

    subroutine dgeqrfp(m, n, A, lda, tau, work, lwork, info)
        integer,  intent(in)    :: m, n, lda, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: tau(*), work(*)
        integer,  intent(out)   :: info
        call dgeqrfp_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine dgeqrfp

    subroutine zgeqrfp(m, n, A, lda, tau, work, lwork, info)
        integer,     intent(in)    :: m, n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: tau(*), work(*)
        integer,     intent(out)   :: info
        call zgeqrfp_quad(m, n, A, lda, tau, work, lwork, info)
    end subroutine zgeqrfp

    subroutine dlatsqr(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
        integer,  intent(in)    :: m, n, mb, nb, lda, ldt, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: T(ldt,*), work(*)
        integer,  intent(out)   :: info
        call dlatsqr_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
    end subroutine dlatsqr

    subroutine zlatsqr(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
        integer,     intent(in)    :: m, n, mb, nb, lda, ldt, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: T(ldt,*), work(*)
        integer,     intent(out)   :: info
        call zlatsqr_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
    end subroutine zlatsqr

    subroutine dorgtsqr(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
        integer,  intent(in)    :: m, n, mb, nb, lda, ldt, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: T(ldt,*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dorgtsqr_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
    end subroutine dorgtsqr

    subroutine zungtsqr(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
        integer,     intent(in)    :: m, n, mb, nb, lda, ldt, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: T(ldt,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zungtsqr_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
    end subroutine zungtsqr

    subroutine dorgtsqr_row(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
        integer,  intent(in)    :: m, n, mb, nb, lda, ldt, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(in)    :: T(ldt,*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dorgtsqr_row_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
    end subroutine dorgtsqr_row

    subroutine zungtsqr_row(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
        integer,     intent(in)    :: m, n, mb, nb, lda, ldt, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: T(ldt,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zungtsqr_row_quad(m, n, mb, nb, A, lda, T, ldt, work, lwork, info)
    end subroutine zungtsqr_row

    subroutine dorhr_col(m, n, nb, A, lda, T, ldt, D, info)
        integer,  intent(in)    :: m, n, nb, lda, ldt
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: T(ldt,*), D(*)
        integer,  intent(out)   :: info
        call dorhr_col_quad(m, n, nb, A, lda, T, ldt, D, info)
    end subroutine dorhr_col

    subroutine zunhr_col(m, n, nb, A, lda, T, ldt, D, info)
        integer,     intent(in)    :: m, n, nb, lda, ldt
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: T(ldt,*), D(*)
        integer,     intent(out)   :: info
        call zunhr_col_quad(m, n, nb, A, lda, T, ldt, D, info)
    end subroutine zunhr_col

    subroutine dgetsqrhrt(m, n, mb1, nb1, nb2, A, lda, T, ldt, work, lwork, info)
        integer,  intent(in)    :: m, n, mb1, nb1, nb2, lda, ldt, lwork
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: T(ldt,*), work(*)
        integer,  intent(out)   :: info
        call dgetsqrhrt_quad(m, n, mb1, nb1, nb2, A, lda, T, ldt, work, lwork, info)
    end subroutine dgetsqrhrt

    subroutine zgetsqrhrt(m, n, mb1, nb1, nb2, A, lda, T, ldt, work, lwork, info)
        integer,     intent(in)    :: m, n, mb1, nb1, nb2, lda, ldt, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: T(ldt,*), work(*)
        integer,     intent(out)   :: info
        call zgetsqrhrt_quad(m, n, mb1, nb1, nb2, A, lda, T, ldt, work, lwork, info)
    end subroutine zgetsqrhrt

    subroutine dgeqp3rk(m, n, nrhs, kmax, abstol, reltol, A, lda, K, maxc2nrmk, relmaxc2nrmk, jpiv, &
                        tau, work, lwork, iwork, info)
        integer,  intent(in)    :: m, n, nrhs, kmax, lda, lwork
        real(ep), intent(in)    :: abstol, reltol
        real(ep), intent(inout) :: A(lda,*)
        integer,  intent(out)   :: K, jpiv(*), iwork(*), info
        real(ep), intent(out)   :: maxc2nrmk, relmaxc2nrmk, tau(*), work(*)
        call dgeqp3rk_quad(m, n, nrhs, kmax, abstol, reltol, A, lda, K, maxc2nrmk, relmaxc2nrmk, &
                           jpiv, tau, work, lwork, iwork, info)
    end subroutine dgeqp3rk

    subroutine zgeqp3rk(m, n, nrhs, kmax, abstol, reltol, A, lda, K, maxc2nrmk, relmaxc2nrmk, jpiv, &
                        tau, work, lwork, rwork, iwork, info)
        integer,     intent(in)    :: m, n, nrhs, kmax, lda, lwork
        real(ep),    intent(in)    :: abstol, reltol
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: K, jpiv(*), iwork(*), info
        real(ep),    intent(out)   :: maxc2nrmk, relmaxc2nrmk, rwork(*)
        complex(ep), intent(out)   :: tau(*), work(*)
        call zgeqp3rk_quad(m, n, nrhs, kmax, abstol, reltol, A, lda, K, maxc2nrmk, relmaxc2nrmk, &
                           jpiv, tau, work, lwork, rwork, iwork, info)
    end subroutine zgeqp3rk

    subroutine dbdsvdx(uplo, jobz, range, n, D, E, vl, vu, il, iu, ns, S, Z, ldz, work, iwork, info)
        character, intent(in)  :: uplo, jobz, range
        integer,   intent(in)  :: n, il, iu, ldz
        real(ep),  intent(in)  :: D(*), E(*), vl, vu
        integer,   intent(out) :: ns, iwork(*), info
        real(ep),  intent(out) :: S(*), Z(ldz,*), work(*)
        call dbdsvdx_quad(uplo, jobz, range, n, D, E, vl, vu, il, iu, ns, S, Z, ldz, work, iwork, &
                          info)
    end subroutine dbdsvdx

    subroutine dbbcsd(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, U1, ldu1, U2, ldu2, &
                      V1t, ldv1t, V2t, ldv2t, B11d, B11e, B12d, B12e, B21d, B21e, B22d, B22e, work, &
                      lwork, info)
        character, intent(in)    :: jobu1, jobu2, jobv1t, jobv2t, trans
        integer,   intent(in)    :: m, p, q, ldu1, ldu2, ldv1t, ldv2t, lwork
        real(ep),  intent(inout) :: theta(*), phi(*)
        real(ep),  intent(inout) :: U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*), V2t(ldv2t,*)
        real(ep),  intent(out)   :: B11d(*), B11e(*), B12d(*), B12e(*)
        real(ep),  intent(out)   :: B21d(*), B21e(*), B22d(*), B22e(*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dbbcsd_quad(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, U1, ldu1, U2, &
                         ldu2, V1t, ldv1t, V2t, ldv2t, B11d, B11e, B12d, B12e, B21d, B21e, B22d, &
                         B22e, work, lwork, info)
    end subroutine dbbcsd

    subroutine zbbcsd(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, U1, ldu1, U2, ldu2, &
                      V1t, ldv1t, V2t, ldv2t, B11d, B11e, B12d, B12e, B21d, B21e, B22d, B22e, rwork, &
                      lrwork, info)
        character,   intent(in)    :: jobu1, jobu2, jobv1t, jobv2t, trans
        integer,     intent(in)    :: m, p, q, ldu1, ldu2, ldv1t, ldv2t, lrwork
        real(ep),    intent(inout) :: theta(*), phi(*)
        complex(ep), intent(inout) :: U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*), V2t(ldv2t,*)
        real(ep),    intent(out)   :: B11d(*), B11e(*), B12d(*), B12e(*)
        real(ep),    intent(out)   :: B21d(*), B21e(*), B22d(*), B22e(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: info
        call zbbcsd_quad(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, U1, ldu1, U2, &
                         ldu2, V1t, ldv1t, V2t, ldv2t, B11d, B11e, B12d, B12e, B21d, B21e, B22d, &
                         B22e, rwork, lrwork, info)
    end subroutine zbbcsd

    subroutine dorbdb(trans, signs, m, p, q, X11, ldx11, X12, ldx12, X21, ldx21, X22, ldx22, theta, &
                      phi, taup1, taup2, tauq1, tauq2, work, lwork, info)
        character, intent(in)    :: trans, signs
        integer,   intent(in)    :: m, p, q, ldx11, ldx12, ldx21, ldx22, lwork
        real(ep),  intent(inout) :: X11(ldx11,*), X12(ldx12,*), X21(ldx21,*), X22(ldx22,*)
        real(ep),  intent(out)   :: theta(*), phi(*), taup1(*), taup2(*), tauq1(*), tauq2(*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dorbdb_quad(trans, signs, m, p, q, X11, ldx11, X12, ldx12, X21, ldx21, X22, ldx22, &
                         theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info)
    end subroutine dorbdb

    subroutine zunbdb(trans, signs, m, p, q, X11, ldx11, X12, ldx12, X21, ldx21, X22, ldx22, theta, &
                      phi, taup1, taup2, tauq1, tauq2, work, lwork, info)
        character,   intent(in)    :: trans, signs
        integer,     intent(in)    :: m, p, q, ldx11, ldx12, ldx21, ldx22, lwork
        complex(ep), intent(inout) :: X11(ldx11,*), X12(ldx12,*), X21(ldx21,*), X22(ldx22,*)
        real(ep),    intent(out)   :: theta(*), phi(*)
        complex(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), tauq2(*), work(*)
        integer,     intent(out)   :: info
        call zunbdb_quad(trans, signs, m, p, q, X11, ldx11, X12, ldx12, X21, ldx21, X22, ldx22, &
                         theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info)
    end subroutine zunbdb

    recursive subroutine dorcsd(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, &
                                X11, ldx11, X12, ldx12, X21, ldx21, X22, ldx22, theta, &
                                U1, ldu1, U2, ldu2, V1t, ldv1t, V2t, ldv2t, &
                                work, lwork, iwork, info)
        character, intent(in)    :: jobu1, jobu2, jobv1t, jobv2t, trans, signs
        integer,   intent(in)    :: m, p, q, ldx11, ldx12, ldx21, ldx22
        integer,   intent(in)    :: ldu1, ldu2, ldv1t, ldv2t, lwork
        real(ep),  intent(inout) :: X11(ldx11,*), X12(ldx12,*), X21(ldx21,*), X22(ldx22,*)
        real(ep),  intent(out)   :: theta(*)
        real(ep),  intent(out)   :: U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*), V2t(ldv2t,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: iwork(*), info
        call dorcsd_quad(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, &
                         X11, ldx11, X12, ldx12, X21, ldx21, X22, ldx22, theta, &
                         U1, ldu1, U2, ldu2, V1t, ldv1t, V2t, ldv2t, &
                         work, lwork, iwork, info)
    end subroutine dorcsd

    recursive subroutine zuncsd(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, &
                                X11, ldx11, X12, ldx12, X21, ldx21, X22, ldx22, theta, &
                                U1, ldu1, U2, ldu2, V1t, ldv1t, V2t, ldv2t, &
                                work, lwork, rwork, lrwork, iwork, info)
        character,   intent(in)    :: jobu1, jobu2, jobv1t, jobv2t, trans, signs
        integer,     intent(in)    :: m, p, q, ldx11, ldx12, ldx21, ldx22
        integer,     intent(in)    :: ldu1, ldu2, ldv1t, ldv2t, lwork, lrwork
        complex(ep), intent(inout) :: X11(ldx11,*), X12(ldx12,*), X21(ldx21,*), X22(ldx22,*)
        real(ep),    intent(out)   :: theta(*)
        complex(ep), intent(out)   :: U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*), V2t(ldv2t,*)
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: iwork(*), info
        call zuncsd_quad(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, &
                         X11, ldx11, X12, ldx12, X21, ldx21, X22, ldx22, theta, &
                         U1, ldu1, U2, ldu2, V1t, ldv1t, V2t, ldv2t, &
                         work, lwork, rwork, lrwork, iwork, info)
    end subroutine zuncsd

    subroutine dorcsd2by1(jobu1, jobu2, jobv1t, m, p, q, X11, ldx11, X21, ldx21, theta, U1, ldu1, U2, &
                          ldu2, V1t, ldv1t, work, lwork, iwork, info)
        character, intent(in)    :: jobu1, jobu2, jobv1t
        integer,   intent(in)    :: m, p, q, ldx11, ldx21, ldu1, ldu2, ldv1t, lwork
        real(ep),  intent(inout) :: X11(ldx11,*), X21(ldx21,*)
        real(ep),  intent(out)   :: theta(*), U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: iwork(*), info
        call dorcsd2by1_quad(jobu1, jobu2, jobv1t, m, p, q, X11, ldx11, X21, ldx21, theta, U1, ldu1, &
                             U2, ldu2, V1t, ldv1t, work, lwork, iwork, info)
    end subroutine dorcsd2by1

    subroutine zuncsd2by1(jobu1, jobu2, jobv1t, m, p, q, X11, ldx11, X21, ldx21, theta, U1, ldu1, U2, &
                          ldu2, V1t, ldv1t, work, lwork, rwork, lrwork, iwork, info)
        character,   intent(in)    :: jobu1, jobu2, jobv1t
        integer,     intent(in)    :: m, p, q, ldx11, ldx21, ldu1, ldu2, ldv1t, lwork, lrwork
        complex(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
        real(ep),    intent(out)   :: theta(*)
        complex(ep), intent(out)   :: U1(ldu1,*), U2(ldu2,*), V1t(ldv1t,*)
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: iwork(*), info
        call zuncsd2by1_quad(jobu1, jobu2, jobv1t, m, p, q, X11, ldx11, X21, ldx21, theta, U1, ldu1, &
                             U2, ldu2, V1t, ldv1t, work, lwork, rwork, lrwork, iwork, info)
    end subroutine zuncsd2by1

    subroutine dggsvd3(jobu, jobv, jobq, m, n, p, k, l, A, lda, B, ldb, alpha, beta, U, ldu, V, ldv, &
                       Q, ldq, work, lwork, iwork, info)
        character, intent(in)    :: jobu, jobv, jobq
        integer,   intent(in)    :: m, n, p, lda, ldb, ldu, ldv, ldq, lwork
        integer,   intent(out)   :: k, l
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: alpha(*), beta(*)
        real(ep),  intent(out)   :: U(ldu,*), V(ldv,*), Q(ldq,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dggsvd3_quad(jobu, jobv, jobq, m, n, p, k, l, A, lda, B, ldb, alpha, beta, U, ldu, V, &
                          ldv, Q, ldq, work, lwork, iwork, info)
    end subroutine dggsvd3

    subroutine zggsvd3(jobu, jobv, jobq, m, n, p, k, l, A, lda, B, ldb, alpha, beta, U, ldu, V, ldv, &
                       Q, ldq, work, lwork, rwork, iwork, info)
        character,   intent(in)    :: jobu, jobv, jobq
        integer,     intent(in)    :: m, n, p, lda, ldb, ldu, ldv, ldq, lwork
        integer,     intent(out)   :: k, l
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),    intent(out)   :: alpha(*), beta(*)
        complex(ep), intent(out)   :: U(ldu,*), V(ldv,*), Q(ldq,*), work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: iwork(*), info
        call zggsvd3_quad(jobu, jobv, jobq, m, n, p, k, l, A, lda, B, ldb, alpha, beta, U, ldu, V, &
                          ldv, Q, ldq, work, lwork, rwork, iwork, info)
    end subroutine zggsvd3

    subroutine dggsvp3(jobu, jobv, jobq, m, p, n, A, lda, B, ldb, tola, tolb, k, l, U, ldu, V, ldv, &
                       Q, ldq, iwork, tau, work, lwork, info)
        character, intent(in)    :: jobu, jobv, jobq
        integer,   intent(in)    :: m, p, n, lda, ldb, ldu, ldv, ldq, lwork
        real(ep),  intent(in)    :: tola, tolb
        integer,   intent(out)   :: k, l
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: U(ldu,*), V(ldv,*), Q(ldq,*), tau(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dggsvp3_quad(jobu, jobv, jobq, m, p, n, A, lda, B, ldb, tola, tolb, k, l, U, ldu, V, &
                          ldv, Q, ldq, iwork, tau, work, lwork, info)
    end subroutine dggsvp3

    subroutine zggsvp3(jobu, jobv, jobq, m, p, n, A, lda, B, ldb, tola, tolb, k, l, U, ldu, V, ldv, &
                       Q, ldq, iwork, rwork, tau, work, lwork, info)
        character,   intent(in)    :: jobu, jobv, jobq
        integer,     intent(in)    :: m, p, n, lda, ldb, ldu, ldv, ldq, lwork
        real(ep),    intent(in)    :: tola, tolb
        integer,     intent(out)   :: k, l
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: U(ldu,*), V(ldv,*), Q(ldq,*), tau(*), work(*)
        real(ep),    intent(out)   :: rwork(*)
        integer,     intent(out)   :: iwork(*), info
        call zggsvp3_quad(jobu, jobv, jobq, m, p, n, A, lda, B, ldb, tola, tolb, k, l, U, ldu, V, &
                          ldv, Q, ldq, iwork, rwork, tau, work, lwork, info)
    end subroutine zggsvp3

    subroutine dorbdb1(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, &
                       info)
        integer,  intent(in)    :: m, p, q, ldx11, ldx21, lwork
        real(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
        real(ep), intent(out)   :: theta(*), phi(*)
        real(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
        integer,  intent(out)   :: info
        call dorbdb1_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, &
                          lwork, info)
    end subroutine dorbdb1

    subroutine dorbdb2(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, &
                       info)
        integer,  intent(in)    :: m, p, q, ldx11, ldx21, lwork
        real(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
        real(ep), intent(out)   :: theta(*), phi(*)
        real(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
        integer,  intent(out)   :: info
        call dorbdb2_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, &
                          lwork, info)
    end subroutine dorbdb2

    subroutine dorbdb3(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, &
                       info)
        integer,  intent(in)    :: m, p, q, ldx11, ldx21, lwork
        real(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
        real(ep), intent(out)   :: theta(*), phi(*)
        real(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
        integer,  intent(out)   :: info
        call dorbdb3_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, &
                          lwork, info)
    end subroutine dorbdb3

    subroutine dorbdb4(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, &
                       work, lwork, info)
        integer,  intent(in)    :: m, p, q, ldx11, ldx21, lwork
        real(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
        real(ep), intent(out)   :: theta(*), phi(*)
        real(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), phantom(*), work(*)
        integer,  intent(out)   :: info
        call dorbdb4_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, &
                          work, lwork, info)
    end subroutine dorbdb4

    subroutine dorbdb5(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, work, lwork, info)
        integer,  intent(in)    :: m1, m2, n, incx1, incx2, ldq1, ldq2, lwork
        real(ep), intent(inout) :: X1(*), X2(*)
        real(ep), intent(in)    :: Q1(ldq1,*), Q2(ldq2,*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dorbdb5_quad(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, work, lwork, info)
    end subroutine dorbdb5

    subroutine dorbdb6(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, work, lwork, info)
        integer,  intent(in)    :: m1, m2, n, incx1, incx2, ldq1, ldq2, lwork
        real(ep), intent(inout) :: X1(*), X2(*)
        real(ep), intent(in)    :: Q1(ldq1,*), Q2(ldq2,*)
        real(ep), intent(out)   :: work(*)
        integer,  intent(out)   :: info
        call dorbdb6_quad(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, work, lwork, info)
    end subroutine dorbdb6

    subroutine zunbdb1(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, &
                       info)
        integer,     intent(in)    :: m, p, q, ldx11, ldx21, lwork
        complex(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
        real(ep),    intent(out)   :: theta(*), phi(*)
        complex(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
        integer,     intent(out)   :: info
        call zunbdb1_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, &
                          lwork, info)
    end subroutine zunbdb1

    subroutine zunbdb2(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, &
                       info)
        integer,     intent(in)    :: m, p, q, ldx11, ldx21, lwork
        complex(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
        real(ep),    intent(out)   :: theta(*), phi(*)
        complex(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
        integer,     intent(out)   :: info
        call zunbdb2_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, &
                          lwork, info)
    end subroutine zunbdb2

    subroutine zunbdb3(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, &
                       info)
        integer,     intent(in)    :: m, p, q, ldx11, ldx21, lwork
        complex(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
        real(ep),    intent(out)   :: theta(*), phi(*)
        complex(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), work(*)
        integer,     intent(out)   :: info
        call zunbdb3_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, work, &
                          lwork, info)
    end subroutine zunbdb3

    subroutine zunbdb4(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, &
                       work, lwork, info)
        integer,     intent(in)    :: m, p, q, ldx11, ldx21, lwork
        complex(ep), intent(inout) :: X11(ldx11,*), X21(ldx21,*)
        real(ep),    intent(out)   :: theta(*), phi(*)
        complex(ep), intent(out)   :: taup1(*), taup2(*), tauq1(*), phantom(*), work(*)
        integer,     intent(out)   :: info
        call zunbdb4_quad(m, p, q, X11, ldx11, X21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, &
                          work, lwork, info)
    end subroutine zunbdb4

    subroutine zunbdb5(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, work, lwork, info)
        integer,     intent(in)    :: m1, m2, n, incx1, incx2, ldq1, ldq2, lwork
        complex(ep), intent(inout) :: X1(*), X2(*)
        complex(ep), intent(in)    :: Q1(ldq1,*), Q2(ldq2,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunbdb5_quad(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, work, lwork, info)
    end subroutine zunbdb5

    subroutine zunbdb6(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, work, lwork, info)
        integer,     intent(in)    :: m1, m2, n, incx1, incx2, ldq1, ldq2, lwork
        complex(ep), intent(inout) :: X1(*), X2(*)
        complex(ep), intent(in)    :: Q1(ldq1,*), Q2(ldq2,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunbdb6_quad(m1, m2, n, X1, incx1, X2, incx2, Q1, ldq1, Q2, ldq2, work, lwork, info)
    end subroutine zunbdb6

    subroutine dorm22(side, trans, m, n, n1, n2, Q, ldq, C, ldc, work, lwork, info)
        character, intent(in)    :: side, trans
        integer,   intent(in)    :: m, n, n1, n2, ldq, ldc, lwork
        real(ep),  intent(in)    :: Q(ldq,*)
        real(ep),  intent(inout) :: C(ldc,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dorm22_quad(side, trans, m, n, n1, n2, Q, ldq, C, ldc, work, lwork, info)
    end subroutine dorm22

    subroutine zunm22(side, trans, m, n, n1, n2, Q, ldq, C, ldc, work, lwork, info)
        character,   intent(in)    :: side, trans
        integer,     intent(in)    :: m, n, n1, n2, ldq, ldc, lwork
        complex(ep), intent(in)    :: Q(ldq,*)
        complex(ep), intent(inout) :: C(ldc,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zunm22_quad(side, trans, m, n, n1, n2, Q, ldq, C, ldc, work, lwork, info)
    end subroutine zunm22

    subroutine dgesvdq(joba, jobp, jobr, jobu, jobv, m, n, A, lda, S, U, ldu, V, ldv, numrank, iwork, &
                       liwork, work, lwork, rwork, lrwork, info)
        character, intent(in)    :: joba, jobp, jobr, jobu, jobv
        integer,   intent(in)    :: m, n, lda, ldu, ldv, liwork, lwork, lrwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: S(*), U(ldu,*), V(ldv,*), work(*), rwork(*)
        integer,   intent(out)   :: numrank, iwork(*), info
        call dgesvdq_quad(joba, jobp, jobr, jobu, jobv, m, n, A, lda, S, U, ldu, V, ldv, numrank, &
                          iwork, liwork, work, lwork, rwork, lrwork, info)
    end subroutine dgesvdq

    subroutine zgesvdq(joba, jobp, jobr, jobu, jobv, m, n, A, lda, S, U, ldu, V, ldv, numrank, iwork, &
                       liwork, cwork, lcwork, rwork, lrwork, info)
        character,   intent(in)    :: joba, jobp, jobr, jobu, jobv
        integer,     intent(in)    :: m, n, lda, ldu, ldv, liwork, lcwork, lrwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: S(*), rwork(*)
        complex(ep), intent(out)   :: U(ldu,*), V(ldv,*), cwork(*)
        integer,     intent(out)   :: numrank, iwork(*), info
        call zgesvdq_quad(joba, jobp, jobr, jobu, jobv, m, n, A, lda, S, U, ldu, V, ldv, numrank, &
                          iwork, liwork, cwork, lcwork, rwork, lrwork, info)
    end subroutine zgesvdq

    subroutine dgesvdx(jobu, jobvt, range, m, n, A, lda, vl, vu, il, iu, ns, S, U, ldu, Vt, ldvt, &
                       work, lwork, iwork, info)
        character, intent(in)    :: jobu, jobvt, range
        integer,   intent(in)    :: m, n, il, iu, lda, ldu, ldvt, lwork
        real(ep),  intent(in)    :: vl, vu
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: ns, iwork(*), info
        real(ep),  intent(out)   :: S(*), U(ldu,*), Vt(ldvt,*), work(*)
        call dgesvdx_quad(jobu, jobvt, range, m, n, A, lda, vl, vu, il, iu, ns, S, U, ldu, Vt, ldvt, &
                          work, lwork, iwork, info)
    end subroutine dgesvdx

    subroutine zgesvdx(jobu, jobvt, range, m, n, A, lda, vl, vu, il, iu, ns, S, U, ldu, Vt, ldvt, &
                       work, lwork, rwork, iwork, info)
        character,   intent(in)    :: jobu, jobvt, range
        integer,     intent(in)    :: m, n, il, iu, lda, ldu, ldvt, lwork
        real(ep),    intent(in)    :: vl, vu
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ns, iwork(*), info
        real(ep),    intent(out)   :: S(*), rwork(*)
        complex(ep), intent(out)   :: U(ldu,*), Vt(ldvt,*), work(*)
        call zgesvdx_quad(jobu, jobvt, range, m, n, A, lda, vl, vu, il, iu, ns, S, U, ldu, Vt, ldvt, &
                          work, lwork, rwork, iwork, info)
    end subroutine zgesvdx

    subroutine dgbbrd(vect, m, n, ncc, kl, ku, AB, ldab, D, E, Q, ldq, Pt, ldpt, C, ldc, work, info)
        character, intent(in)    :: vect
        integer,   intent(in)    :: m, n, ncc, kl, ku, ldab, ldq, ldpt, ldc
        real(ep),  intent(inout) :: AB(ldab,*), C(ldc,*)
        real(ep),  intent(out)   :: D(*), E(*), Q(ldq,*), Pt(ldpt,*), work(*)
        integer,   intent(out)   :: info
        call dgbbrd_quad(vect, m, n, ncc, kl, ku, AB, ldab, D, E, Q, ldq, Pt, ldpt, C, ldc, work, &
                         info)
    end subroutine dgbbrd

    subroutine zgbbrd(vect, m, n, ncc, kl, ku, AB, ldab, D, E, Q, ldq, Pt, ldpt, C, ldc, work, rwork, &
                      info)
        character,   intent(in)    :: vect
        integer,     intent(in)    :: m, n, ncc, kl, ku, ldab, ldq, ldpt, ldc
        complex(ep), intent(inout) :: AB(ldab,*), C(ldc,*)
        real(ep),    intent(out)   :: D(*), E(*), rwork(*)
        complex(ep), intent(out)   :: Q(ldq,*), Pt(ldpt,*), work(*)
        integer,     intent(out)   :: info
        call zgbbrd_quad(vect, m, n, ncc, kl, ku, AB, ldab, D, E, Q, ldq, Pt, ldpt, C, ldc, work, &
                         rwork, info)
    end subroutine zgbbrd

    subroutine dstegr(jobz, range, n, D, E, vl, vu, il, iu, abstol, m, W, Z, ldz, isuppz, work, &
                      lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, range
        integer,   intent(in)    :: n, il, iu, ldz, lwork, liwork
        real(ep),  intent(in)    :: vl, vu, abstol
        real(ep),  intent(inout) :: D(*), E(*)
        integer,   intent(out)   :: m, isuppz(*), iwork(*), info
        real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
        call dstegr_quad(jobz, range, n, D, E, vl, vu, il, iu, abstol, m, W, Z, ldz, isuppz, work, &
                         lwork, iwork, liwork, info)
    end subroutine dstegr

    subroutine zstegr(jobz, range, n, D, E, vl, vu, il, iu, abstol, m, W, Z, ldz, isuppz, work, &
                      lwork, iwork, liwork, info)
        character,   intent(in)    :: jobz, range
        integer,     intent(in)    :: n, il, iu, ldz, lwork, liwork
        real(ep),    intent(in)    :: vl, vu, abstol
        real(ep),    intent(inout) :: D(*), E(*)
        integer,     intent(out)   :: m, isuppz(*), iwork(*), info
        real(ep),    intent(out)   :: W(*), work(*)
        complex(ep), intent(out)   :: Z(ldz,*)
        call zstegr_quad(jobz, range, n, D, E, vl, vu, il, iu, abstol, m, W, Z, ldz, isuppz, work, &
                         lwork, iwork, liwork, info)
    end subroutine zstegr

    subroutine dstemr(jobz, range, n, D, E, vl, vu, il, iu, m, W, Z, ldz, nzc, isuppz, tryrac, work, &
                      lwork, iwork, liwork, info)
        character, intent(in)    :: jobz, range
        integer,   intent(in)    :: n, il, iu, ldz, nzc, lwork, liwork
        logical,   intent(inout) :: tryrac
        real(ep),  intent(in)    :: vl, vu
        real(ep),  intent(inout) :: D(*), E(*)
        integer,   intent(out)   :: m, isuppz(*), iwork(*), info
        real(ep),  intent(out)   :: W(*), Z(ldz,*), work(*)
        call dstemr_quad(jobz, range, n, D, E, vl, vu, il, iu, m, W, Z, ldz, nzc, isuppz, tryrac, &
                         work, lwork, iwork, liwork, info)
    end subroutine dstemr

    subroutine zstemr(jobz, range, n, D, E, vl, vu, il, iu, m, W, Z, ldz, nzc, isuppz, tryrac, work, &
                      lwork, iwork, liwork, info)
        character,   intent(in)    :: jobz, range
        integer,     intent(in)    :: n, il, iu, ldz, nzc, lwork, liwork
        logical,     intent(inout) :: tryrac
        real(ep),    intent(in)    :: vl, vu
        real(ep),    intent(inout) :: D(*), E(*)
        integer,     intent(out)   :: m, isuppz(*), iwork(*), info
        real(ep),    intent(out)   :: W(*), work(*)
        complex(ep), intent(out)   :: Z(ldz,*)
        call zstemr_quad(jobz, range, n, D, E, vl, vu, il, iu, m, W, Z, ldz, nzc, isuppz, tryrac, &
                         work, lwork, iwork, liwork, info)
    end subroutine zstemr

    subroutine dstein(n, D, E, m, W, iblock, isplit, Z, ldz, work, iwork, ifail, info)
        integer,  intent(in)    :: n, m, ldz
        real(ep), intent(in)    :: D(*), E(*), W(*)
        integer,  intent(in)    :: iblock(*), isplit(*)
        real(ep), intent(out)   :: Z(ldz,*), work(*)
        integer,  intent(out)   :: iwork(*), ifail(*), info
        call dstein_quad(n, D, E, m, W, iblock, isplit, Z, ldz, work, iwork, ifail, info)
    end subroutine dstein

    subroutine zstein(n, D, E, m, W, iblock, isplit, Z, ldz, work, iwork, ifail, info)
        integer,     intent(in)    :: n, m, ldz
        real(ep),    intent(in)    :: D(*), E(*), W(*)
        integer,     intent(in)    :: iblock(*), isplit(*)
        complex(ep), intent(out)   :: Z(ldz,*)
        real(ep),    intent(out)   :: work(*)
        integer,     intent(out)   :: iwork(*), ifail(*), info
        call zstein_quad(n, D, E, m, W, iblock, isplit, Z, ldz, work, iwork, ifail, info)
    end subroutine zstein

    subroutine dsytrf_rook(uplo, n, A, lda, ipiv, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: ipiv(*), info
        real(ep),  intent(out)   :: work(*)
        call dsytrf_rook_quad(uplo, n, A, lda, ipiv, work, lwork, info)
    end subroutine dsytrf_rook

    subroutine dsytf2_rook(uplo, n, A, lda, ipiv, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: ipiv(*), info
        call dsytf2_rook_quad(uplo, n, A, lda, ipiv, info)
    end subroutine dsytf2_rook

    subroutine dsytrs_rook(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dsytrs_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine dsytrs_rook

    subroutine dsytri_rook(uplo, n, A, lda, ipiv, work, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dsytri_rook_quad(uplo, n, A, lda, ipiv, work, info)
    end subroutine dsytri_rook

    subroutine dsycon_rook(uplo, n, A, lda, ipiv, anorm, rcond, work, iwork, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n, lda
        real(ep),  intent(in)  :: A(lda,*), anorm
        integer,   intent(in)  :: ipiv(*)
        real(ep),  intent(out) :: rcond
        real(ep),  intent(out) :: work(*)
        integer,   intent(out) :: iwork(*), info
        call dsycon_rook_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, iwork, info)
    end subroutine dsycon_rook

    subroutine dsysv_rook(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        integer,   intent(out)   :: ipiv(*), info
        real(ep),  intent(out)   :: work(*)
        call dsysv_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine dsysv_rook

    subroutine zhetrf_rook(uplo, n, A, lda, ipiv, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zhetrf_rook_quad(uplo, n, A, lda, ipiv, work, lwork, info)
    end subroutine zhetrf_rook

    subroutine zhetf2_rook(uplo, n, A, lda, ipiv, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ipiv(*), info
        call zhetf2_rook_quad(uplo, n, A, lda, ipiv, info)
    end subroutine zhetf2_rook

    subroutine zhetrs_rook(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(in)    :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zhetrs_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine zhetrs_rook

    subroutine zhetri_rook(uplo, n, A, lda, ipiv, work, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zhetri_rook_quad(uplo, n, A, lda, ipiv, work, info)
    end subroutine zhetri_rook

    subroutine zhecon_rook(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        integer,     intent(in)  :: ipiv(*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zhecon_rook_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
    end subroutine zhecon_rook

    subroutine zhesv_rook(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zhesv_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine zhesv_rook

    subroutine dsytrf_rk(uplo, n, A, lda, e, ipiv, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: e(*), work(*)
        integer,   intent(out)   :: ipiv(*), info
        call dsytrf_rk_quad(uplo, n, A, lda, e, ipiv, work, lwork, info)
    end subroutine dsytrf_rk

    subroutine dsytf2_rk(uplo, n, A, lda, e, ipiv, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: e(*)
        integer,   intent(out)   :: ipiv(*), info
        call dsytf2_rk_quad(uplo, n, A, lda, e, ipiv, info)
    end subroutine dsytf2_rk

    subroutine dsysv_rk(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        real(ep),  intent(out)   :: e(*), work(*)
        integer,   intent(out)   :: ipiv(*), info
        call dsysv_rk_quad(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, work, lwork, info)
    end subroutine dsysv_rk

    subroutine dsytri_3x(uplo, n, A, lda, e, ipiv, work, nb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, nb
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(in)    :: e(*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(out)   :: work(n+nb+1, nb+3)
        integer,   intent(out)   :: info
        call dsytri_3x_quad(uplo, n, A, lda, e, ipiv, work, nb, info)
    end subroutine dsytri_3x

    subroutine dsytri_3(uplo, n, A, lda, e, ipiv, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(in)    :: e(*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dsytri_3_quad(uplo, n, A, lda, e, ipiv, work, lwork, info)
    end subroutine dsytri_3

    subroutine zhetrf_rk(uplo, n, A, lda, e, ipiv, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: e(*), work(*)
        integer,     intent(out)   :: ipiv(*), info
        call zhetrf_rk_quad(uplo, n, A, lda, e, ipiv, work, lwork, info)
    end subroutine zhetrf_rk

    subroutine zhetf2_rk(uplo, n, A, lda, e, ipiv, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: e(*)
        integer,     intent(out)   :: ipiv(*), info
        call zhetf2_rk_quad(uplo, n, A, lda, e, ipiv, info)
    end subroutine zhetf2_rk

    subroutine zhesv_rk(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: e(*), work(*)
        integer,     intent(out)   :: ipiv(*), info
        call zhesv_rk_quad(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, work, lwork, info)
    end subroutine zhesv_rk

    subroutine zhetri_3x(uplo, n, A, lda, e, ipiv, work, nb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, nb
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: e(*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(n+nb+1, nb+3)
        integer,     intent(out)   :: info
        call zhetri_3x_quad(uplo, n, A, lda, e, ipiv, work, nb, info)
    end subroutine zhetri_3x

    subroutine zsytrf_rk(uplo, n, A, lda, e, ipiv, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: e(*), work(*)
        integer,     intent(out)   :: ipiv(*), info
        call zsytrf_rk_quad(uplo, n, A, lda, e, ipiv, work, lwork, info)
    end subroutine zsytrf_rk

    subroutine zsytf2_rk(uplo, n, A, lda, e, ipiv, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: e(*)
        integer,     intent(out)   :: ipiv(*), info
        call zsytf2_rk_quad(uplo, n, A, lda, e, ipiv, info)
    end subroutine zsytf2_rk

    subroutine zsysv_rk(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        complex(ep), intent(out)   :: e(*), work(*)
        integer,     intent(out)   :: ipiv(*), info
        call zsysv_rk_quad(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, work, lwork, info)
    end subroutine zsysv_rk

    subroutine zsytri_3x(uplo, n, A, lda, e, ipiv, work, nb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, nb
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: e(*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(n+nb+1, nb+3)
        integer,     intent(out)   :: info
        call zsytri_3x_quad(uplo, n, A, lda, e, ipiv, work, nb, info)
    end subroutine zsytri_3x

    subroutine zsytri_3(uplo, n, A, lda, e, ipiv, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(in)    :: e(*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zsytri_3_quad(uplo, n, A, lda, e, ipiv, work, lwork, info)
    end subroutine zsytri_3

    subroutine dsytri2(uplo, n, A, lda, ipiv, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dsytri2_quad(uplo, n, A, lda, ipiv, work, lwork, info)
    end subroutine dsytri2

    subroutine zsytri2(uplo, n, A, lda, ipiv, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zsytri2_quad(uplo, n, A, lda, ipiv, work, lwork, info)
    end subroutine zsytri2

    subroutine dsytrs2(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dsytrs2_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, info)
    end subroutine dsytrs2

    subroutine zsytrs2(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(in)    :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zsytrs2_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, info)
    end subroutine zsytrs2

    subroutine dsytrs_3(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb
        real(ep),  intent(in)    :: A(lda,*), e(*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dsytrs_3_quad(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, info)
    end subroutine dsytrs_3

    subroutine zsytrs_3(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(in)    :: A(lda,*), e(*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zsytrs_3_quad(uplo, n, nrhs, A, lda, e, ipiv, B, ldb, info)
    end subroutine zsytrs_3

    subroutine dsycon_3(uplo, n, A, lda, e, ipiv, anorm, rcond, work, iwork, info)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n, lda
        real(ep),  intent(in)  :: A(lda,*), e(*), anorm
        integer,   intent(in)  :: ipiv(*)
        real(ep),  intent(out) :: rcond, work(*)
        integer,   intent(out) :: iwork(*), info
        call dsycon_3_quad(uplo, n, A, lda, e, ipiv, anorm, rcond, work, iwork, info)
    end subroutine dsycon_3

    subroutine zsycon_3(uplo, n, A, lda, e, ipiv, anorm, rcond, work, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*), e(*)
        real(ep),    intent(in)  :: anorm
        integer,     intent(in)  :: ipiv(*)
        real(ep),    intent(out) :: rcond
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zsycon_3_quad(uplo, n, A, lda, e, ipiv, anorm, rcond, work, info)
    end subroutine zsycon_3

    subroutine dsytrf_aa(uplo, n, A, lda, ipiv, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, lwork
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(out)   :: ipiv(*), info
        real(ep),  intent(out)   :: work(*)
        call dsytrf_aa_quad(uplo, n, A, lda, ipiv, work, lwork, info)
    end subroutine dsytrf_aa

    subroutine dsytrs_aa(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb, lwork
        real(ep),  intent(in)    :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(inout) :: B(ldb,*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dsytrs_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine dsytrs_aa

    subroutine dsysv_aa(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*)
        integer,   intent(out)   :: ipiv(*), info
        real(ep),  intent(out)   :: work(*)
        call dsysv_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine dsysv_aa

    subroutine dsytrf_aa_2stage(uplo, n, A, lda, TB, ltb, ipiv, ipiv2, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, ltb, lwork
        real(ep),  intent(inout) :: A(lda,*), TB(*)
        integer,   intent(out)   :: ipiv(*), ipiv2(*), info
        real(ep),  intent(out)   :: work(*)
        call dsytrf_aa_2stage_quad(uplo, n, A, lda, TB, ltb, ipiv, ipiv2, work, lwork, info)
    end subroutine dsytrf_aa_2stage

    subroutine dsytrs_aa_2stage(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ltb, ldb
        real(ep),  intent(in)    :: A(lda,*), TB(*)
        integer,   intent(in)    :: ipiv(*), ipiv2(*)
        real(ep),  intent(inout) :: B(ldb,*)
        integer,   intent(out)   :: info
        call dsytrs_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, info)
    end subroutine dsytrs_aa_2stage

    subroutine dsysv_aa_2stage(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, work, lwork, &
                               info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, nrhs, lda, ltb, ldb, lwork
        real(ep),  intent(inout) :: A(lda,*), B(ldb,*), TB(*)
        integer,   intent(out)   :: ipiv(*), ipiv2(*), info
        real(ep),  intent(out)   :: work(*)
        call dsysv_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, work, lwork, &
                                  info)
    end subroutine dsysv_aa_2stage

    subroutine zhetrf_aa(uplo, n, A, lda, ipiv, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zhetrf_aa_quad(uplo, n, A, lda, ipiv, work, lwork, info)
    end subroutine zhetrf_aa

    subroutine zhetrs_aa(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
        complex(ep), intent(in)    :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zhetrs_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine zhetrs_aa

    subroutine zhesv_aa(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zhesv_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine zhesv_aa

    subroutine zhetrf_aa_2stage(uplo, n, A, lda, TB, ltb, ipiv, ipiv2, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, ltb, lwork
        complex(ep), intent(inout) :: A(lda,*), TB(*)
        integer,     intent(out)   :: ipiv(*), ipiv2(*), info
        complex(ep), intent(out)   :: work(*)
        call zhetrf_aa_2stage_quad(uplo, n, A, lda, TB, ltb, ipiv, ipiv2, work, lwork, info)
    end subroutine zhetrf_aa_2stage

    subroutine zhetrs_aa_2stage(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ltb, ldb
        complex(ep), intent(in)    :: A(lda,*), TB(*)
        integer,     intent(in)    :: ipiv(*), ipiv2(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zhetrs_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, info)
    end subroutine zhetrs_aa_2stage

    subroutine zhesv_aa_2stage(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, work, lwork, &
                               info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ltb, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*), TB(*)
        integer,     intent(out)   :: ipiv(*), ipiv2(*), info
        complex(ep), intent(out)   :: work(*)
        call zhesv_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, work, lwork, &
                                  info)
    end subroutine zhesv_aa_2stage

    subroutine zsytrf_aa(uplo, n, A, lda, ipiv, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zsytrf_aa_quad(uplo, n, A, lda, ipiv, work, lwork, info)
    end subroutine zsytrf_aa

    subroutine zsytrs_aa(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
        complex(ep), intent(in)    :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zsytrs_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine zsytrs_aa

    subroutine zsysv_aa(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zsysv_aa_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine zsysv_aa

    subroutine zsytrf_aa_2stage(uplo, n, A, lda, TB, ltb, ipiv, ipiv2, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, ltb, lwork
        complex(ep), intent(inout) :: A(lda,*), TB(*)
        integer,     intent(out)   :: ipiv(*), ipiv2(*), info
        complex(ep), intent(out)   :: work(*)
        call zsytrf_aa_2stage_quad(uplo, n, A, lda, TB, ltb, ipiv, ipiv2, work, lwork, info)
    end subroutine zsytrf_aa_2stage

    subroutine zsytrs_aa_2stage(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ltb, ldb
        complex(ep), intent(in)    :: A(lda,*), TB(*)
        integer,     intent(in)    :: ipiv(*), ipiv2(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zsytrs_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, info)
    end subroutine zsytrs_aa_2stage

    subroutine zsysv_aa_2stage(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, work, lwork, &
                               info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ltb, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*), TB(*)
        integer,     intent(out)   :: ipiv(*), ipiv2(*), info
        complex(ep), intent(out)   :: work(*)
        call zsysv_aa_2stage_quad(uplo, n, nrhs, A, lda, TB, ltb, ipiv, ipiv2, B, ldb, work, lwork, &
                                  info)
    end subroutine zsysv_aa_2stage

    subroutine zsytrf(uplo, n, A, lda, ipiv, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zsytrf_quad(uplo, n, A, lda, ipiv, work, lwork, info)
    end subroutine zsytrf

    subroutine zsytrs(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(in)    :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zsytrs_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine zsytrs

    subroutine zsysv(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zsysv_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine zsysv

    subroutine zsycon(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        integer,     intent(in)  :: ipiv(*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zsycon_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
    end subroutine zsycon

    subroutine zsytri(uplo, n, A, lda, ipiv, work, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zsytri_quad(uplo, n, A, lda, ipiv, work, info)
    end subroutine zsytri

    subroutine zsytrf_rook(uplo, n, A, lda, ipiv, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, lwork
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zsytrf_rook_quad(uplo, n, A, lda, ipiv, work, lwork, info)
    end subroutine zsytrf_rook

    subroutine zsytf2_rook(uplo, n, A, lda, ipiv, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(out)   :: ipiv(*), info
        call zsytf2_rook_quad(uplo, n, A, lda, ipiv, info)
    end subroutine zsytf2_rook

    subroutine zsytrs_rook(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb
        complex(ep), intent(in)    :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: B(ldb,*)
        integer,     intent(out)   :: info
        call zsytrs_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, info)
    end subroutine zsytrs_rook

    subroutine zsytri_rook(uplo, n, A, lda, ipiv, work, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zsytri_rook_quad(uplo, n, A, lda, ipiv, work, info)
    end subroutine zsytri_rook

    subroutine zsycon_rook(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n, lda
        complex(ep), intent(in)  :: A(lda,*)
        integer,     intent(in)  :: ipiv(*)
        real(ep),    intent(in)  :: anorm
        real(ep),    intent(out) :: rcond
        complex(ep), intent(out) :: work(*)
        integer,     intent(out) :: info
        call zsycon_rook_quad(uplo, n, A, lda, ipiv, anorm, rcond, work, info)
    end subroutine zsycon_rook

    subroutine zsysv_rook(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, lda, ldb, lwork
        complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
        integer,     intent(out)   :: ipiv(*), info
        complex(ep), intent(out)   :: work(*)
        call zsysv_rook_quad(uplo, n, nrhs, A, lda, ipiv, B, ldb, work, lwork, info)
    end subroutine zsysv_rook

    subroutine dsytri(uplo, n, A, lda, ipiv, work, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(out)   :: work(*)
        integer,   intent(out)   :: info
        call dsytri_quad(uplo, n, A, lda, ipiv, work, info)
    end subroutine dsytri

    subroutine zhetri(uplo, n, A, lda, ipiv, work, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zhetri_quad(uplo, n, A, lda, ipiv, work, info)
    end subroutine zhetri

    subroutine dsytri2x(uplo, n, A, lda, ipiv, work, nb, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, nb
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(out)   :: work(n+nb+1, nb+3)
        integer,   intent(out)   :: info
        call dsytri2x_quad(uplo, n, A, lda, ipiv, work, nb, info)
    end subroutine dsytri2x

    subroutine zhetri2x(uplo, n, A, lda, ipiv, work, nb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, nb
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(n+nb+1, nb+3)
        integer,     intent(out)   :: info
        call zhetri2x_quad(uplo, n, A, lda, ipiv, work, nb, info)
    end subroutine zhetri2x

    subroutine zsytri2x(uplo, n, A, lda, ipiv, work, nb, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, nb
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: work(n+nb+1, nb+3)
        integer,     intent(out)   :: info
        call zsytri2x_quad(uplo, n, A, lda, ipiv, work, nb, info)
    end subroutine zsytri2x

    subroutine zsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta, A(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        call zsymv_quad(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine zsymv

    subroutine zsyr(uplo, n, alpha, x, incx, A, lda)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, incx
        complex(ep), intent(in)    :: alpha, x(*)
        complex(ep), intent(inout) :: A(lda,*)
        call zsyr_quad(uplo, n, alpha, x, incx, A, lda)
    end subroutine zsyr

    subroutine zspmv(uplo, n, alpha, AP, x, incx, beta, y, incy)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, incx, incy
        complex(ep), intent(in)    :: alpha, beta, AP(*), x(*)
        complex(ep), intent(inout) :: y(*)
        call zspmv_quad(uplo, n, alpha, AP, x, incx, beta, y, incy)
    end subroutine zspmv

    subroutine zspr(uplo, n, alpha, x, incx, AP)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, incx
        complex(ep), intent(in)    :: alpha, x(*)
        complex(ep), intent(inout) :: AP(*)
        call zspr_quad(uplo, n, alpha, x, incx, AP)
    end subroutine zspr

    subroutine dsyswapr(uplo, n, A, lda, i1, i2)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, i1, i2
        real(ep),  intent(inout) :: A(lda,*)
        call dsyswapr_quad(uplo, n, A, lda, i1, i2)
    end subroutine dsyswapr

    subroutine zheswapr(uplo, n, A, lda, i1, i2)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, i1, i2
        complex(ep), intent(inout) :: A(lda,*)
        call zheswapr_quad(uplo, n, A, lda, i1, i2)
    end subroutine zheswapr

    subroutine zsyswapr(uplo, n, A, lda, i1, i2)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, i1, i2
        complex(ep), intent(inout) :: A(lda,*)
        call zsyswapr_quad(uplo, n, A, lda, i1, i2)
    end subroutine zsyswapr

    subroutine dsyconv(uplo, way, n, A, lda, ipiv, e, info)
        character, intent(in)    :: uplo, way
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(out)   :: e(*)
        integer,   intent(out)   :: info
        call dsyconv_quad(uplo, way, n, A, lda, ipiv, e, info)
    end subroutine dsyconv

    subroutine dsyconvf(uplo, way, n, A, lda, e, ipiv, info)
        character, intent(in)    :: uplo, way
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*), e(*)
        integer,   intent(inout) :: ipiv(*)
        integer,   intent(out)   :: info
        call dsyconvf_quad(uplo, way, n, A, lda, e, ipiv, info)
    end subroutine dsyconvf

    subroutine dsyconvf_rook(uplo, way, n, A, lda, e, ipiv, info)
        character, intent(in)    :: uplo, way
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*), e(*)
        integer,   intent(in)    :: ipiv(*)
        integer,   intent(out)   :: info
        call dsyconvf_rook_quad(uplo, way, n, A, lda, e, ipiv, info)
    end subroutine dsyconvf_rook

    subroutine zsyconv(uplo, way, n, A, lda, ipiv, e, info)
        character,   intent(in)    :: uplo, way
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(out)   :: e(*)
        integer,     intent(out)   :: info
        call zsyconv_quad(uplo, way, n, A, lda, ipiv, e, info)
    end subroutine zsyconv

    subroutine zsyconvf(uplo, way, n, A, lda, e, ipiv, info)
        character,   intent(in)    :: uplo, way
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*), e(*)
        integer,     intent(inout) :: ipiv(*)
        integer,     intent(out)   :: info
        call zsyconvf_quad(uplo, way, n, A, lda, e, ipiv, info)
    end subroutine zsyconvf

    subroutine zsyconvf_rook(uplo, way, n, A, lda, e, ipiv, info)
        character,   intent(in)    :: uplo, way
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*), e(*)
        integer,     intent(in)    :: ipiv(*)
        integer,     intent(out)   :: info
        call zsyconvf_rook_quad(uplo, way, n, A, lda, e, ipiv, info)
    end subroutine zsyconvf_rook

    subroutine zhprfs(uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, ferr, berr, work, rwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, nrhs, ldb, ldx
        complex(ep), intent(in)    :: AP(*), AFP(*), B(ldb,*)
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zhprfs_quad(uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, ferr, berr, work, rwork, info)
    end subroutine zhprfs

    subroutine zhpsvx(fact, uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, rcond, ferr, berr, work, &
                      rwork, info)
        character,   intent(in)    :: fact, uplo
        integer,     intent(in)    :: n, nrhs, ldb, ldx
        complex(ep), intent(in)    :: AP(*), B(ldb,*)
        complex(ep), intent(inout) :: AFP(*)
        integer,     intent(inout) :: ipiv(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zhpsvx_quad(fact, uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, rcond, ferr, berr, work, &
                         rwork, info)
    end subroutine zhpsvx

    subroutine zspsvx(fact, uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, rcond, ferr, berr, work, &
                      rwork, info)
        character,   intent(in)    :: fact, uplo
        integer,     intent(in)    :: n, nrhs, ldb, ldx
        complex(ep), intent(in)    :: AP(*), B(ldb,*)
        complex(ep), intent(inout) :: AFP(*)
        integer,     intent(inout) :: ipiv(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zspsvx_quad(fact, uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, rcond, ferr, berr, work, &
                         rwork, info)
    end subroutine zspsvx

    subroutine dspsvx(fact, uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, rcond, ferr, berr, work, &
                      iwork, info)
        character, intent(in)    :: fact, uplo
        integer,   intent(in)    :: n, nrhs, ldb, ldx
        real(ep),  intent(in)    :: AP(*), B(ldb,*)
        real(ep),  intent(inout) :: AFP(*)
        integer,   intent(inout) :: ipiv(*)
        real(ep),  intent(inout) :: X(ldx,*)
        real(ep),  intent(out)   :: rcond, ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dspsvx_quad(fact, uplo, n, nrhs, AP, AFP, ipiv, B, ldb, X, ldx, rcond, ferr, berr, work, &
                         iwork, info)
    end subroutine dspsvx

    subroutine dgtsv(n, nrhs, dl, d, du, B, ldb, info)
        integer,  intent(in)    :: n, nrhs, ldb
        real(ep), intent(inout) :: dl(*), d(*), du(*), B(ldb,*)
        integer,  intent(out)   :: info
        call dgtsv_quad(n, nrhs, dl, d, du, B, ldb, info)
    end subroutine dgtsv

    subroutine dptsv(n, nrhs, d, e, B, ldb, info)
        integer,  intent(in)    :: n, nrhs, ldb
        real(ep), intent(inout) :: d(*), e(*), B(ldb,*)
        integer,  intent(out)   :: info
        call dptsv_quad(n, nrhs, d, e, B, ldb, info)
    end subroutine dptsv

    subroutine dgtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, B, ldb, X, ldx, &
                      rcond, ferr, berr, work, iwork, info)
        character, intent(in)    :: fact, trans
        integer,   intent(in)    :: n, nrhs, ldb, ldx
        real(ep),  intent(in)    :: dl(*), d(*), du(*), B(ldb,*)
        real(ep),  intent(inout) :: dlf(*), df(*), duf(*), du2(*)
        integer,   intent(inout) :: ipiv(*)
        real(ep),  intent(inout) :: X(ldx,*)
        real(ep),  intent(out)   :: rcond, ferr(*), berr(*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgtsvx_quad(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, B, ldb, X, ldx, &
                         rcond, ferr, berr, work, iwork, info)
    end subroutine dgtsvx

    subroutine zgtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, B, ldb, X, ldx, &
                      rcond, ferr, berr, work, rwork, info)
        character,   intent(in)    :: fact, trans
        integer,     intent(in)    :: n, nrhs, ldb, ldx
        complex(ep), intent(in)    :: dl(*), d(*), du(*), B(ldb,*)
        complex(ep), intent(inout) :: dlf(*), df(*), duf(*), du2(*)
        integer,     intent(inout) :: ipiv(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zgtsvx_quad(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, B, ldb, X, ldx, &
                         rcond, ferr, berr, work, rwork, info)
    end subroutine zgtsvx

    subroutine dptsvx(fact, n, nrhs, d, e, df, ef, B, ldb, X, ldx, rcond, ferr, berr, work, info)
        character, intent(in)    :: fact
        integer,   intent(in)    :: n, nrhs, ldb, ldx
        real(ep),  intent(in)    :: d(*), e(*), B(ldb,*)
        real(ep),  intent(inout) :: df(*), ef(*)
        real(ep),  intent(inout) :: X(ldx,*)
        real(ep),  intent(out)   :: rcond, ferr(*), berr(*), work(*)
        integer,   intent(out)   :: info
        call dptsvx_quad(fact, n, nrhs, d, e, df, ef, B, ldb, X, ldx, rcond, ferr, berr, work, info)
    end subroutine dptsvx

    subroutine zptsvx(fact, n, nrhs, d, e, df, ef, B, ldb, X, ldx, rcond, ferr, berr, work, rwork, &
                      info)
        character,   intent(in)    :: fact
        integer,     intent(in)    :: n, nrhs, ldb, ldx
        real(ep),    intent(in)    :: d(*)
        complex(ep), intent(in)    :: e(*), B(ldb,*)
        real(ep),    intent(inout) :: df(*)
        complex(ep), intent(inout) :: ef(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(out)   :: rcond, ferr(*), berr(*), rwork(*)
        complex(ep), intent(out)   :: work(*)
        integer,     intent(out)   :: info
        call zptsvx_quad(fact, n, nrhs, d, e, df, ef, B, ldb, X, ldx, rcond, ferr, berr, work, rwork, &
                         info)
    end subroutine zptsvx

    subroutine dtrsna(job, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, S, SEP, mm, m, work, ldwork, &
                      iwork, info)
        character, intent(in)    :: job, howmny
        logical,   intent(in)    :: sel(*)
        integer,   intent(in)    :: n, ldt, ldvl, ldvr, mm, ldwork
        real(ep),  intent(in)    :: T(ldt,*), VL(ldvl,*), VR(ldvr,*)
        real(ep),  intent(out)   :: S(*), SEP(*), work(ldwork,*)
        integer,   intent(out)   :: m, iwork(*), info
        call dtrsna_quad(job, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, S, SEP, mm, m, work, &
                         ldwork, iwork, info)
    end subroutine dtrsna

    subroutine ztrsna(job, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, S, SEP, mm, m, work, ldwork, &
                      rwork, info)
        character,   intent(in)    :: job, howmny
        logical,     intent(in)    :: sel(*)
        integer,     intent(in)    :: n, ldt, ldvl, ldvr, mm, ldwork
        complex(ep), intent(in)    :: T(ldt,*), VL(ldvl,*), VR(ldvr,*)
        real(ep),    intent(out)   :: S(*), SEP(*), rwork(*)
        complex(ep), intent(out)   :: work(ldwork,*)
        integer,     intent(out)   :: m, info
        call ztrsna_quad(job, howmny, sel, n, T, ldt, VL, ldvl, VR, ldvr, S, SEP, mm, m, work, &
                         ldwork, rwork, info)
    end subroutine ztrsna

    subroutine dsytrd_2stage(vect, uplo, n, A, lda, D, E, tau, hous2, lhous2, work, lwork, info)
        character, intent(in)    :: vect, uplo
        integer,   intent(in)    :: n, lda, lhous2, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: D(*), E(*), tau(*), hous2(*), work(*)
        integer,   intent(out)   :: info
        call dsytrd_2stage_quad(vect, uplo, n, A, lda, D, E, tau, hous2, lhous2, work, lwork, info)
    end subroutine dsytrd_2stage

    subroutine zhetrd_2stage(vect, uplo, n, A, lda, D, E, tau, hous2, lhous2, work, lwork, info)
        character,   intent(in)    :: vect, uplo
        integer,     intent(in)    :: n, lda, lhous2, lwork
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: D(*), E(*)
        complex(ep), intent(out)   :: tau(*), hous2(*), work(*)
        integer,     intent(out)   :: info
        call zhetrd_2stage_quad(vect, uplo, n, A, lda, D, E, tau, hous2, lhous2, work, lwork, info)
    end subroutine zhetrd_2stage

    subroutine dsytrd_sy2sb(uplo, n, kd, A, lda, AB, ldab, tau, work, lwork, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, kd, lda, ldab, lwork
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: AB(ldab,*), tau(*), work(*)
        integer,   intent(out)   :: info
        call dsytrd_sy2sb_quad(uplo, n, kd, A, lda, AB, ldab, tau, work, lwork, info)
    end subroutine dsytrd_sy2sb

    subroutine zhetrd_he2hb(uplo, n, kd, A, lda, AB, ldab, tau, work, lwork, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, kd, lda, ldab, lwork
        complex(ep), intent(inout) :: A(lda,*)
        complex(ep), intent(out)   :: AB(ldab,*), tau(*), work(*)
        integer,     intent(out)   :: info
        call zhetrd_he2hb_quad(uplo, n, kd, A, lda, AB, ldab, tau, work, lwork, info)
    end subroutine zhetrd_he2hb

    subroutine dgedmd(jobs, jobz, jobr, jobf, whtsvd, m, n, X, ldx, Y, ldy, nrnk, tol, k, reig, &
                      imeig, Z, ldz, res, B, ldb, W, ldw, S, lds, work, lwork, iwork, liwork, info)
        character, intent(in)    :: jobs, jobz, jobr, jobf
        integer,   intent(in)    :: whtsvd, m, n, ldx, ldy, nrnk, ldz, ldb, &
        ldw, lds, lwork, liwork
        real(ep),  intent(in)    :: tol
        real(ep),  intent(inout) :: X(ldx,*), Y(ldy,*)
        integer,   intent(out)   :: k, info
        real(ep),  intent(out)   :: reig(*), imeig(*), Z(ldz,*), res(*),    &
        B(ldb,*), W(ldw,*), S(lds,*), work(*)
        integer,   intent(out)   :: iwork(*)
        call dgedmd_quad(jobs, jobz, jobr, jobf, whtsvd, m, n, X, ldx, Y, ldy, nrnk, tol, k, reig, &
                         imeig, Z, ldz, res, B, ldb, W, ldw, S, lds, work, lwork, iwork, liwork, &
                         info)
    end subroutine dgedmd

    subroutine zgedmd(jobs, jobz, jobr, jobf, whtsvd, m, n, X, ldx, Y, ldy, nrnk, tol, k, eigs, Z, &
                      ldz, res, B, ldb, W, ldw, S, lds, zwork, lzwork, rwork, lrwork, iwork, liwork, &
                      info)
        character,   intent(in)    :: jobs, jobz, jobr, jobf
        integer,     intent(in)    :: whtsvd, m, n, ldx, ldy, nrnk, ldz, ldb, &
        ldw, lds, lzwork, lrwork, liwork
        real(ep),    intent(in)    :: tol
        complex(ep), intent(inout) :: X(ldx,*), Y(ldy,*)
        integer,     intent(out)   :: k, info
        complex(ep), intent(out)   :: eigs(*), Z(ldz,*), B(ldb,*), W(ldw,*), &
        S(lds,*), zwork(*)
        real(ep),    intent(out)   :: res(*), rwork(*)
        integer,     intent(out)   :: iwork(*)
        call zgedmd_quad(jobs, jobz, jobr, jobf, whtsvd, m, n, X, ldx, Y, ldy, nrnk, tol, k, eigs, Z, &
                         ldz, res, B, ldb, W, ldw, S, lds, zwork, lzwork, rwork, lrwork, iwork, &
                         liwork, info)
    end subroutine zgedmd

    subroutine dgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, F, ldf, X, ldx, Y, ldy, &
                       nrnk, tol, k, reig, imeig, Z, ldz, res, B, ldb, V, ldv, S, lds, work, lwork, &
                       iwork, liwork, info)
        character, intent(in)    :: jobs, jobz, jobr, jobq, jobt, jobf
        integer,   intent(in)    :: whtsvd, m, n, ldf, ldx, ldy, nrnk, ldz,  &
        ldb, ldv, lds, lwork, liwork
        real(ep),  intent(in)    :: tol
        real(ep),  intent(inout) :: F(ldf,*)
        integer,   intent(out)   :: k, info
        real(ep),  intent(out)   :: X(ldx,*), Y(ldy,*), Z(ldz,*), B(ldb,*),  &
        V(ldv,*), S(lds,*), reig(*), imeig(*),  &
        res(*), work(*)
        integer,   intent(out)   :: iwork(*)
        call dgedmdq_quad(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, F, ldf, X, ldx, Y, ldy, &
                          nrnk, tol, k, reig, imeig, Z, ldz, res, B, ldb, V, ldv, S, lds, work, &
                          lwork, iwork, liwork, info)
    end subroutine dgedmdq

    subroutine zgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, F, ldf, X, ldx, Y, ldy, &
                       nrnk, tol, k, eigs, Z, ldz, res, B, ldb, V, ldv, S, lds, zwork, lzwork, work, &
                       lwork, iwork, liwork, info)
        character,   intent(in)    :: jobs, jobz, jobr, jobq, jobt, jobf
        integer,     intent(in)    :: whtsvd, m, n, ldf, ldx, ldy, nrnk,     &
        ldz, ldb, ldv, lds, lzwork, lwork,    &
        liwork
        real(ep),    intent(in)    :: tol
        complex(ep), intent(inout) :: F(ldf,*)
        integer,     intent(out)   :: k, info
        complex(ep), intent(out)   :: X(ldx,*), Y(ldy,*), Z(ldz,*), B(ldb,*), &
        V(ldv,*), S(lds,*), eigs(*), zwork(*)
        real(ep),    intent(out)   :: res(*), work(*)
        integer,     intent(out)   :: iwork(*)
        call zgedmdq_quad(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, F, ldf, X, ldx, Y, ldy, &
                          nrnk, tol, k, eigs, Z, ldz, res, B, ldb, V, ldv, S, lds, zwork, lzwork, &
                          work, lwork, iwork, liwork, info)
    end subroutine zgedmdq

    subroutine dgesvxx(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, R, C, B, ldb, X, ldx, &
                       rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, &
                       params, work, iwork, info)
        character, intent(in)    :: fact, trans
        character, intent(inout) :: equed
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        real(ep),  intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*), R(*), C(*), params(*)
        integer,   intent(inout) :: ipiv(*)
        real(ep),  intent(out)   :: X(ldx,*), rcond, rpvgrw, berr(*),       &
        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgesvxx_quad(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, R, C, B, ldb, X, ldx, &
                          rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, &
                          params, work, iwork, info)
    end subroutine dgesvxx

    subroutine zgesvxx(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, R, C, B, ldb, X, ldx, &
                       rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, &
                       params, work, rwork, info)
        character,   intent(in)    :: fact, trans
        character,   intent(inout) :: equed
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
        real(ep),    intent(inout) :: R(*), C(*), params(*)
        integer,     intent(inout) :: ipiv(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, rpvgrw, berr(*),               &
        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), rwork(*)
        integer,     intent(out)   :: info
        call zgesvxx_quad(fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, R, C, B, ldb, X, ldx, &
                          rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, &
                          params, work, rwork, info)
    end subroutine zgesvxx

    subroutine dgerfsx(trans, equed, n, nrhs, A, lda, AF, ldaf, ipiv, R, C, B, ldb, X, ldx, rcond, &
                       berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, &
                       info)
        character, intent(in)    :: trans, equed
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*), R(*), C(*)
        real(ep),  intent(inout) :: X(ldx,*), params(*)
        real(ep),  intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*),  &
        err_bnds_comp(nrhs,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgerfsx_quad(trans, equed, n, nrhs, A, lda, AF, ldaf, ipiv, R, C, B, ldb, X, ldx, rcond, &
                          berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, &
                          iwork, info)
    end subroutine dgerfsx

    subroutine zgerfsx(trans, equed, n, nrhs, A, lda, AF, ldaf, ipiv, R, C, B, ldb, X, ldx, rcond, &
                       berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, &
                       info)
        character,   intent(in)    :: trans, equed
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
        real(ep),    intent(in)    :: R(*), C(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(inout) :: params(*)
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
        err_bnds_comp(nrhs,*), rwork(*)
        integer,     intent(out)   :: info
        call zgerfsx_quad(trans, equed, n, nrhs, A, lda, AF, ldaf, ipiv, R, C, B, ldb, X, ldx, rcond, &
                          berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, &
                          rwork, info)
    end subroutine zgerfsx

    subroutine dgbsvxx(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, equed, R, C, B, ldb, &
                       X, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, &
                       nparams, params, work, iwork, info)
        character, intent(in)    :: fact, trans
        character, intent(inout) :: equed
        integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx, n_err_bnds, nparams
        real(ep),  intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*), R(*), C(*), params(*)
        integer,   intent(inout) :: ipiv(*)
        real(ep),  intent(out)   :: X(ldx,*), rcond, rpvgrw, berr(*),       &
        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgbsvxx_quad(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, equed, R, C, B, &
                          ldb, X, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, &
                          nparams, params, work, iwork, info)
    end subroutine dgbsvxx

    subroutine zgbsvxx(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, equed, R, C, B, ldb, &
                       X, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, &
                       nparams, params, work, rwork, info)
        character,   intent(in)    :: fact, trans
        character,   intent(inout) :: equed
        integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx, n_err_bnds, nparams
        complex(ep), intent(inout) :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
        real(ep),    intent(inout) :: R(*), C(*), params(*)
        integer,     intent(inout) :: ipiv(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, rpvgrw, berr(*),               &
        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), rwork(*)
        integer,     intent(out)   :: info
        call zgbsvxx_quad(fact, trans, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, equed, R, C, B, &
                          ldb, X, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, &
                          nparams, params, work, rwork, info)
    end subroutine zgbsvxx

    subroutine dgbrfsx(trans, equed, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, R, C, B, ldb, X, &
                       ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, &
                       work, iwork, info)
        character, intent(in)    :: trans, equed
        integer,   intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx, n_err_bnds, nparams
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*), R(*), C(*)
        real(ep),  intent(inout) :: X(ldx,*), params(*)
        real(ep),  intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
        err_bnds_comp(nrhs,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dgbrfsx_quad(trans, equed, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, R, C, B, ldb, X, &
                          ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, &
                          params, work, iwork, info)
    end subroutine dgbrfsx

    subroutine zgbrfsx(trans, equed, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, R, C, B, ldb, X, &
                       ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, &
                       work, rwork, info)
        character,   intent(in)    :: trans, equed
        integer,     intent(in)    :: n, kl, ku, nrhs, ldab, ldafb, ldb, ldx, n_err_bnds, nparams
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(in)    :: AB(ldab,*), AFB(ldafb,*), B(ldb,*)
        real(ep),    intent(in)    :: R(*), C(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(inout) :: params(*)
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
        err_bnds_comp(nrhs,*), rwork(*)
        integer,     intent(out)   :: info
        call zgbrfsx_quad(trans, equed, n, kl, ku, nrhs, AB, ldab, AFB, ldafb, ipiv, R, C, B, ldb, X, &
                          ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, &
                          params, work, rwork, info)
    end subroutine zgbrfsx

    subroutine dposvxx(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S, B, ldb, X, ldx, rcond, &
                       rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, &
                       iwork, info)
        character, intent(in)    :: fact, uplo
        character, intent(inout) :: equed
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        real(ep),  intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*), S(*), params(*)
        real(ep),  intent(out)   :: X(ldx,*), rcond, rpvgrw, berr(*),       &
        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dposvxx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S, B, ldb, X, ldx, rcond, &
                          rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, &
                          work, iwork, info)
    end subroutine dposvxx

    subroutine zposvxx(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S, B, ldb, X, ldx, rcond, &
                       rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, &
                       rwork, info)
        character,   intent(in)    :: fact, uplo
        character,   intent(inout) :: equed
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
        real(ep),    intent(inout) :: S(*), params(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, rpvgrw, berr(*),               &
        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), rwork(*)
        integer,     intent(out)   :: info
        call zposvxx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, equed, S, B, ldb, X, ldx, rcond, &
                          rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, &
                          work, rwork, info)
    end subroutine zposvxx

    subroutine dporfsx(uplo, equed, n, nrhs, A, lda, AF, ldaf, S, B, ldb, X, ldx, rcond, berr, &
                       n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info)
        character, intent(in)    :: uplo, equed
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), S(*), B(ldb,*)
        real(ep),  intent(inout) :: X(ldx,*), params(*)
        real(ep),  intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*),  &
        err_bnds_comp(nrhs,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dporfsx_quad(uplo, equed, n, nrhs, A, lda, AF, ldaf, S, B, ldb, X, ldx, rcond, berr, &
                          n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, &
                          info)
    end subroutine dporfsx

    subroutine zporfsx(uplo, equed, n, nrhs, A, lda, AF, ldaf, S, B, ldb, X, ldx, rcond, berr, &
                       n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info)
        character,   intent(in)    :: uplo, equed
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
        real(ep),    intent(in)    :: S(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(inout) :: params(*)
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
        err_bnds_comp(nrhs,*), rwork(*)
        integer,     intent(out)   :: info
        call zporfsx_quad(uplo, equed, n, nrhs, A, lda, AF, ldaf, S, B, ldb, X, ldx, rcond, berr, &
                          n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, &
                          info)
    end subroutine zporfsx

    subroutine dsysvxx(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, equed, S, B, ldb, X, ldx, rcond, &
                       rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, &
                       iwork, info)
        character, intent(in)    :: fact, uplo
        character, intent(inout) :: equed
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        real(ep),  intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*), S(*), params(*)
        integer,   intent(inout) :: ipiv(*)
        real(ep),  intent(out)   :: X(ldx,*), rcond, rpvgrw, berr(*),       &
        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsysvxx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, equed, S, B, ldb, X, ldx, &
                          rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, &
                          params, work, iwork, info)
    end subroutine dsysvxx

    subroutine zsysvxx(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, equed, S, B, ldb, X, ldx, rcond, &
                       rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, &
                       rwork, info)
        character,   intent(in)    :: fact, uplo
        character,   intent(inout) :: equed
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
        real(ep),    intent(inout) :: S(*), params(*)
        integer,     intent(inout) :: ipiv(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, rpvgrw, berr(*),               &
        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), rwork(*)
        integer,     intent(out)   :: info
        call zsysvxx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, equed, S, B, ldb, X, ldx, &
                          rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, &
                          params, work, rwork, info)
    end subroutine zsysvxx

    subroutine dsyrfsx(uplo, equed, n, nrhs, A, lda, AF, ldaf, ipiv, S, B, ldb, X, ldx, rcond, berr, &
                       n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork, info)
        character, intent(in)    :: uplo, equed
        integer,   intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        integer,   intent(in)    :: ipiv(*)
        real(ep),  intent(in)    :: A(lda,*), AF(ldaf,*), S(*), B(ldb,*)
        real(ep),  intent(inout) :: X(ldx,*), params(*)
        real(ep),  intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
        err_bnds_comp(nrhs,*), work(*)
        integer,   intent(out)   :: iwork(*), info
        call dsyrfsx_quad(uplo, equed, n, nrhs, A, lda, AF, ldaf, ipiv, S, B, ldb, X, ldx, rcond, &
                          berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, &
                          iwork, info)
    end subroutine dsyrfsx

    subroutine zsyrfsx(uplo, equed, n, nrhs, A, lda, AF, ldaf, ipiv, S, B, ldb, X, ldx, rcond, berr, &
                       n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info)
        character,   intent(in)    :: uplo, equed
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
        real(ep),    intent(in)    :: S(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(inout) :: params(*)
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
        err_bnds_comp(nrhs,*), rwork(*)
        integer,     intent(out)   :: info
        call zsyrfsx_quad(uplo, equed, n, nrhs, A, lda, AF, ldaf, ipiv, S, B, ldb, X, ldx, rcond, &
                          berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, &
                          rwork, info)
    end subroutine zsyrfsx

    subroutine zhesvxx(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, equed, S, B, ldb, X, ldx, rcond, &
                       rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, &
                       rwork, info)
        character,   intent(in)    :: fact, uplo
        character,   intent(inout) :: equed
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        complex(ep), intent(inout) :: A(lda,*), AF(ldaf,*), B(ldb,*)
        real(ep),    intent(inout) :: S(*), params(*)
        integer,     intent(inout) :: ipiv(*)
        complex(ep), intent(out)   :: X(ldx,*), work(*)
        real(ep),    intent(out)   :: rcond, rpvgrw, berr(*),               &
        err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*), rwork(*)
        integer,     intent(out)   :: info
        call zhesvxx_quad(fact, uplo, n, nrhs, A, lda, AF, ldaf, ipiv, equed, S, B, ldb, X, ldx, &
                          rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, &
                          params, work, rwork, info)
    end subroutine zhesvxx

    subroutine zherfsx(uplo, equed, n, nrhs, A, lda, AF, ldaf, ipiv, S, B, ldb, X, ldx, rcond, berr, &
                       n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info)
        character,   intent(in)    :: uplo, equed
        integer,     intent(in)    :: n, nrhs, lda, ldaf, ldb, ldx, n_err_bnds, nparams
        integer,     intent(in)    :: ipiv(*)
        complex(ep), intent(in)    :: A(lda,*), AF(ldaf,*), B(ldb,*)
        real(ep),    intent(in)    :: S(*)
        complex(ep), intent(inout) :: X(ldx,*)
        real(ep),    intent(inout) :: params(*)
        complex(ep), intent(out)   :: work(*)
        real(ep),    intent(out)   :: rcond, berr(*), err_bnds_norm(nrhs,*), &
        err_bnds_comp(nrhs,*), rwork(*)
        integer,     intent(out)   :: info
        call zherfsx_quad(uplo, equed, n, nrhs, A, lda, AF, ldaf, ipiv, S, B, ldb, X, ldx, rcond, &
                          berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, &
                          rwork, info)
    end subroutine zherfsx


end module ref_quad_lapack
