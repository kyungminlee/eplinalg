! Explicit interfaces to refblas_quad — the vendored Netlib BLAS
! compiled with gfortran's -freal-8-real-16 so REAL(KIND=8) and
! DOUBLE PRECISION entities are promoted to REAL(KIND=16) without
! routine renaming. The D/Z routines below operate at quad precision
! and serve as the universal reference for the differential precision
! tests.
!
! Calling these via implicit typing would be wrong: a caller without
! the explicit interface would default DDOT etc. to DOUBLE PRECISION
! (KIND=8), mismatching the promoted KIND=16 symbol. Always `use
! ref_quad_blas`.

! Module name is pbblas_ref_quad_blas (not ref_quad_blas) to avoid
! colliding with tests/blas/common/ref_quad_blas. Both targets write
! their .mod files into the same ${PROJECT_BINARY_DIR}/fmod directory,
! so the module names must differ — same convention as
! pblas_prec_report vs prec_report.
module pbblas_ref_quad_blas
    use prec_kinds, only: ep
    implicit none

    interface
        ! ── Level 1 — real ───────────────────────────────────────────
        function ddot(n, x, incx, y, incy) result(r)
            import :: ep
            integer,  intent(in) :: n, incx, incy
            real(ep), intent(in) :: x(*), y(*)
            real(ep) :: r
        end function ddot

        function dasum(n, x, incx) result(r)
            import :: ep
            integer,  intent(in) :: n, incx
            real(ep), intent(in) :: x(*)
            real(ep) :: r
        end function dasum

        function dnrm2(n, x, incx) result(r)
            import :: ep
            integer,  intent(in) :: n, incx
            real(ep), intent(in) :: x(*)
            real(ep) :: r
        end function dnrm2

        function idamax(n, x, incx) result(r)
            import :: ep
            integer,  intent(in) :: n, incx
            real(ep), intent(in) :: x(*)
            integer :: r
        end function idamax

        subroutine daxpy(n, alpha, x, incx, y, incy)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(in)    :: alpha, x(*)
            real(ep), intent(inout) :: y(*)
        end subroutine daxpy

        subroutine dcopy(n, x, incx, y, incy)
            import :: ep
            integer,  intent(in)  :: n, incx, incy
            real(ep), intent(in)  :: x(*)
            real(ep), intent(out) :: y(*)
        end subroutine dcopy

        subroutine dscal(n, alpha, x, incx)
            import :: ep
            integer,  intent(in)    :: n, incx
            real(ep), intent(in)    :: alpha
            real(ep), intent(inout) :: x(*)
        end subroutine dscal

        subroutine dswap(n, x, incx, y, incy)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(inout) :: x(*), y(*)
        end subroutine dswap

        subroutine drot(n, x, incx, y, incy, c, s)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(inout) :: x(*), y(*)
            real(ep), intent(in)    :: c, s
        end subroutine drot

        subroutine drotg(a, b, c, s)
            import :: ep
            real(ep), intent(inout) :: a, b
            real(ep), intent(out)   :: c, s
        end subroutine drotg

        subroutine drotm(n, x, incx, y, incy, param)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(inout) :: x(*), y(*)
            real(ep), intent(in)    :: param(5)
        end subroutine drotm

        subroutine drotmg(d1, d2, x1, y1, param)
            import :: ep
            real(ep), intent(inout) :: d1, d2, x1
            real(ep), intent(in)    :: y1
            real(ep), intent(out)   :: param(5)
        end subroutine drotmg

        ! ── Level 1 — complex ────────────────────────────────────────
        function zdotc(n, x, incx, y, incy) result(r)
            import :: ep
            integer,     intent(in) :: n, incx, incy
            complex(ep), intent(in) :: x(*), y(*)
            complex(ep) :: r
        end function zdotc

        function zdotu(n, x, incx, y, incy) result(r)
            import :: ep
            integer,     intent(in) :: n, incx, incy
            complex(ep), intent(in) :: x(*), y(*)
            complex(ep) :: r
        end function zdotu

        function dzasum(n, x, incx) result(r)
            import :: ep
            integer,     intent(in) :: n, incx
            complex(ep), intent(in) :: x(*)
            real(ep) :: r
        end function dzasum

        subroutine zaxpy(n, alpha, x, incx, y, incy)
            import :: ep
            integer,     intent(in)    :: n, incx, incy
            complex(ep), intent(in)    :: alpha, x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine zaxpy

        subroutine zscal(n, alpha, x, incx)
            import :: ep
            integer,     intent(in)    :: n, incx
            complex(ep), intent(in)    :: alpha
            complex(ep), intent(inout) :: x(*)
        end subroutine zscal

        subroutine zdscal(n, alpha, x, incx)
            import :: ep
            integer,     intent(in)    :: n, incx
            real(ep),    intent(in)    :: alpha
            complex(ep), intent(inout) :: x(*)
        end subroutine zdscal

        subroutine zcopy(n, x, incx, y, incy)
            import :: ep
            integer,     intent(in)  :: n, incx, incy
            complex(ep), intent(in)  :: x(*)
            complex(ep), intent(out) :: y(*)
        end subroutine zcopy

        subroutine zswap(n, x, incx, y, incy)
            import :: ep
            integer,     intent(in)    :: n, incx, incy
            complex(ep), intent(inout) :: x(*), y(*)
        end subroutine zswap

        function izamax(n, x, incx) result(r)
            import :: ep
            integer,     intent(in) :: n, incx
            complex(ep), intent(in) :: x(*)
            integer :: r
        end function izamax

        function dznrm2(n, x, incx) result(r)
            import :: ep
            integer,     intent(in) :: n, incx
            complex(ep), intent(in) :: x(*)
            real(ep) :: r
        end function dznrm2

        ! ── Level 2 — real ───────────────────────────────────────────
        subroutine dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine dgemv

        subroutine dgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, kl, ku, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine dgbmv

        subroutine dger(m, n, alpha, x, incx, y, incy, A, lda)
            import :: ep
            integer,  intent(in)    :: m, n, incx, incy, lda
            real(ep), intent(in)    :: alpha, x(*), y(*)
            real(ep), intent(inout) :: A(lda,*)
        end subroutine dger

        subroutine dsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine dsymv

        subroutine dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, incx, incy
            real(ep),  intent(in)    :: alpha, beta, ap(*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine dspmv

        subroutine dsbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, k, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine dsbmv

        subroutine dtrmv(uplo, trans, diag, n, A, lda, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, lda, incx
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine dtrmv

        subroutine dtbmv(uplo, trans, diag, n, k, A, lda, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, k, lda, incx
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine dtbmv

        subroutine dtpmv(uplo, trans, diag, n, ap, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, incx
            real(ep),  intent(in)    :: ap(*)
            real(ep),  intent(inout) :: x(*)
        end subroutine dtpmv

        subroutine dtrsv(uplo, trans, diag, n, A, lda, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, lda, incx
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine dtrsv

        subroutine dsyr(uplo, n, alpha, x, incx, A, lda)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, incx, lda
            real(ep),  intent(in)    :: alpha, x(*)
            real(ep),  intent(inout) :: A(lda,*)
        end subroutine dsyr

        subroutine dsyr2(uplo, n, alpha, x, incx, y, incy, A, lda)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, incx, incy, lda
            real(ep),  intent(in)    :: alpha, x(*), y(*)
            real(ep),  intent(inout) :: A(lda,*)
        end subroutine dsyr2

        ! ── Level 2 — complex ────────────────────────────────────────
        subroutine zgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, lda, incx, incy
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine zgemv

        subroutine zhemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, incx, incy
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine zhemv

        subroutine zgerc(m, n, alpha, x, incx, y, incy, A, lda)
            import :: ep
            integer,     intent(in)    :: m, n, incx, incy, lda
            complex(ep), intent(in)    :: alpha, x(*), y(*)
            complex(ep), intent(inout) :: A(lda,*)
        end subroutine zgerc

        subroutine zgeru(m, n, alpha, x, incx, y, incy, A, lda)
            import :: ep
            integer,     intent(in)    :: m, n, incx, incy, lda
            complex(ep), intent(in)    :: alpha, x(*), y(*)
            complex(ep), intent(inout) :: A(lda,*)
        end subroutine zgeru

        subroutine zher(uplo, n, alpha, x, incx, A, lda)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, incx, lda
            real(ep),    intent(in)    :: alpha
            complex(ep), intent(in)    :: x(*)
            complex(ep), intent(inout) :: A(lda,*)
        end subroutine zher

        subroutine zher2(uplo, n, alpha, x, incx, y, incy, A, lda)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, incx, incy, lda
            complex(ep), intent(in)    :: alpha, x(*), y(*)
            complex(ep), intent(inout) :: A(lda,*)
        end subroutine zher2

        subroutine ztrmv(uplo, trans, diag, n, A, lda, x, incx)
            import :: ep
            character,   intent(in)    :: uplo, trans, diag
            integer,     intent(in)    :: n, lda, incx
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: x(*)
        end subroutine ztrmv

        subroutine ztrsv(uplo, trans, diag, n, A, lda, x, incx)
            import :: ep
            character,   intent(in)    :: uplo, trans, diag
            integer,     intent(in)    :: n, lda, incx
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: x(*)
        end subroutine ztrsv

        ! ── Level 3 — real ───────────────────────────────────────────
        subroutine dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character, intent(in)    :: transa, transb
            integer,   intent(in)    :: m, n, k, lda, ldb, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine dgemm

        subroutine dsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character, intent(in)    :: side, uplo
            integer,   intent(in)    :: m, n, lda, ldb, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine dsymm

        subroutine dsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
            import :: ep
            character, intent(in)    :: uplo, trans
            integer,   intent(in)    :: n, k, lda, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine dsyrk

        subroutine dsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character, intent(in)    :: uplo, trans
            integer,   intent(in)    :: n, k, lda, ldb, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine dsyr2k

        subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: ep
            character, intent(in)    :: side, uplo, transa, diag
            integer,   intent(in)    :: m, n, lda, ldb
            real(ep),  intent(in)    :: alpha
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
        end subroutine dtrmm

        subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: ep
            character, intent(in)    :: side, uplo, transa, diag
            integer,   intent(in)    :: m, n, lda, ldb
            real(ep),  intent(in)    :: alpha
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
        end subroutine dtrsm

        ! ── Level 3 — complex ────────────────────────────────────────
        subroutine zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: transa, transb
            integer,     intent(in)    :: m, n, k, lda, ldb, ldc
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine zgemm

        subroutine zhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: side, uplo
            integer,     intent(in)    :: m, n, lda, ldb, ldc
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine zhemm

        subroutine zherk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: uplo, trans
            integer,     intent(in)    :: n, k, lda, ldc
            real(ep),    intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine zherk

        subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: ep
            character,   intent(in)    :: side, uplo, transa, diag
            integer,     intent(in)    :: m, n, lda, ldb
            complex(ep), intent(in)    :: alpha
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: B(ldb,*)
        end subroutine ztrsm

        subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: ep
            character,   intent(in)    :: side, uplo, transa, diag
            integer,     intent(in)    :: m, n, lda, ldb
            complex(ep), intent(in)    :: alpha
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: B(ldb,*)
        end subroutine ztrmm

        subroutine zsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: side, uplo
            integer,     intent(in)    :: m, n, lda, ldb, ldc
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine zsymm

        subroutine zsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: uplo, trans
            integer,     intent(in)    :: n, k, lda, ldc
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine zsyrk

        subroutine zsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: uplo, trans
            integer,     intent(in)    :: n, k, lda, ldb, ldc
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine zsyr2k

        subroutine zher2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: uplo, trans
            integer,     intent(in)    :: n, k, lda, ldb, ldc
            complex(ep), intent(in)    :: alpha
            real(ep),    intent(in)    :: beta
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine zher2k
    end interface

end module pbblas_ref_quad_blas
