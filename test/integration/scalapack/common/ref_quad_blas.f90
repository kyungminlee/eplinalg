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

module ref_quad_blas
    use prec_kinds, only: ep
    implicit none

    interface
        ! ── Level 1 — real ───────────────────────────────────────────
        function ddot_quad(n, x, incx, y, incy) result(r)
            import :: ep
            integer,  intent(in) :: n, incx, incy
            real(ep), intent(in) :: x(*), y(*)
            real(ep) :: r
        end function ddot_quad

        function dasum_quad(n, x, incx) result(r)
            import :: ep
            integer,  intent(in) :: n, incx
            real(ep), intent(in) :: x(*)
            real(ep) :: r
        end function dasum_quad

        function dnrm2_quad(n, x, incx) result(r)
            import :: ep
            integer,  intent(in) :: n, incx
            real(ep), intent(in) :: x(*)
            real(ep) :: r
        end function dnrm2_quad

        function idamax_quad(n, x, incx) result(r)
            import :: ep
            integer,  intent(in) :: n, incx
            real(ep), intent(in) :: x(*)
            integer :: r
        end function idamax_quad

        subroutine daxpy_quad(n, alpha, x, incx, y, incy)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(in)    :: alpha, x(*)
            real(ep), intent(inout) :: y(*)
        end subroutine daxpy_quad

        subroutine dcopy_quad(n, x, incx, y, incy)
            import :: ep
            integer,  intent(in)  :: n, incx, incy
            real(ep), intent(in)  :: x(*)
            real(ep), intent(out) :: y(*)
        end subroutine dcopy_quad

        subroutine dscal_quad(n, alpha, x, incx)
            import :: ep
            integer,  intent(in)    :: n, incx
            real(ep), intent(in)    :: alpha
            real(ep), intent(inout) :: x(*)
        end subroutine dscal_quad

        subroutine dswap_quad(n, x, incx, y, incy)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(inout) :: x(*), y(*)
        end subroutine dswap_quad

        subroutine drot_quad(n, x, incx, y, incy, c, s)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(inout) :: x(*), y(*)
            real(ep), intent(in)    :: c, s
        end subroutine drot_quad

        subroutine drotg_quad(a, b, c, s)
            import :: ep
            real(ep), intent(inout) :: a, b
            real(ep), intent(out)   :: c, s
        end subroutine drotg_quad

        subroutine drotm_quad(n, x, incx, y, incy, param)
            import :: ep
            integer,  intent(in)    :: n, incx, incy
            real(ep), intent(inout) :: x(*), y(*)
            real(ep), intent(in)    :: param(5)
        end subroutine drotm_quad

        subroutine drotmg_quad(d1, d2, x1, y1, param)
            import :: ep
            real(ep), intent(inout) :: d1, d2, x1
            real(ep), intent(in)    :: y1
            real(ep), intent(out)   :: param(5)
        end subroutine drotmg_quad

        ! ── Level 1 — complex ────────────────────────────────────────
        function zdotc_quad(n, x, incx, y, incy) result(r)
            import :: ep
            integer,     intent(in) :: n, incx, incy
            complex(ep), intent(in) :: x(*), y(*)
            complex(ep) :: r
        end function zdotc_quad

        function zdotu_quad(n, x, incx, y, incy) result(r)
            import :: ep
            integer,     intent(in) :: n, incx, incy
            complex(ep), intent(in) :: x(*), y(*)
            complex(ep) :: r
        end function zdotu_quad

        function dzasum_quad(n, x, incx) result(r)
            import :: ep
            integer,     intent(in) :: n, incx
            complex(ep), intent(in) :: x(*)
            real(ep) :: r
        end function dzasum_quad

        subroutine zaxpy_quad(n, alpha, x, incx, y, incy)
            import :: ep
            integer,     intent(in)    :: n, incx, incy
            complex(ep), intent(in)    :: alpha, x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine zaxpy_quad

        subroutine zscal_quad(n, alpha, x, incx)
            import :: ep
            integer,     intent(in)    :: n, incx
            complex(ep), intent(in)    :: alpha
            complex(ep), intent(inout) :: x(*)
        end subroutine zscal_quad

        ! ── Level 2 — real ───────────────────────────────────────────
        subroutine dgemv_quad(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine dgemv_quad

        subroutine dgbmv_quad(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: trans
            integer,   intent(in)    :: m, n, kl, ku, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine dgbmv_quad

        subroutine dger_quad(m, n, alpha, x, incx, y, incy, A, lda)
            import :: ep
            integer,  intent(in)    :: m, n, incx, incy, lda
            real(ep), intent(in)    :: alpha, x(*), y(*)
            real(ep), intent(inout) :: A(lda,*)
        end subroutine dger_quad

        subroutine dsymv_quad(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine dsymv_quad

        subroutine dspmv_quad(uplo, n, alpha, ap, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, incx, incy
            real(ep),  intent(in)    :: alpha, beta, ap(*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine dspmv_quad

        subroutine dsbmv_quad(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in)    :: uplo
            integer,   intent(in)    :: n, k, lda, incx, incy
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine dsbmv_quad

        subroutine dtrmv_quad(uplo, trans, diag, n, A, lda, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, lda, incx
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine dtrmv_quad

        subroutine dtbmv_quad(uplo, trans, diag, n, k, A, lda, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, k, lda, incx
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine dtbmv_quad

        subroutine dtpmv_quad(uplo, trans, diag, n, ap, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, incx
            real(ep),  intent(in)    :: ap(*)
            real(ep),  intent(inout) :: x(*)
        end subroutine dtpmv_quad

        subroutine dtrsv_quad(uplo, trans, diag, n, A, lda, x, incx)
            import :: ep
            character, intent(in)    :: uplo, trans, diag
            integer,   intent(in)    :: n, lda, incx
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine dtrsv_quad

        ! ── Level 2 — complex ────────────────────────────────────────
        subroutine zgemv_quad(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character,   intent(in)    :: trans
            integer,     intent(in)    :: m, n, lda, incx, incy
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine zgemv_quad

        subroutine zhemv_quad(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
            import :: ep
            character,   intent(in)    :: uplo
            integer,     intent(in)    :: n, lda, incx, incy
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), x(*)
            complex(ep), intent(inout) :: y(*)
        end subroutine zhemv_quad

        subroutine zgerc_quad(m, n, alpha, x, incx, y, incy, A, lda)
            import :: ep
            integer,     intent(in)    :: m, n, incx, incy, lda
            complex(ep), intent(in)    :: alpha, x(*), y(*)
            complex(ep), intent(inout) :: A(lda,*)
        end subroutine zgerc_quad

        ! ── Level 3 — real ───────────────────────────────────────────
        subroutine dgemm_quad(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character, intent(in)    :: transa, transb
            integer,   intent(in)    :: m, n, k, lda, ldb, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine dgemm_quad

        subroutine dsymm_quad(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character, intent(in)    :: side, uplo
            integer,   intent(in)    :: m, n, lda, ldb, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine dsymm_quad

        subroutine dsyrk_quad(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
            import :: ep
            character, intent(in)    :: uplo, trans
            integer,   intent(in)    :: n, k, lda, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine dsyrk_quad

        subroutine dsyr2k_quad(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character, intent(in)    :: uplo, trans
            integer,   intent(in)    :: n, k, lda, ldb, ldc
            real(ep),  intent(in)    :: alpha, beta
            real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
            real(ep),  intent(inout) :: C(ldc,*)
        end subroutine dsyr2k_quad

        subroutine dtrmm_quad(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: ep
            character, intent(in)    :: side, uplo, transa, diag
            integer,   intent(in)    :: m, n, lda, ldb
            real(ep),  intent(in)    :: alpha
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
        end subroutine dtrmm_quad

        subroutine dtrsm_quad(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: ep
            character, intent(in)    :: side, uplo, transa, diag
            integer,   intent(in)    :: m, n, lda, ldb
            real(ep),  intent(in)    :: alpha
            real(ep),  intent(in)    :: A(lda,*)
            real(ep),  intent(inout) :: B(ldb,*)
        end subroutine dtrsm_quad

        ! ── Level 3 — complex ────────────────────────────────────────
        subroutine zgemm_quad(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: transa, transb
            integer,     intent(in)    :: m, n, k, lda, ldb, ldc
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine zgemm_quad

        subroutine zhemm_quad(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: side, uplo
            integer,     intent(in)    :: m, n, lda, ldb, ldc
            complex(ep), intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine zhemm_quad

        subroutine zherk_quad(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
            import :: ep
            character,   intent(in)    :: uplo, trans
            integer,     intent(in)    :: n, k, lda, ldc
            real(ep),    intent(in)    :: alpha, beta
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: C(ldc,*)
        end subroutine zherk_quad

        subroutine ztrsm_quad(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
            import :: ep
            character,   intent(in)    :: side, uplo, transa, diag
            integer,     intent(in)    :: m, n, lda, ldb
            complex(ep), intent(in)    :: alpha
            complex(ep), intent(in)    :: A(lda,*)
            complex(ep), intent(inout) :: B(ldb,*)
        end subroutine ztrsm_quad
    end interface

contains

    function ddot(n, x, incx, y, incy) result(r)
        integer,  intent(in) :: n, incx, incy
        real(ep), intent(in) :: x(*), y(*)
        real(ep) :: r
        r = ddot_quad(n, x, incx, y, incy)
    end function ddot

    function dasum(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        real(ep) :: r
        r = dasum_quad(n, x, incx)
    end function dasum

    function dnrm2(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        real(ep) :: r
        r = dnrm2_quad(n, x, incx)
    end function dnrm2

    function idamax(n, x, incx) result(r)
        integer,  intent(in) :: n, incx
        real(ep), intent(in) :: x(*)
        integer :: r
        r = idamax_quad(n, x, incx)
    end function idamax

    subroutine daxpy(n, alpha, x, incx, y, incy)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(in)    :: alpha, x(*)
        real(ep), intent(inout) :: y(*)
        call daxpy_quad(n, alpha, x, incx, y, incy)
    end subroutine daxpy

    subroutine dcopy(n, x, incx, y, incy)
        integer,  intent(in)  :: n, incx, incy
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: y(*)
        call dcopy_quad(n, x, incx, y, incy)
    end subroutine dcopy

    subroutine dscal(n, alpha, x, incx)
        integer,  intent(in)    :: n, incx
        real(ep), intent(in)    :: alpha
        real(ep), intent(inout) :: x(*)
        call dscal_quad(n, alpha, x, incx)
    end subroutine dscal

    subroutine dswap(n, x, incx, y, incy)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        call dswap_quad(n, x, incx, y, incy)
    end subroutine dswap

    subroutine drot(n, x, incx, y, incy, c, s)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        real(ep), intent(in)    :: c, s
        call drot_quad(n, x, incx, y, incy, c, s)
    end subroutine drot

    subroutine drotg(a, b, c, s)
        real(ep), intent(inout) :: a, b
        real(ep), intent(out)   :: c, s
        call drotg_quad(a, b, c, s)
    end subroutine drotg

    subroutine drotm(n, x, incx, y, incy, param)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(inout) :: x(*), y(*)
        real(ep), intent(in)    :: param(5)
        call drotm_quad(n, x, incx, y, incy, param)
    end subroutine drotm

    subroutine drotmg(d1, d2, x1, y1, param)
        real(ep), intent(inout) :: d1, d2, x1
        real(ep), intent(in)    :: y1
        real(ep), intent(out)   :: param(5)
        call drotmg_quad(d1, d2, x1, y1, param)
    end subroutine drotmg

    function zdotc(n, x, incx, y, incy) result(r)
        integer,     intent(in) :: n, incx, incy
        complex(ep), intent(in) :: x(*), y(*)
        complex(ep) :: r
        r = zdotc_quad(n, x, incx, y, incy)
    end function zdotc

    function zdotu(n, x, incx, y, incy) result(r)
        integer,     intent(in) :: n, incx, incy
        complex(ep), intent(in) :: x(*), y(*)
        complex(ep) :: r
        r = zdotu_quad(n, x, incx, y, incy)
    end function zdotu

    function dzasum(n, x, incx) result(r)
        integer,     intent(in) :: n, incx
        complex(ep), intent(in) :: x(*)
        real(ep) :: r
        r = dzasum_quad(n, x, incx)
    end function dzasum

    subroutine zaxpy(n, alpha, x, incx, y, incy)
        integer,     intent(in)    :: n, incx, incy
        complex(ep), intent(in)    :: alpha, x(*)
        complex(ep), intent(inout) :: y(*)
        call zaxpy_quad(n, alpha, x, incx, y, incy)
    end subroutine zaxpy

    subroutine zscal(n, alpha, x, incx)
        integer,     intent(in)    :: n, incx
        complex(ep), intent(in)    :: alpha
        complex(ep), intent(inout) :: x(*)
        call zscal_quad(n, alpha, x, incx)
    end subroutine zscal

    subroutine dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        call dgemv_quad(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine dgemv

    subroutine dgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: trans
        integer,   intent(in)    :: m, n, kl, ku, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        call dgbmv_quad(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine dgbmv

    subroutine dger(m, n, alpha, x, incx, y, incy, A, lda)
        integer,  intent(in)    :: m, n, incx, incy, lda
        real(ep), intent(in)    :: alpha, x(*), y(*)
        real(ep), intent(inout) :: A(lda,*)
        call dger_quad(m, n, alpha, x, incx, y, incy, A, lda)
    end subroutine dger

    subroutine dsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        call dsymv_quad(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine dsymv

    subroutine dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, incx, incy
        real(ep),  intent(in)    :: alpha, beta, ap(*), x(*)
        real(ep),  intent(inout) :: y(*)
        call dspmv_quad(uplo, n, alpha, ap, x, incx, beta, y, incy)
    end subroutine dspmv

    subroutine dsbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, k, lda, incx, incy
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), x(*)
        real(ep),  intent(inout) :: y(*)
        call dsbmv_quad(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine dsbmv

    subroutine dtrmv(uplo, trans, diag, n, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        call dtrmv_quad(uplo, trans, diag, n, A, lda, x, incx)
    end subroutine dtrmv

    subroutine dtbmv(uplo, trans, diag, n, k, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, k, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        call dtbmv_quad(uplo, trans, diag, n, k, A, lda, x, incx)
    end subroutine dtbmv

    subroutine dtpmv(uplo, trans, diag, n, ap, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, incx
        real(ep),  intent(in)    :: ap(*)
        real(ep),  intent(inout) :: x(*)
        call dtpmv_quad(uplo, trans, diag, n, ap, x, incx)
    end subroutine dtpmv

    subroutine dtrsv(uplo, trans, diag, n, A, lda, x, incx)
        character, intent(in)    :: uplo, trans, diag
        integer,   intent(in)    :: n, lda, incx
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: x(*)
        call dtrsv_quad(uplo, trans, diag, n, A, lda, x, incx)
    end subroutine dtrsv

    subroutine zgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
        character,   intent(in)    :: trans
        integer,     intent(in)    :: m, n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        call zgemv_quad(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine zgemv

    subroutine zhemv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        call zhemv_quad(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
    end subroutine zhemv

    subroutine zgerc(m, n, alpha, x, incx, y, incy, A, lda)
        integer,     intent(in)    :: m, n, incx, incy, lda
        complex(ep), intent(in)    :: alpha, x(*), y(*)
        complex(ep), intent(inout) :: A(lda,*)
        call zgerc_quad(m, n, alpha, x, incx, y, incy, A, lda)
    end subroutine zgerc

    subroutine dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: transa, transb
        integer,   intent(in)    :: m, n, k, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        call dgemm_quad(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    end subroutine dgemm

    subroutine dsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: side, uplo
        integer,   intent(in)    :: m, n, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        call dsymm_quad(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
    end subroutine dsymm

    subroutine dsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
        character, intent(in)    :: uplo, trans
        integer,   intent(in)    :: n, k, lda, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: C(ldc,*)
        call dsyrk_quad(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
    end subroutine dsyrk

    subroutine dsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character, intent(in)    :: uplo, trans
        integer,   intent(in)    :: n, k, lda, ldb, ldc
        real(ep),  intent(in)    :: alpha, beta
        real(ep),  intent(in)    :: A(lda,*), B(ldb,*)
        real(ep),  intent(inout) :: C(ldc,*)
        call dsyr2k_quad(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    end subroutine dsyr2k

    subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character, intent(in)    :: side, uplo, transa, diag
        integer,   intent(in)    :: m, n, lda, ldb
        real(ep),  intent(in)    :: alpha
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        call dtrmm_quad(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
    end subroutine dtrmm

    subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character, intent(in)    :: side, uplo, transa, diag
        integer,   intent(in)    :: m, n, lda, ldb
        real(ep),  intent(in)    :: alpha
        real(ep),  intent(in)    :: A(lda,*)
        real(ep),  intent(inout) :: B(ldb,*)
        call dtrsm_quad(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
    end subroutine dtrsm

    subroutine zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
        character,   intent(in)    :: transa, transb
        integer,     intent(in)    :: m, n, k, lda, ldb, ldc
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
        complex(ep), intent(inout) :: C(ldc,*)
        call zgemm_quad(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
    end subroutine zgemm

    subroutine zhemm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
        character,   intent(in)    :: side, uplo
        integer,     intent(in)    :: m, n, lda, ldb, ldc
        complex(ep), intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*), B(ldb,*)
        complex(ep), intent(inout) :: C(ldc,*)
        call zhemm_quad(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc)
    end subroutine zhemm

    subroutine zherk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
        character,   intent(in)    :: uplo, trans
        integer,     intent(in)    :: n, k, lda, ldc
        real(ep),    intent(in)    :: alpha, beta
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: C(ldc,*)
        call zherk_quad(uplo, trans, n, k, alpha, A, lda, beta, C, ldc)
    end subroutine zherk

    subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
        character,   intent(in)    :: side, uplo, transa, diag
        integer,     intent(in)    :: m, n, lda, ldb
        complex(ep), intent(in)    :: alpha
        complex(ep), intent(in)    :: A(lda,*)
        complex(ep), intent(inout) :: B(ldb,*)
        call ztrsm_quad(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb)
    end subroutine ztrsm


end module ref_quad_blas
