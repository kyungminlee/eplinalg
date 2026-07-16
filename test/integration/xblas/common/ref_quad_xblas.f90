module ref_quad_xblas
    !! Quad-precision reference implementations for XBLAS routines.
    !!
    !! For routines with a direct standard-BLAS analogue (gemv, gbmv,
    !! symv, sbmv, spmv, hemv, hbmv, hpmv, trmv, tpmv, tbsv, trsv,
    !! gemm, symm, hemm) we forward to the ``-freal-8-real-16``-promoted
    !! Netlib routine in ``refblas_quad`` and apply XBLAS's extra
    !! alpha-scaling for the triangular-only-input routines.  For
    !! XBLAS-only operations (axpby, sum, waxpby, gemv2/gbmv2/symv2/
    !! hemv2 head-tail forms, ge_sum_mv, complex symmetric matvec) we
    !! hand-code at quad precision.
    use prec_kinds, only: ep
    implicit none
    private

    public :: ref_blas_daxpby_x, ref_blas_zaxpby_x
    public :: ref_blas_ddot_x,   ref_blas_zdot_x
    public :: ref_blas_dsum_x,   ref_blas_zsum_x
    public :: ref_blas_dwaxpby_x, ref_blas_zwaxpby_x

    public :: ref_blas_dgemv_x,  ref_blas_zgemv_x
    public :: ref_blas_dgemv2_x, ref_blas_zgemv2_x
    public :: ref_blas_dgbmv_x,  ref_blas_zgbmv_x
    public :: ref_blas_dgbmv2_x, ref_blas_zgbmv2_x
    public :: ref_blas_dsymv_x,  ref_blas_zsymv_x
    public :: ref_blas_dsymv2_x, ref_blas_zsymv2_x
    public :: ref_blas_dsbmv_x,  ref_blas_zsbmv_x
    public :: ref_blas_dspmv_x,  ref_blas_zspmv_x
    public :: ref_blas_zhemv_x,  ref_blas_zhemv2_x
    public :: ref_blas_zhbmv_x,  ref_blas_zhpmv_x
    public :: ref_blas_dge_sum_mv_x, ref_blas_zge_sum_mv_x
    public :: ref_blas_dtrmv_x,  ref_blas_ztrmv_x
    public :: ref_blas_dtpmv_x,  ref_blas_ztpmv_x
    public :: ref_blas_dtbsv_x,  ref_blas_ztbsv_x
    public :: ref_blas_dtrsv_x,  ref_blas_ztrsv_x

    public :: ref_blas_dgemm_x,  ref_blas_zgemm_x
    public :: ref_blas_dsymm_x,  ref_blas_zsymm_x
    public :: ref_blas_zhemm_x

    integer, parameter :: BLAS_NO_CONJ      = 192
    integer, parameter :: BLAS_CONJ_REF     = 191
    integer, parameter :: BLAS_NO_TRANS_REF = 111
    integer, parameter :: BLAS_TRANS_REF    = 112
    integer, parameter :: BLAS_CONJ_TRANS_REF = 113
    integer, parameter :: BLAS_UPPER_REF    = 121
    integer, parameter :: BLAS_LOWER_REF    = 122
    integer, parameter :: BLAS_NON_UNIT_DIAG_REF = 131
    integer, parameter :: BLAS_LEFT_SIDE_REF = 141

    ! ── refblas_quad declarations ───────────────────────────────────
    interface
        subroutine dgemv_quad(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in) :: trans
            integer,   intent(in) :: m, n, lda, incx, incy
            real(ep),  intent(in) :: alpha, beta, a(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine zgemv_quad(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import :: ep
            character,    intent(in) :: trans
            integer,      intent(in) :: m, n, lda, incx, incy
            complex(ep),  intent(in) :: alpha, beta, a(lda,*), x(*)
            complex(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine dgbmv_quad(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in) :: trans
            integer,   intent(in) :: m, n, kl, ku, lda, incx, incy
            real(ep),  intent(in) :: alpha, beta, a(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine zgbmv_quad(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
            import :: ep
            character,    intent(in) :: trans
            integer,      intent(in) :: m, n, kl, ku, lda, incx, incy
            complex(ep),  intent(in) :: alpha, beta, a(lda,*), x(*)
            complex(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine dsymv_quad(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in) :: uplo
            integer,   intent(in) :: n, lda, incx, incy
            real(ep),  intent(in) :: alpha, beta, a(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine zhemv_quad(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
            import :: ep
            character,    intent(in) :: uplo
            integer,      intent(in) :: n, lda, incx, incy
            complex(ep),  intent(in) :: alpha, beta, a(lda,*), x(*)
            complex(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine dsbmv_quad(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
            import :: ep
            character, intent(in) :: uplo
            integer,   intent(in) :: n, k, lda, incx, incy
            real(ep),  intent(in) :: alpha, beta, a(lda,*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine zhbmv_quad(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
            import :: ep
            character,    intent(in) :: uplo
            integer,      intent(in) :: n, k, lda, incx, incy
            complex(ep),  intent(in) :: alpha, beta, a(lda,*), x(*)
            complex(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine dspmv_quad(uplo, n, alpha, ap, x, incx, beta, y, incy)
            import :: ep
            character, intent(in) :: uplo
            integer,   intent(in) :: n, incx, incy
            real(ep),  intent(in) :: alpha, beta, ap(*), x(*)
            real(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine zhpmv_quad(uplo, n, alpha, ap, x, incx, beta, y, incy)
            import :: ep
            character,    intent(in) :: uplo
            integer,      intent(in) :: n, incx, incy
            complex(ep),  intent(in) :: alpha, beta, ap(*), x(*)
            complex(ep),  intent(inout) :: y(*)
        end subroutine
        subroutine dtrmv_quad(uplo, trans, diag, n, a, lda, x, incx)
            import :: ep
            character, intent(in) :: uplo, trans, diag
            integer,   intent(in) :: n, lda, incx
            real(ep),  intent(in) :: a(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine
        subroutine ztrmv_quad(uplo, trans, diag, n, a, lda, x, incx)
            import :: ep
            character,    intent(in) :: uplo, trans, diag
            integer,      intent(in) :: n, lda, incx
            complex(ep),  intent(in) :: a(lda,*)
            complex(ep),  intent(inout) :: x(*)
        end subroutine
        subroutine dtpmv_quad(uplo, trans, diag, n, ap, x, incx)
            import :: ep
            character, intent(in) :: uplo, trans, diag
            integer,   intent(in) :: n, incx
            real(ep),  intent(in) :: ap(*)
            real(ep),  intent(inout) :: x(*)
        end subroutine
        subroutine ztpmv_quad(uplo, trans, diag, n, ap, x, incx)
            import :: ep
            character,    intent(in) :: uplo, trans, diag
            integer,      intent(in) :: n, incx
            complex(ep),  intent(in) :: ap(*)
            complex(ep),  intent(inout) :: x(*)
        end subroutine
        subroutine dtbsv_quad(uplo, trans, diag, n, k, a, lda, x, incx)
            import :: ep
            character, intent(in) :: uplo, trans, diag
            integer,   intent(in) :: n, k, lda, incx
            real(ep),  intent(in) :: a(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine
        subroutine ztbsv_quad(uplo, trans, diag, n, k, a, lda, x, incx)
            import :: ep
            character,    intent(in) :: uplo, trans, diag
            integer,      intent(in) :: n, k, lda, incx
            complex(ep),  intent(in) :: a(lda,*)
            complex(ep),  intent(inout) :: x(*)
        end subroutine
        subroutine dtrsv_quad(uplo, trans, diag, n, a, lda, x, incx)
            import :: ep
            character, intent(in) :: uplo, trans, diag
            integer,   intent(in) :: n, lda, incx
            real(ep),  intent(in) :: a(lda,*)
            real(ep),  intent(inout) :: x(*)
        end subroutine
        subroutine ztrsv_quad(uplo, trans, diag, n, a, lda, x, incx)
            import :: ep
            character,    intent(in) :: uplo, trans, diag
            integer,      intent(in) :: n, lda, incx
            complex(ep),  intent(in) :: a(lda,*)
            complex(ep),  intent(inout) :: x(*)
        end subroutine
        subroutine dgemm_quad(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: ep
            character, intent(in) :: transa, transb
            integer,   intent(in) :: m, n, k, lda, ldb, ldc
            real(ep),  intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            real(ep),  intent(inout) :: c(ldc,*)
        end subroutine
        subroutine zgemm_quad(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: ep
            character,    intent(in) :: transa, transb
            integer,      intent(in) :: m, n, k, lda, ldb, ldc
            complex(ep),  intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            complex(ep),  intent(inout) :: c(ldc,*)
        end subroutine
        subroutine dsymm_quad(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: ep
            character, intent(in) :: side, uplo
            integer,   intent(in) :: m, n, lda, ldb, ldc
            real(ep),  intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            real(ep),  intent(inout) :: c(ldc,*)
        end subroutine
        subroutine zsymm_quad(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: ep
            character,    intent(in) :: side, uplo
            integer,      intent(in) :: m, n, lda, ldb, ldc
            complex(ep),  intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            complex(ep),  intent(inout) :: c(ldc,*)
        end subroutine
        subroutine zhemm_quad(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: ep
            character,    intent(in) :: side, uplo
            integer,      intent(in) :: m, n, lda, ldb, ldc
            complex(ep),  intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            complex(ep),  intent(inout) :: c(ldc,*)
        end subroutine
    end interface

contains

    ! ── XBLAS-enum → BLAS-character helpers ─────────────────────────
    pure function trans_char(t) result(c)
        integer, intent(in) :: t
        character           :: c
        select case (t)
        case (111); c = 'N'
        case (112); c = 'T'
        case (113); c = 'C'
        case default; c = 'N'
        end select
    end function

    pure function uplo_char(u) result(c)
        integer, intent(in) :: u
        character           :: c
        select case (u)
        case (121); c = 'U'
        case (122); c = 'L'
        case default; c = 'U'
        end select
    end function

    pure function diag_char(d) result(c)
        integer, intent(in) :: d
        character           :: c
        select case (d)
        case (131); c = 'N'
        case (132); c = 'U'
        case default; c = 'N'
        end select
    end function

    pure function side_char(s) result(c)
        integer, intent(in) :: s
        character           :: c
        select case (s)
        case (141); c = 'L'
        case (142); c = 'R'
        case default; c = 'L'
        end select
    end function

    pure function ix_start(n, inc) result(ix)
        integer, intent(in) :: n, inc
        integer             :: ix
        ix = 1
        if (inc < 0) ix = 1 - (n - 1) * inc
    end function

    ! ════════════════════════════════════════════════════════════════
    ! Level 1
    ! ════════════════════════════════════════════════════════════════

    pure subroutine ref_blas_daxpby_x(n, alpha, x, incx, beta, y, incy)
        integer,  intent(in)    :: n, incx, incy
        real(ep), intent(in)    :: alpha, beta, x(*)
        real(ep), intent(inout) :: y(*)
        integer :: i, ix, iy
        ix = ix_start(n, incx); iy = ix_start(n, incy)
        do i = 1, n
            y(iy) = alpha * x(ix) + beta * y(iy)
            ix = ix + incx; iy = iy + incy
        end do
    end subroutine

    pure subroutine ref_blas_zaxpby_x(n, alpha, x, incx, beta, y, incy)
        integer,     intent(in)    :: n, incx, incy
        complex(ep), intent(in)    :: alpha, beta, x(*)
        complex(ep), intent(inout) :: y(*)
        integer :: i, ix, iy
        ix = ix_start(n, incx); iy = ix_start(n, incy)
        do i = 1, n
            y(iy) = alpha * x(ix) + beta * y(iy)
            ix = ix + incx; iy = iy + incy
        end do
    end subroutine

    pure subroutine ref_blas_ddot_x(conj, n, alpha, x, incx, beta, y, incy, r)
        integer,  intent(in)    :: conj, n, incx, incy
        real(ep), intent(in)    :: alpha, beta, x(*), y(*)
        real(ep), intent(inout) :: r
        real(ep) :: acc
        integer  :: i, ix, iy
        acc = 0.0_ep
        ix = ix_start(n, incx); iy = ix_start(n, incy)
        do i = 1, n
            acc = acc + x(ix) * y(iy)
            ix = ix + incx; iy = iy + incy
        end do
        r = alpha * acc + beta * r
    end subroutine

    pure subroutine ref_blas_zdot_x(conj, n, alpha, x, incx, beta, y, incy, r)
        integer,     intent(in)    :: conj, n, incx, incy
        complex(ep), intent(in)    :: alpha, beta, x(*), y(*)
        complex(ep), intent(inout) :: r
        complex(ep) :: acc
        integer     :: i, ix, iy
        acc = (0.0_ep, 0.0_ep)
        ix = ix_start(n, incx); iy = ix_start(n, incy)
        if (conj == BLAS_CONJ_REF) then
            do i = 1, n
                acc = acc + conjg(x(ix)) * y(iy)
                ix = ix + incx; iy = iy + incy
            end do
        else
            do i = 1, n
                acc = acc + x(ix) * y(iy)
                ix = ix + incx; iy = iy + incy
            end do
        end if
        r = alpha * acc + beta * r
    end subroutine

    pure subroutine ref_blas_dsum_x(n, x, incx, s)
        integer,  intent(in)  :: n, incx
        real(ep), intent(in)  :: x(*)
        real(ep), intent(out) :: s
        integer :: i, ix
        s = 0.0_ep
        ix = ix_start(n, incx)
        do i = 1, n
            s = s + x(ix)
            ix = ix + incx
        end do
    end subroutine

    pure subroutine ref_blas_zsum_x(n, x, incx, s)
        integer,     intent(in)  :: n, incx
        complex(ep), intent(in)  :: x(*)
        complex(ep), intent(out) :: s
        integer :: i, ix
        s = (0.0_ep, 0.0_ep)
        ix = ix_start(n, incx)
        do i = 1, n
            s = s + x(ix)
            ix = ix + incx
        end do
    end subroutine

    pure subroutine ref_blas_dwaxpby_x(n, alpha, x, incx, beta, y, incy, w, incw)
        integer,  intent(in)  :: n, incx, incy, incw
        real(ep), intent(in)  :: alpha, beta, x(*), y(*)
        real(ep), intent(out) :: w(*)
        integer :: i, ix, iy, iw
        ix = ix_start(n, incx); iy = ix_start(n, incy); iw = ix_start(n, incw)
        do i = 1, n
            w(iw) = alpha * x(ix) + beta * y(iy)
            ix = ix + incx; iy = iy + incy; iw = iw + incw
        end do
    end subroutine

    pure subroutine ref_blas_zwaxpby_x(n, alpha, x, incx, beta, y, incy, w, incw)
        integer,     intent(in)  :: n, incx, incy, incw
        complex(ep), intent(in)  :: alpha, beta, x(*), y(*)
        complex(ep), intent(out) :: w(*)
        integer :: i, ix, iy, iw
        ix = ix_start(n, incx); iy = ix_start(n, incy); iw = ix_start(n, incw)
        do i = 1, n
            w(iw) = alpha * x(ix) + beta * y(iy)
            ix = ix + incx; iy = iy + incy; iw = iw + incw
        end do
    end subroutine

    ! ════════════════════════════════════════════════════════════════
    ! Level 2 — gemv / gemv2
    ! ════════════════════════════════════════════════════════════════

    subroutine ref_blas_dgemv_x(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        integer,  intent(in)    :: trans, m, n, lda, incx, incy
        real(ep), intent(in)    :: alpha, beta, a(lda,*), x(*)
        real(ep), intent(inout) :: y(*)
        call dgemv_quad(trans_char(trans), m, n, alpha, a, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine ref_blas_zgemv_x(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        integer,     intent(in)    :: trans, m, n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        call zgemv_quad(trans_char(trans), m, n, alpha, a, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine ref_blas_dgemv2_x(trans, m, n, alpha, a, lda, head_x, tail_x, &
                                 incx, beta, y, incy)
        integer,  intent(in)    :: trans, m, n, lda, incx, incy
        real(ep), intent(in)    :: alpha, beta, a(lda,*), head_x(*), tail_x(*)
        real(ep), intent(inout) :: y(*)
        ! y := alpha*A*head_x + (1*y + beta*y_init - 1*y) ... easier:
        ! y := alpha*A*(head + tail) + beta*y
        ! Run dgemv twice: y := alpha*A*head + beta*y, then y := alpha*A*tail + 1*y.
        call dgemv_quad(trans_char(trans), m, n, alpha, a, lda, head_x, incx, beta, y, incy)
        call dgemv_quad(trans_char(trans), m, n, alpha, a, lda, tail_x, incx, 1.0_ep, y, incy)
    end subroutine

    subroutine ref_blas_zgemv2_x(trans, m, n, alpha, a, lda, head_x, tail_x, &
                                 incx, beta, y, incy)
        integer,     intent(in)    :: trans, m, n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), head_x(*), tail_x(*)
        complex(ep), intent(inout) :: y(*)
        call zgemv_quad(trans_char(trans), m, n, alpha, a, lda, head_x, incx, beta, y, incy)
        call zgemv_quad(trans_char(trans), m, n, alpha, a, lda, tail_x, incx, &
                   (1.0_ep, 0.0_ep), y, incy)
    end subroutine

    ! ════════════════════════════════════════════════════════════════
    ! Level 2 — gbmv / gbmv2
    ! ════════════════════════════════════════════════════════════════

    subroutine ref_blas_dgbmv_x(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
        integer,  intent(in)    :: trans, m, n, kl, ku, lda, incx, incy
        real(ep), intent(in)    :: alpha, beta, a(lda,*), x(*)
        real(ep), intent(inout) :: y(*)
        call dgbmv_quad(trans_char(trans), m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine ref_blas_zgbmv_x(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
        integer,     intent(in)    :: trans, m, n, kl, ku, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        call zgbmv_quad(trans_char(trans), m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine ref_blas_dgbmv2_x(trans, m, n, kl, ku, alpha, a, lda, head_x, tail_x, &
                                 incx, beta, y, incy)
        integer,  intent(in)    :: trans, m, n, kl, ku, lda, incx, incy
        real(ep), intent(in)    :: alpha, beta, a(lda,*), head_x(*), tail_x(*)
        real(ep), intent(inout) :: y(*)
        call dgbmv_quad(trans_char(trans), m, n, kl, ku, alpha, a, lda, head_x, incx, beta, y, incy)
        call dgbmv_quad(trans_char(trans), m, n, kl, ku, alpha, a, lda, tail_x, incx, 1.0_ep, y, incy)
    end subroutine

    subroutine ref_blas_zgbmv2_x(trans, m, n, kl, ku, alpha, a, lda, head_x, tail_x, &
                                 incx, beta, y, incy)
        integer,     intent(in)    :: trans, m, n, kl, ku, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), head_x(*), tail_x(*)
        complex(ep), intent(inout) :: y(*)
        call zgbmv_quad(trans_char(trans), m, n, kl, ku, alpha, a, lda, head_x, incx, beta, y, incy)
        call zgbmv_quad(trans_char(trans), m, n, kl, ku, alpha, a, lda, tail_x, incx, &
                   (1.0_ep, 0.0_ep), y, incy)
    end subroutine

    ! ════════════════════════════════════════════════════════════════
    ! Level 2 — symv (real and complex), symv2, sbmv, spmv
    ! ════════════════════════════════════════════════════════════════

    subroutine ref_blas_dsymv_x(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
        integer,  intent(in)    :: uplo, n, lda, incx, incy
        real(ep), intent(in)    :: alpha, beta, a(lda,*), x(*)
        real(ep), intent(inout) :: y(*)
        call dsymv_quad(uplo_char(uplo), n, alpha, a, lda, x, incx, beta, y, incy)
    end subroutine

    ! Complex symmetric matvec — no BLAS analogue (BLAS has only Hermitian).
    ! Hand-coded: y := alpha * A * x + beta * y, where A is symmetric (A = A^T)
    ! with only the indicated triangle stored.
    pure subroutine ref_blas_zsymv_x(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
        integer,     intent(in)    :: uplo, n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        complex(ep) :: acc, aij
        integer     :: i, j, ix, iy, jx
        ix = ix_start(n, incx); iy = ix_start(n, incy)
        do i = 1, n
            y(iy + (i - 1) * incy) = beta * y(iy + (i - 1) * incy)
        end do
        do i = 1, n
            acc = (0.0_ep, 0.0_ep)
            jx = ix
            do j = 1, n
                if (uplo == BLAS_UPPER_REF) then
                    if (i <= j) then; aij = a(i, j); else; aij = a(j, i); end if
                else
                    if (i >= j) then; aij = a(i, j); else; aij = a(j, i); end if
                end if
                acc = acc + aij * x(jx)
                jx = jx + incx
            end do
            y(iy + (i - 1) * incy) = y(iy + (i - 1) * incy) + alpha * acc
        end do
    end subroutine

    subroutine ref_blas_dsymv2_x(uplo, n, alpha, a, lda, x_head, x_tail, incx, &
                                 beta, y, incy)
        integer,  intent(in)    :: uplo, n, lda, incx, incy
        real(ep), intent(in)    :: alpha, beta, a(lda,*), x_head(*), x_tail(*)
        real(ep), intent(inout) :: y(*)
        call dsymv_quad(uplo_char(uplo), n, alpha, a, lda, x_head, incx, beta, y, incy)
        call dsymv_quad(uplo_char(uplo), n, alpha, a, lda, x_tail, incx, 1.0_ep, y, incy)
    end subroutine

    pure subroutine ref_blas_zsymv2_x(uplo, n, alpha, a, lda, x_head, x_tail, incx, &
                                      beta, y, incy)
        integer,     intent(in)    :: uplo, n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), x_head(*), x_tail(*)
        complex(ep), intent(inout) :: y(*)
        complex(ep) :: acc, aij, xij
        integer     :: i, j, ix, iy, jx
        ix = ix_start(n, incx); iy = ix_start(n, incy)
        do i = 1, n
            y(iy + (i - 1) * incy) = beta * y(iy + (i - 1) * incy)
        end do
        do i = 1, n
            acc = (0.0_ep, 0.0_ep)
            jx = ix
            do j = 1, n
                if (uplo == BLAS_UPPER_REF) then
                    if (i <= j) then; aij = a(i, j); else; aij = a(j, i); end if
                else
                    if (i >= j) then; aij = a(i, j); else; aij = a(j, i); end if
                end if
                xij = x_head(jx) + x_tail(jx)
                acc = acc + aij * xij
                jx = jx + incx
            end do
            y(iy + (i - 1) * incy) = y(iy + (i - 1) * incy) + alpha * acc
        end do
    end subroutine

    subroutine ref_blas_dsbmv_x(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
        integer,  intent(in)    :: uplo, n, k, lda, incx, incy
        real(ep), intent(in)    :: alpha, beta, a(lda,*), x(*)
        real(ep), intent(inout) :: y(*)
        call dsbmv_quad(uplo_char(uplo), n, k, alpha, a, lda, x, incx, beta, y, incy)
    end subroutine

    ! Complex symmetric banded — hand-code (no zsbmv in BLAS).
    pure subroutine ref_blas_zsbmv_x(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
        integer,     intent(in)    :: uplo, n, k, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        complex(ep) :: acc, aij
        integer     :: i, j, ix, iy, jx
        ix = ix_start(n, incx); iy = ix_start(n, incy)
        do i = 1, n
            y(iy + (i - 1) * incy) = beta * y(iy + (i - 1) * incy)
        end do
        do i = 1, n
            acc = (0.0_ep, 0.0_ep)
            jx = ix
            do j = 1, n
                if (abs(i - j) > k) then
                    aij = (0.0_ep, 0.0_ep)
                else if (uplo == BLAS_UPPER_REF) then
                    ! Upper banded: row stored at a(k+1+i-j, j) for i<=j;
                    ! symmetric mirror for i>j
                    if (i <= j) then
                        aij = a(k + 1 + i - j, j)
                    else
                        aij = a(k + 1 + j - i, i)
                    end if
                else
                    ! Lower banded: row stored at a(1+i-j, j) for i>=j
                    if (i >= j) then
                        aij = a(1 + i - j, j)
                    else
                        aij = a(1 + j - i, i)
                    end if
                end if
                acc = acc + aij * x(jx)
                jx = jx + incx
            end do
            y(iy + (i - 1) * incy) = y(iy + (i - 1) * incy) + alpha * acc
        end do
    end subroutine

    subroutine ref_blas_dspmv_x(uplo, n, alpha, ap, x, incx, beta, y, incy)
        integer,  intent(in)    :: uplo, n, incx, incy
        real(ep), intent(in)    :: alpha, beta, ap(*), x(*)
        real(ep), intent(inout) :: y(*)
        call dspmv_quad(uplo_char(uplo), n, alpha, ap, x, incx, beta, y, incy)
    end subroutine

    pure subroutine ref_blas_zspmv_x(uplo, n, alpha, ap, x, incx, beta, y, incy)
        integer,     intent(in)    :: uplo, n, incx, incy
        complex(ep), intent(in)    :: alpha, beta, ap(*), x(*)
        complex(ep), intent(inout) :: y(*)
        complex(ep) :: acc, aij
        integer     :: i, j, ix, iy, jx, idx
        ix = ix_start(n, incx); iy = ix_start(n, incy)
        do i = 1, n
            y(iy + (i - 1) * incy) = beta * y(iy + (i - 1) * incy)
        end do
        do i = 1, n
            acc = (0.0_ep, 0.0_ep)
            jx = ix
            do j = 1, n
                if (uplo == BLAS_UPPER_REF) then
                    ! Upper packed (column-major): a(i,j) at idx = i + j*(j-1)/2 for i<=j
                    if (i <= j) then
                        idx = i + j * (j - 1) / 2
                    else
                        idx = j + i * (i - 1) / 2
                    end if
                else
                    ! Lower packed: a(i,j) at idx = i + (j-1)*(2*n-j)/2 for i>=j
                    if (i >= j) then
                        idx = i + (j - 1) * (2 * n - j) / 2
                    else
                        idx = j + (i - 1) * (2 * n - i) / 2
                    end if
                end if
                aij = ap(idx)
                acc = acc + aij * x(jx)
                jx = jx + incx
            end do
            y(iy + (i - 1) * incy) = y(iy + (i - 1) * incy) + alpha * acc
        end do
    end subroutine

    ! ════════════════════════════════════════════════════════════════
    ! Level 2 — hemv / hemv2 / hbmv / hpmv (complex hermitian)
    ! ════════════════════════════════════════════════════════════════

    subroutine ref_blas_zhemv_x(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
        integer,     intent(in)    :: uplo, n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        call zhemv_quad(uplo_char(uplo), n, alpha, a, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine ref_blas_zhemv2_x(uplo, n, alpha, a, lda, x_head, x_tail, incx, &
                                 beta, y, incy)
        integer,     intent(in)    :: uplo, n, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), x_head(*), x_tail(*)
        complex(ep), intent(inout) :: y(*)
        call zhemv_quad(uplo_char(uplo), n, alpha, a, lda, x_head, incx, beta, y, incy)
        call zhemv_quad(uplo_char(uplo), n, alpha, a, lda, x_tail, incx, &
                   (1.0_ep, 0.0_ep), y, incy)
    end subroutine

    subroutine ref_blas_zhbmv_x(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
        integer,     intent(in)    :: uplo, n, k, lda, incx, incy
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), x(*)
        complex(ep), intent(inout) :: y(*)
        call zhbmv_quad(uplo_char(uplo), n, k, alpha, a, lda, x, incx, beta, y, incy)
    end subroutine

    subroutine ref_blas_zhpmv_x(uplo, n, alpha, ap, x, incx, beta, y, incy)
        integer,     intent(in)    :: uplo, n, incx, incy
        complex(ep), intent(in)    :: alpha, beta, ap(*), x(*)
        complex(ep), intent(inout) :: y(*)
        call zhpmv_quad(uplo_char(uplo), n, alpha, ap, x, incx, beta, y, incy)
    end subroutine

    ! ════════════════════════════════════════════════════════════════
    ! Level 2 — ge_sum_mv (y := alpha*A*x + beta*B*x)
    ! ════════════════════════════════════════════════════════════════

    subroutine ref_blas_dge_sum_mv_x(m, n, alpha, a, lda, x, incx, beta, b, ldb, y, incy)
        integer,  intent(in)    :: m, n, lda, ldb, incx, incy
        real(ep), intent(in)    :: alpha, beta, a(lda,*), b(ldb,*), x(*)
        real(ep), intent(inout) :: y(*)
        ! y := alpha*A*x   (with beta=0 init)
        call dgemv_quad('N', m, n, alpha, a, lda, x, incx, 0.0_ep, y, incy)
        ! y := beta*B*x + y
        call dgemv_quad('N', m, n, beta, b, ldb, x, incx, 1.0_ep, y, incy)
    end subroutine

    subroutine ref_blas_zge_sum_mv_x(m, n, alpha, a, lda, x, incx, beta, b, ldb, y, incy)
        integer,     intent(in)    :: m, n, lda, ldb, incx, incy
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), b(ldb,*), x(*)
        complex(ep), intent(inout) :: y(*)
        call zgemv_quad('N', m, n, alpha, a, lda, x, incx, (0.0_ep, 0.0_ep), y, incy)
        call zgemv_quad('N', m, n, beta,  b, ldb, x, incx, (1.0_ep, 0.0_ep), y, incy)
    end subroutine

    ! ════════════════════════════════════════════════════════════════
    ! Level 2 — triangular
    ! XBLAS adds an alpha pre-scaling. Reference: call BLAS to compute
    ! op(T)*x (or solve), then scale by alpha. (Standard BLAS trmv/trsv
    ! has no alpha argument.)
    ! ════════════════════════════════════════════════════════════════

    subroutine ref_blas_dtrmv_x(uplo, trans, diag, n, alpha, t, ldt, x, incx)
        integer,  intent(in)    :: uplo, trans, diag, n, ldt, incx
        real(ep), intent(in)    :: alpha, t(ldt,*)
        real(ep), intent(inout) :: x(*)
        integer :: lenx, i, ix
        lenx = 1 + (n - 1) * abs(incx)
        call dtrmv_quad(uplo_char(uplo), trans_char(trans), diag_char(diag), &
                   n, t, ldt, x, incx)
        ix = ix_start(n, incx)
        do i = 1, n
            x(ix) = alpha * x(ix)
            ix = ix + incx
        end do
    end subroutine

    subroutine ref_blas_ztrmv_x(uplo, trans, diag, n, alpha, t, ldt, x, incx)
        integer,     intent(in)    :: uplo, trans, diag, n, ldt, incx
        complex(ep), intent(in)    :: alpha, t(ldt,*)
        complex(ep), intent(inout) :: x(*)
        integer :: i, ix
        call ztrmv_quad(uplo_char(uplo), trans_char(trans), diag_char(diag), &
                   n, t, ldt, x, incx)
        ix = ix_start(n, incx)
        do i = 1, n
            x(ix) = alpha * x(ix)
            ix = ix + incx
        end do
    end subroutine

    subroutine ref_blas_dtpmv_x(uplo, trans, diag, n, alpha, tp, x, incx)
        integer,  intent(in)    :: uplo, trans, diag, n, incx
        real(ep), intent(in)    :: alpha, tp(*)
        real(ep), intent(inout) :: x(*)
        integer :: i, ix
        call dtpmv_quad(uplo_char(uplo), trans_char(trans), diag_char(diag), &
                   n, tp, x, incx)
        ix = ix_start(n, incx)
        do i = 1, n
            x(ix) = alpha * x(ix)
            ix = ix + incx
        end do
    end subroutine

    subroutine ref_blas_ztpmv_x(uplo, trans, diag, n, alpha, tp, x, incx)
        integer,     intent(in)    :: uplo, trans, diag, n, incx
        complex(ep), intent(in)    :: alpha, tp(*)
        complex(ep), intent(inout) :: x(*)
        integer :: i, ix
        call ztpmv_quad(uplo_char(uplo), trans_char(trans), diag_char(diag), &
                   n, tp, x, incx)
        ix = ix_start(n, incx)
        do i = 1, n
            x(ix) = alpha * x(ix)
            ix = ix + incx
        end do
    end subroutine

    subroutine ref_blas_dtbsv_x(uplo, trans, diag, n, k, alpha, t, ldt, x, incx)
        integer,  intent(in)    :: uplo, trans, diag, n, k, ldt, incx
        real(ep), intent(in)    :: alpha, t(ldt,*)
        real(ep), intent(inout) :: x(*)
        integer :: i, ix
        ! XBLAS: x := alpha * op(T)^-1 * x
        ! Equivalent: scale x by alpha, then solve.
        ix = ix_start(n, incx)
        do i = 1, n
            x(ix) = alpha * x(ix)
            ix = ix + incx
        end do
        call dtbsv_quad(uplo_char(uplo), trans_char(trans), diag_char(diag), &
                   n, k, t, ldt, x, incx)
    end subroutine

    subroutine ref_blas_ztbsv_x(uplo, trans, diag, n, k, alpha, t, ldt, x, incx)
        integer,     intent(in)    :: uplo, trans, diag, n, k, ldt, incx
        complex(ep), intent(in)    :: alpha, t(ldt,*)
        complex(ep), intent(inout) :: x(*)
        integer :: i, ix
        ix = ix_start(n, incx)
        do i = 1, n
            x(ix) = alpha * x(ix)
            ix = ix + incx
        end do
        call ztbsv_quad(uplo_char(uplo), trans_char(trans), diag_char(diag), &
                   n, k, t, ldt, x, incx)
    end subroutine

    subroutine ref_blas_dtrsv_x(uplo, trans, diag, n, alpha, t, ldt, x, incx)
        integer,  intent(in)    :: uplo, trans, diag, n, ldt, incx
        real(ep), intent(in)    :: alpha, t(ldt,*)
        real(ep), intent(inout) :: x(*)
        integer :: i, ix
        ix = ix_start(n, incx)
        do i = 1, n
            x(ix) = alpha * x(ix)
            ix = ix + incx
        end do
        call dtrsv_quad(uplo_char(uplo), trans_char(trans), diag_char(diag), &
                   n, t, ldt, x, incx)
    end subroutine

    subroutine ref_blas_ztrsv_x(uplo, trans, diag, n, alpha, t, ldt, x, incx)
        integer,     intent(in)    :: uplo, trans, diag, n, ldt, incx
        complex(ep), intent(in)    :: alpha, t(ldt,*)
        complex(ep), intent(inout) :: x(*)
        integer :: i, ix
        ix = ix_start(n, incx)
        do i = 1, n
            x(ix) = alpha * x(ix)
            ix = ix + incx
        end do
        call ztrsv_quad(uplo_char(uplo), trans_char(trans), diag_char(diag), &
                   n, t, ldt, x, incx)
    end subroutine

    ! ════════════════════════════════════════════════════════════════
    ! Level 3
    ! ════════════════════════════════════════════════════════════════

    subroutine ref_blas_dgemm_x(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                                beta, c, ldc)
        integer,  intent(in)    :: transa, transb, m, n, k, lda, ldb, ldc
        real(ep), intent(in)    :: alpha, beta, a(lda,*), b(ldb,*)
        real(ep), intent(inout) :: c(ldc,*)
        call dgemm_quad(trans_char(transa), trans_char(transb), m, n, k, &
                   alpha, a, lda, b, ldb, beta, c, ldc)
    end subroutine

    subroutine ref_blas_zgemm_x(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                                beta, c, ldc)
        integer,     intent(in)    :: transa, transb, m, n, k, lda, ldb, ldc
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), b(ldb,*)
        complex(ep), intent(inout) :: c(ldc,*)
        call zgemm_quad(trans_char(transa), trans_char(transb), m, n, k, &
                   alpha, a, lda, b, ldb, beta, c, ldc)
    end subroutine

    subroutine ref_blas_dsymm_x(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
        integer,  intent(in)    :: side, uplo, m, n, lda, ldb, ldc
        real(ep), intent(in)    :: alpha, beta, a(lda,*), b(ldb,*)
        real(ep), intent(inout) :: c(ldc,*)
        call dsymm_quad(side_char(side), uplo_char(uplo), m, n, &
                   alpha, a, lda, b, ldb, beta, c, ldc)
    end subroutine

    subroutine ref_blas_zsymm_x(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
        integer,     intent(in)    :: side, uplo, m, n, lda, ldb, ldc
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), b(ldb,*)
        complex(ep), intent(inout) :: c(ldc,*)
        call zsymm_quad(side_char(side), uplo_char(uplo), m, n, &
                   alpha, a, lda, b, ldb, beta, c, ldc)
    end subroutine

    subroutine ref_blas_zhemm_x(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
        integer,     intent(in)    :: side, uplo, m, n, lda, ldb, ldc
        complex(ep), intent(in)    :: alpha, beta, a(lda,*), b(ldb,*)
        complex(ep), intent(inout) :: c(ldc,*)
        call zhemm_quad(side_char(side), uplo_char(uplo), m, n, &
                   alpha, a, lda, b, ldb, beta, c, ldc)
    end subroutine

end module ref_quad_xblas
