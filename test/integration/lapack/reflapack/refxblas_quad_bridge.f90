! refxblas_quad_bridge — quad-precision implementations of the XBLAS
! Level-2 entry points called by the upstream LAPACK xx-family
! iterative-refinement helpers (dla_gerfsx_extended, dla_gbrfsx_extended,
! dla_porfsx_extended, dla_syrfsx_extended, and the corresponding z and
! complex-Hermitian variants).
!
! Compiled with -freal-8-real-16 so DOUBLE PRECISION / COMPLEX*16 are
! promoted to REAL/COMPLEX(KIND=16), matching the quad-precision call
! sites inside reflapack_quad. The implementations forward to the
! quad-promoted Netlib BLAS (refblas_quad) for everything that has a
! standard BLAS analogue and hand-code the complex-symmetric matvec
! (which has no upstream BLAS routine).
!
! At quad precision XBLAS's "extra-precise" residual refinement is the
! identity — quad working precision is already high enough that the
! head/tail decomposition collapses (head + tail in DD-on-quad becomes
! a single quad). For *2 routines we evaluate y := alpha*A*(head+tail)
! + beta*y by chaining two BLAS calls (one with beta scale, one with
! beta=1).
!
! Symbol naming and parameter list match the f2c bridge in
! ``extern/xblas-1.0.248/src/<routine>/BLAS_<routine>_x-f2c.c``: the
! Fortran-callable surface OMITS the C interface's leading ``order``
! parameter (the f2c bridge hardcodes ``blas_colmajor``) — every
! Fortran caller in LAPACK is column-major by convention. The trailing
! ``prec`` enum is accepted but ignored — the working precision is
! already quad. gfortran's default (no BIND(C)) emit produces
! lowercase + trailing-underscore symbols (blas_dgemv_x_,
! blas_zhemv2_x_, ...) matching reflapack_quad's call-site emission.

! ── Real-side bridges ───────────────────────────────────────────────

subroutine blas_dgemv_x(trans, m, n, alpha, A, lda, x, incx, beta, y, incy, prec)
    integer, intent(in) :: trans, m, n, lda, incx, incy, prec
    double precision, intent(in)    :: alpha, beta
    double precision, intent(in)    :: A(lda, *), x(*)
    double precision, intent(inout) :: y(*)
    character :: t
    select case (trans)
    case (111); t = 'N'
    case (112); t = 'T'
    case default; t = 'C'
    end select
    call dgemv(t, m, n, alpha, A, lda, x, incx, beta, y, incy)
end subroutine

subroutine blas_dgemv2_x(trans, m, n, alpha, A, lda, head_x, tail_x, &
                        incx, beta, y, incy, prec)
    integer, intent(in) :: trans, m, n, lda, incx, incy, prec
    double precision, intent(in)    :: alpha, beta
    double precision, intent(in)    :: A(lda, *), head_x(*), tail_x(*)
    double precision, intent(inout) :: y(*)
    character :: t
    select case (trans)
    case (111); t = 'N'
    case (112); t = 'T'
    case default; t = 'C'
    end select
    call dgemv(t, m, n, alpha, A, lda, head_x, incx, beta, y, incy)
    call dgemv(t, m, n, alpha, A, lda, tail_x, incx, 1.0d0, y, incy)
end subroutine

subroutine blas_dgbmv_x(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy, prec)
    integer, intent(in) :: trans, m, n, kl, ku, lda, incx, incy, prec
    double precision, intent(in)    :: alpha, beta
    double precision, intent(in)    :: A(lda, *), x(*)
    double precision, intent(inout) :: y(*)
    character :: t
    select case (trans)
    case (111); t = 'N'
    case (112); t = 'T'
    case default; t = 'C'
    end select
    call dgbmv(t, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
end subroutine

subroutine blas_dgbmv2_x(trans, m, n, kl, ku, alpha, A, lda, head_x, tail_x, &
                        incx, beta, y, incy, prec)
    integer, intent(in) :: trans, m, n, kl, ku, lda, incx, incy, prec
    double precision, intent(in)    :: alpha, beta
    double precision, intent(in)    :: A(lda, *), head_x(*), tail_x(*)
    double precision, intent(inout) :: y(*)
    character :: t
    select case (trans)
    case (111); t = 'N'
    case (112); t = 'T'
    case default; t = 'C'
    end select
    call dgbmv(t, m, n, kl, ku, alpha, A, lda, head_x, incx, beta, y, incy)
    call dgbmv(t, m, n, kl, ku, alpha, A, lda, tail_x, incx, 1.0d0, y, incy)
end subroutine

subroutine blas_dsymv_x(uplo, n, alpha, A, lda, x, incx, beta, y, incy, prec)
    integer, intent(in) :: uplo, n, lda, incx, incy, prec
    double precision, intent(in)    :: alpha, beta
    double precision, intent(in)    :: A(lda, *), x(*)
    double precision, intent(inout) :: y(*)
    character :: u
    if (uplo == 121) then
        u = 'U'
    else
        u = 'L'
    end if
    call dsymv(u, n, alpha, A, lda, x, incx, beta, y, incy)
end subroutine

subroutine blas_dsymv2_x(uplo, n, alpha, A, lda, x_head, x_tail, incx, &
                        beta, y, incy, prec)
    integer, intent(in) :: uplo, n, lda, incx, incy, prec
    double precision, intent(in)    :: alpha, beta
    double precision, intent(in)    :: A(lda, *), x_head(*), x_tail(*)
    double precision, intent(inout) :: y(*)
    character :: u
    if (uplo == 121) then
        u = 'U'
    else
        u = 'L'
    end if
    call dsymv(u, n, alpha, A, lda, x_head, incx, beta, y, incy)
    call dsymv(u, n, alpha, A, lda, x_tail, incx, 1.0d0, y, incy)
end subroutine

! ── Complex-side bridges ────────────────────────────────────────────

subroutine blas_zgemv_x(trans, m, n, alpha, A, lda, x, incx, beta, y, incy, prec)
    integer, intent(in) :: trans, m, n, lda, incx, incy, prec
    complex(kind=kind((1d0,0d0))), intent(in)    :: alpha, beta
    complex(kind=kind((1d0,0d0))), intent(in)    :: A(lda, *), x(*)
    complex(kind=kind((1d0,0d0))), intent(inout) :: y(*)
    character :: t
    select case (trans)
    case (111); t = 'N'
    case (112); t = 'T'
    case default; t = 'C'
    end select
    call zgemv(t, m, n, alpha, A, lda, x, incx, beta, y, incy)
end subroutine

subroutine blas_zgemv2_x(trans, m, n, alpha, A, lda, head_x, tail_x, &
                        incx, beta, y, incy, prec)
    integer, intent(in) :: trans, m, n, lda, incx, incy, prec
    complex(kind=kind((1d0,0d0))), intent(in)    :: alpha, beta
    complex(kind=kind((1d0,0d0))), intent(in)    :: A(lda, *), head_x(*), tail_x(*)
    complex(kind=kind((1d0,0d0))), intent(inout) :: y(*)
    character :: t
    select case (trans)
    case (111); t = 'N'
    case (112); t = 'T'
    case default; t = 'C'
    end select
    call zgemv(t, m, n, alpha, A, lda, head_x, incx, beta, y, incy)
    call zgemv(t, m, n, alpha, A, lda, tail_x, incx, (1.0d0, 0.0d0), y, incy)
end subroutine

subroutine blas_zgbmv_x(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy, prec)
    integer, intent(in) :: trans, m, n, kl, ku, lda, incx, incy, prec
    complex(kind=kind((1d0,0d0))), intent(in)    :: alpha, beta
    complex(kind=kind((1d0,0d0))), intent(in)    :: A(lda, *), x(*)
    complex(kind=kind((1d0,0d0))), intent(inout) :: y(*)
    character :: t
    select case (trans)
    case (111); t = 'N'
    case (112); t = 'T'
    case default; t = 'C'
    end select
    call zgbmv(t, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy)
end subroutine

subroutine blas_zgbmv2_x(trans, m, n, kl, ku, alpha, A, lda, head_x, tail_x, &
                        incx, beta, y, incy, prec)
    integer, intent(in) :: trans, m, n, kl, ku, lda, incx, incy, prec
    complex(kind=kind((1d0,0d0))), intent(in)    :: alpha, beta
    complex(kind=kind((1d0,0d0))), intent(in)    :: A(lda, *), head_x(*), tail_x(*)
    complex(kind=kind((1d0,0d0))), intent(inout) :: y(*)
    character :: t
    select case (trans)
    case (111); t = 'N'
    case (112); t = 'T'
    case default; t = 'C'
    end select
    call zgbmv(t, m, n, kl, ku, alpha, A, lda, head_x, incx, beta, y, incy)
    call zgbmv(t, m, n, kl, ku, alpha, A, lda, tail_x, incx, (1.0d0, 0.0d0), y, incy)
end subroutine

! Complex symmetric matvec — no BLAS analogue (BLAS has only Hermitian).
! Hand-code: y := alpha * A * x + beta * y. XBLAS callers from LAPACK
! always pass unit increments, so the increment handling is just a
! safety net.
subroutine blas_zsymv_x(uplo, n, alpha, A, lda, x, incx, beta, y, incy, prec)
    integer, intent(in) :: uplo, n, lda, incx, incy, prec
    complex(kind=kind((1d0,0d0))), intent(in)    :: alpha, beta
    complex(kind=kind((1d0,0d0))), intent(in)    :: A(lda, *), x(*)
    complex(kind=kind((1d0,0d0))), intent(inout) :: y(*)
    complex(kind=kind((1d0,0d0))) :: acc, aij
    integer :: i, j, ix, iy, jx
    ix = 1; if (incx < 0) ix = 1 - (n - 1) * incx
    iy = 1; if (incy < 0) iy = 1 - (n - 1) * incy
    do i = 1, n
        y(iy + (i - 1) * incy) = beta * y(iy + (i - 1) * incy)
    end do
    do i = 1, n
        acc = (0.0d0, 0.0d0)
        jx = ix
        do j = 1, n
            if (uplo == 121) then
                if (i <= j) then
                    aij = A(i, j)
                else
                    aij = A(j, i)
                end if
            else
                if (i >= j) then
                    aij = A(i, j)
                else
                    aij = A(j, i)
                end if
            end if
            acc = acc + aij * x(jx)
            jx = jx + incx
        end do
        y(iy + (i - 1) * incy) = y(iy + (i - 1) * incy) + alpha * acc
    end do
end subroutine

subroutine blas_zsymv2_x(uplo, n, alpha, A, lda, x_head, x_tail, incx, &
                        beta, y, incy, prec)
    integer, intent(in) :: uplo, n, lda, incx, incy, prec
    complex(kind=kind((1d0,0d0))), intent(in)    :: alpha, beta
    complex(kind=kind((1d0,0d0))), intent(in)    :: A(lda, *), x_head(*), x_tail(*)
    complex(kind=kind((1d0,0d0))), intent(inout) :: y(*)
    complex(kind=kind((1d0,0d0))) :: acc, aij, xij
    integer :: i, j, ix, iy, jx
    ix = 1; if (incx < 0) ix = 1 - (n - 1) * incx
    iy = 1; if (incy < 0) iy = 1 - (n - 1) * incy
    do i = 1, n
        y(iy + (i - 1) * incy) = beta * y(iy + (i - 1) * incy)
    end do
    do i = 1, n
        acc = (0.0d0, 0.0d0)
        jx = ix
        do j = 1, n
            if (uplo == 121) then
                if (i <= j) then
                    aij = A(i, j)
                else
                    aij = A(j, i)
                end if
            else
                if (i >= j) then
                    aij = A(i, j)
                else
                    aij = A(j, i)
                end if
            end if
            xij = x_head(jx) + x_tail(jx)
            acc = acc + aij * xij
            jx = jx + incx
        end do
        y(iy + (i - 1) * incy) = y(iy + (i - 1) * incy) + alpha * acc
    end do
end subroutine

subroutine blas_zhemv_x(uplo, n, alpha, A, lda, x, incx, beta, y, incy, prec)
    integer, intent(in) :: uplo, n, lda, incx, incy, prec
    complex(kind=kind((1d0,0d0))), intent(in)    :: alpha, beta
    complex(kind=kind((1d0,0d0))), intent(in)    :: A(lda, *), x(*)
    complex(kind=kind((1d0,0d0))), intent(inout) :: y(*)
    character :: u
    if (uplo == 121) then
        u = 'U'
    else
        u = 'L'
    end if
    call zhemv(u, n, alpha, A, lda, x, incx, beta, y, incy)
end subroutine

subroutine blas_zhemv2_x(uplo, n, alpha, A, lda, x_head, x_tail, incx, &
                        beta, y, incy, prec)
    integer, intent(in) :: uplo, n, lda, incx, incy, prec
    complex(kind=kind((1d0,0d0))), intent(in)    :: alpha, beta
    complex(kind=kind((1d0,0d0))), intent(in)    :: A(lda, *), x_head(*), x_tail(*)
    complex(kind=kind((1d0,0d0))), intent(inout) :: y(*)
    character :: u
    if (uplo == 121) then
        u = 'U'
    else
        u = 'L'
    end if
    call zhemv(u, n, alpha, A, lda, x_head, incx, beta, y, incy)
    call zhemv(u, n, alpha, A, lda, x_tail, incx, (1.0d0, 0.0d0), y, incy)
end subroutine
