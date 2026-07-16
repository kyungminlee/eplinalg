! Minimal interfaces to reflapack_quad — just the dense linear-solve
! routines the MUMPS tests use to compute ground-truth solutions.
!
! reflapack_quad is the vendored Netlib LAPACK compiled with gfortran's
! -freal-8-real-16, so DOUBLE PRECISION / REAL(KIND=8) entities are
! promoted to REAL(KIND=16) without any routine renaming. The
! interfaces below pin the argument types so callers don't accidentally
! drop into implicit-typing-induced KIND=8 pass-by-reference.

module ref_quad_lapack_solve
    use prec_kinds, only: ep
    implicit none

    interface
        subroutine dgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            integer,  intent(in)    :: n, nrhs, lda, ldb
            real(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,  intent(out)   :: ipiv(*), info
        end subroutine dgesv

        subroutine zgesv(n, nrhs, A, lda, ipiv, B, ldb, info)
            import :: ep
            integer,     intent(in)    :: n, nrhs, lda, ldb
            complex(ep), intent(inout) :: A(lda,*), B(ldb,*)
            integer,     intent(out)   :: ipiv(*), info
        end subroutine zgesv
    end interface
end module ref_quad_lapack_solve
