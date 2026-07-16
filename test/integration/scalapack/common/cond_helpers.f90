module cond_helpers
    ! Compute exact 1-norm condition numbers kappa_1(A) = ||A||_1 *
    ! ||A^-1||_1 at quad precision via reference LAPACK routines on
    ! rank 0. Used to verify Hager-Higham 1-norm estimators
    ! (PD/PZ-GECON, PD/PZ-POCON, and the rcond field of *POSVX /
    ! *GESVX) which only guarantee kappa_true/3 <= kappa_est
    ! <= kappa_true.
    use prec_kinds,      only: ep
    use ref_quad_lapack, only: dgetrf, dgetri, dpotrf, dpotri, dlange, &
                               dlansy, &
                               zgetrf, zgetri, zpotrf, zpotri, zlange, &
                               zlanhe, zlansy
    implicit none
    private
    public :: true_kappa1_general, true_kappa1_general_z
    public :: true_kappa1_posdef,  true_kappa1_posdef_z

contains

    ! General-square 1-norm condition number from a quad LU + inverse.
    ! A is overwritten; pass a copy.
    subroutine true_kappa1_general(n, A, lda, kappa, info)
        integer,  intent(in)    :: n, lda
        real(ep), intent(inout) :: A(lda,*)
        real(ep), intent(out)   :: kappa
        integer,  intent(out)   :: info
        real(ep), allocatable :: Ainv(:,:), work(:), wn(:)
        integer,  allocatable :: ipiv(:)
        real(ep) :: anorm, ainv_norm
        integer  :: lwork

        kappa = -1.0_ep
        allocate(Ainv(n, n), wn(n), ipiv(n))
        Ainv(1:n, 1:n) = A(1:n, 1:n)

        anorm = dlange('1', n, n, A, lda, wn)

        call dgetrf(n, n, Ainv, n, ipiv, info)
        if (info /= 0) then
            deallocate(Ainv, wn, ipiv)
            return
        end if
        lwork = max(1, 4 * n)
        allocate(work(lwork))
        call dgetri(n, Ainv, n, ipiv, work, lwork, info)
        if (info /= 0) then
            deallocate(Ainv, wn, ipiv, work)
            return
        end if
        ainv_norm = dlange('1', n, n, Ainv, n, wn)
        kappa = anorm * ainv_norm
        deallocate(Ainv, wn, ipiv, work)
    end subroutine true_kappa1_general

    subroutine true_kappa1_general_z(n, A, lda, kappa, info)
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: kappa
        integer,     intent(out)   :: info
        complex(ep), allocatable :: Ainv(:,:), work(:)
        real(ep),    allocatable :: wn(:)
        integer,     allocatable :: ipiv(:)
        real(ep) :: anorm, ainv_norm
        integer  :: lwork

        kappa = -1.0_ep
        allocate(Ainv(n, n), wn(n), ipiv(n))
        Ainv(1:n, 1:n) = A(1:n, 1:n)

        anorm = zlange('1', n, n, A, lda, wn)

        call zgetrf(n, n, Ainv, n, ipiv, info)
        if (info /= 0) then
            deallocate(Ainv, wn, ipiv)
            return
        end if
        lwork = max(1, 4 * n)
        allocate(work(lwork))
        call zgetri(n, Ainv, n, ipiv, work, lwork, info)
        if (info /= 0) then
            deallocate(Ainv, wn, ipiv, work)
            return
        end if
        ainv_norm = zlange('1', n, n, Ainv, n, wn)
        kappa = anorm * ainv_norm
        deallocate(Ainv, wn, ipiv, work)
    end subroutine true_kappa1_general_z

    ! SPD 1-norm condition number via Cholesky + dpotri. UPLO selects
    ! which triangle of A is referenced; the other triangle is filled
    ! by reflection before the symmetric norm is computed.
    subroutine true_kappa1_posdef(uplo, n, A, lda, kappa, info)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, lda
        real(ep),  intent(inout) :: A(lda,*)
        real(ep),  intent(out)   :: kappa
        integer,   intent(out)   :: info
        real(ep), allocatable :: Ainv(:,:), wn(:)
        real(ep) :: anorm, ainv_norm
        integer :: i, j

        kappa = -1.0_ep
        allocate(Ainv(n, n), wn(n))
        Ainv(1:n, 1:n) = A(1:n, 1:n)

        anorm = dlansy('1', uplo, n, A, lda, wn)

        call dpotrf(uplo, n, Ainv, n, info)
        if (info /= 0) then
            deallocate(Ainv, wn)
            return
        end if
        call dpotri(uplo, n, Ainv, n, info)
        if (info /= 0) then
            deallocate(Ainv, wn)
            return
        end if
        ! Mirror so dlange sees the full matrix.
        if (uplo == 'U' .or. uplo == 'u') then
            do j = 1, n
                do i = j + 1, n
                    Ainv(i, j) = Ainv(j, i)
                end do
            end do
        else
            do j = 1, n
                do i = 1, j - 1
                    Ainv(i, j) = Ainv(j, i)
                end do
            end do
        end if
        ainv_norm = dlange('1', n, n, Ainv, n, wn)
        kappa = anorm * ainv_norm
        deallocate(Ainv, wn)
    end subroutine true_kappa1_posdef

    subroutine true_kappa1_posdef_z(uplo, n, A, lda, kappa, info)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, lda
        complex(ep), intent(inout) :: A(lda,*)
        real(ep),    intent(out)   :: kappa
        integer,     intent(out)   :: info
        complex(ep), allocatable :: Ainv(:,:)
        real(ep),    allocatable :: wn(:)
        real(ep) :: anorm, ainv_norm
        integer :: i, j

        kappa = -1.0_ep
        allocate(Ainv(n, n), wn(n))
        Ainv(1:n, 1:n) = A(1:n, 1:n)

        anorm = zlanhe('1', uplo, n, A, lda, wn)

        call zpotrf(uplo, n, Ainv, n, info)
        if (info /= 0) then
            deallocate(Ainv, wn)
            return
        end if
        call zpotri(uplo, n, Ainv, n, info)
        if (info /= 0) then
            deallocate(Ainv, wn)
            return
        end if
        if (uplo == 'U' .or. uplo == 'u') then
            do j = 1, n
                do i = j + 1, n
                    Ainv(i, j) = conjg(Ainv(j, i))
                end do
            end do
        else
            do j = 1, n
                do i = 1, j - 1
                    Ainv(i, j) = conjg(Ainv(j, i))
                end do
            end do
        end if
        ainv_norm = zlange('1', n, n, Ainv, n, wn)
        kappa = anorm * ainv_norm
        deallocate(Ainv, wn)
    end subroutine true_kappa1_posdef_z

end module cond_helpers
