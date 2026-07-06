module test_data
    use prec_kinds, only: ep
    implicit none
    private
    public :: gen_vector_quad, gen_matrix_quad
    public :: gen_vector_complex, gen_matrix_complex
    public :: gen_spd_matrix_quad, gen_hermitian_matrix_quad
    public :: gen_symmetric_matrix_quad, gen_hpd_matrix_quad
    public :: pack_sym_band_quad, pack_herm_band_quad
    public :: pack_sym_packed_quad, pack_herm_packed_quad
    public :: gen_complex_symmetric_quad
    public :: gen_orthogonal_quad

contains

    ! Seed Fortran's intrinsic RNG with a value derived from `seed` so
    ! every test gets reproducible inputs that differ across cases.
    subroutine seed_rng(seed)
        integer, intent(in) :: seed
        integer :: n, i
        integer, allocatable :: s(:)
        call random_seed(size=n)
        allocate(s(n))
        do i = 1, n
            s(i) = seed + i * 1000003
        end do
        call random_seed(put=s)
    end subroutine seed_rng

    subroutine gen_vector_quad(n, x, seed)
        integer,   intent(in)  :: n, seed
        real(ep),  intent(out), allocatable :: x(:)
        real(8),   allocatable :: r(:)
        integer :: i

        call seed_rng(seed)
        allocate(r(n), x(n))
        call random_number(r)
        do i = 1, n
            x(i) = real(r(i), ep) * 2.0_ep - 1.0_ep
        end do
    end subroutine gen_vector_quad

    subroutine gen_matrix_quad(m, n, A, seed)
        integer,  intent(in)  :: m, n, seed
        real(ep), intent(out), allocatable :: A(:,:)
        real(8),  allocatable :: R(:,:)
        integer :: i, j

        call seed_rng(seed)
        allocate(R(m, n), A(m, n))
        call random_number(R)
        do j = 1, n
            do i = 1, m
                A(i, j) = real(R(i, j), ep) * 2.0_ep - 1.0_ep
            end do
        end do
    end subroutine gen_matrix_quad

    subroutine gen_vector_complex(n, x, seed)
        integer,     intent(in)  :: n, seed
        complex(ep), intent(out), allocatable :: x(:)
        real(ep),    allocatable :: re(:), im(:)

        call gen_vector_quad(n, re, seed)
        call gen_vector_quad(n, im, seed + 1)
        allocate(x(n))
        x = cmplx(re, im, ep)
    end subroutine gen_vector_complex

    subroutine gen_matrix_complex(m, n, A, seed)
        integer,     intent(in)  :: m, n, seed
        complex(ep), intent(out), allocatable :: A(:,:)
        real(ep),    allocatable :: R(:,:), I_(:,:)

        call gen_matrix_quad(m, n, R,  seed)
        call gen_matrix_quad(m, n, I_, seed + 1)
        allocate(A(m, n))
        A = cmplx(R, I_, ep)
    end subroutine gen_matrix_complex

    ! Symmetric positive-definite matrix: A = X*X^T + n*I. The diagonal
    ! shift guarantees SPD even for rank-deficient X. Output is n x n.
    subroutine gen_spd_matrix_quad(n, A, seed)
        integer,  intent(in)  :: n, seed
        real(ep), intent(out), allocatable :: A(:,:)
        real(ep), allocatable :: X(:,:)
        integer :: i

        call gen_matrix_quad(n, n, X, seed)
        allocate(A(n, n))
        A = matmul(X, transpose(X))
        do i = 1, n
            A(i, i) = A(i, i) + real(n, ep)
        end do
    end subroutine gen_spd_matrix_quad

    ! Real symmetric matrix: A = (X + X^T)/2.
    subroutine gen_symmetric_matrix_quad(n, A, seed)
        integer,  intent(in)  :: n, seed
        real(ep), intent(out), allocatable :: A(:,:)
        real(ep), allocatable :: X(:,:)

        call gen_matrix_quad(n, n, X, seed)
        allocate(A(n, n))
        A = 0.5_ep * (X + transpose(X))
    end subroutine gen_symmetric_matrix_quad

    ! Hermitian matrix: A = (X + X^H)/2.
    subroutine gen_hermitian_matrix_quad(n, A, seed)
        integer,     intent(in)  :: n, seed
        complex(ep), intent(out), allocatable :: A(:,:)
        complex(ep), allocatable :: X(:,:)

        call gen_matrix_complex(n, n, X, seed)
        allocate(A(n, n))
        A = 0.5_ep * (X + transpose(conjg(X)))
    end subroutine gen_hermitian_matrix_quad

    ! Hermitian positive-definite: A = X*X^H + n*I.
    subroutine gen_hpd_matrix_quad(n, A, seed)
        integer,     intent(in)  :: n, seed
        complex(ep), intent(out), allocatable :: A(:,:)
        complex(ep), allocatable :: X(:,:)
        integer :: i

        call gen_matrix_complex(n, n, X, seed)
        allocate(A(n, n))
        A = matmul(X, transpose(conjg(X)))
        do i = 1, n
            A(i, i) = cmplx(real(A(i, i), ep) + real(n, ep), 0.0_ep, ep)
        end do
    end subroutine gen_hpd_matrix_quad

    ! Complex symmetric (A = A^T, not Hermitian): A = (X + X^T)/2.
    subroutine gen_complex_symmetric_quad(n, A, seed)
        integer,     intent(in)  :: n, seed
        complex(ep), intent(out), allocatable :: A(:,:)
        complex(ep), allocatable :: X(:,:)

        call gen_matrix_complex(n, n, X, seed)
        allocate(A(n, n))
        A = 0.5_ep * (X + transpose(X))
    end subroutine gen_complex_symmetric_quad

    ! Pack a symmetric matrix into LAPACK banded storage AB(kd+1, n).
    ! Off-band entries of the source are first zeroed (so the dense
    ! reference and the banded reference describe the same matrix).
    subroutine pack_sym_band_quad(uplo, n, kd, A, AB)
        character, intent(in)    :: uplo
        integer,   intent(in)    :: n, kd
        real(ep),  intent(inout) :: A(n, n)
        real(ep),  intent(out)   :: AB(kd+1, n)
        integer :: i, j

        do j = 1, n
            do i = 1, n
                if (abs(i - j) > kd) A(i, j) = 0.0_ep
            end do
        end do
        AB = 0.0_ep
        if (uplo == 'U' .or. uplo == 'u') then
            do j = 1, n
                do i = max(1, j - kd), j
                    AB(kd + 1 + i - j, j) = A(i, j)
                end do
            end do
        else
            do j = 1, n
                do i = j, min(n, j + kd)
                    AB(1 + i - j, j) = A(i, j)
                end do
            end do
        end if
    end subroutine pack_sym_band_quad

    subroutine pack_herm_band_quad(uplo, n, kd, A, AB)
        character,   intent(in)    :: uplo
        integer,     intent(in)    :: n, kd
        complex(ep), intent(inout) :: A(n, n)
        complex(ep), intent(out)   :: AB(kd+1, n)
        integer :: i, j

        do j = 1, n
            do i = 1, n
                if (abs(i - j) > kd) A(i, j) = (0.0_ep, 0.0_ep)
            end do
        end do
        AB = (0.0_ep, 0.0_ep)
        if (uplo == 'U' .or. uplo == 'u') then
            do j = 1, n
                do i = max(1, j - kd), j
                    AB(kd + 1 + i - j, j) = A(i, j)
                end do
            end do
        else
            do j = 1, n
                do i = j, min(n, j + kd)
                    AB(1 + i - j, j) = A(i, j)
                end do
            end do
        end if
    end subroutine pack_herm_band_quad

    ! Pack a symmetric matrix into LAPACK packed storage AP(n*(n+1)/2).
    subroutine pack_sym_packed_quad(uplo, n, A, AP)
        character, intent(in)  :: uplo
        integer,   intent(in)  :: n
        real(ep),  intent(in)  :: A(n, n)
        real(ep),  intent(out) :: AP(*)
        integer :: i, j, k

        k = 0
        if (uplo == 'U' .or. uplo == 'u') then
            do j = 1, n
                do i = 1, j
                    k = k + 1
                    AP(k) = A(i, j)
                end do
            end do
        else
            do j = 1, n
                do i = j, n
                    k = k + 1
                    AP(k) = A(i, j)
                end do
            end do
        end if
    end subroutine pack_sym_packed_quad

    subroutine pack_herm_packed_quad(uplo, n, A, AP)
        character,   intent(in)  :: uplo
        integer,     intent(in)  :: n
        complex(ep), intent(in)  :: A(n, n)
        complex(ep), intent(out) :: AP(*)
        integer :: i, j, k

        k = 0
        if (uplo == 'U' .or. uplo == 'u') then
            do j = 1, n
                do i = 1, j
                    k = k + 1
                    AP(k) = A(i, j)
                end do
            end do
        else
            do j = 1, n
                do i = j, n
                    k = k + 1
                    AP(k) = A(i, j)
                end do
            end do
        end if
    end subroutine pack_herm_packed_quad

    ! Random orthogonal n x n matrix via modified Gram-Schmidt on a
    ! uniform random seed matrix. Sign-fixed so the diagonal of the
    ! implicit upper-triangular factor R is non-negative — eliminates
    ! seed-to-seed sign drift in U columns.
    subroutine gen_orthogonal_quad(n, U, seed)
        integer,  intent(in)  :: n, seed
        real(ep), intent(out), allocatable :: U(:,:)
        real(ep), allocatable :: Z(:,:)
        real(ep) :: nrm, dotp, sgn
        integer :: j, k

        call gen_matrix_quad(n, n, Z, seed)
        allocate(U(n, n))
        do j = 1, n
            U(:, j) = Z(:, j)
            do k = 1, j - 1
                dotp = dot_product(U(:, k), U(:, j))
                U(:, j) = U(:, j) - dotp * U(:, k)
            end do
            nrm = sqrt(dot_product(U(:, j), U(:, j)))
            if (nrm > tiny(1.0_ep)) then
                U(:, j) = U(:, j) / nrm
            end if
            sgn = 1.0_ep
            if (U(maxloc(abs(U(:, j)), dim=1), j) < 0.0_ep) sgn = -1.0_ep
            U(:, j) = sgn * U(:, j)
        end do
    end subroutine gen_orthogonal_quad

end module test_data
