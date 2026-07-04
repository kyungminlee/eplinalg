! Quad-precision test-data generators for the MUMPS differential
! precision tests. Produces dense problems plus the matching MUMPS
! triplet (IRN/JCN/A) representation. All outputs are REAL(KIND=ep)
! so they can be passed to either the migrated quad-precision MUMPS
! (qmumps via the target wrapper) or the reflapack_quad reference.
!
! The matrices are small enough for the tests we care about (n ≤ 64)
! that we keep every entry — we don't actually exercise sparsity.

module test_data_mumps
    use prec_kinds, only: ep
    implicit none
    private
    public :: gen_dense_problem
    public :: gen_dense_problem_z
    public :: gen_spd_dense_problem
    public :: gen_hpd_dense_problem
    public :: gen_general_sym_problem
    public :: dense_to_triplet
    public :: dense_to_triplet_z
    public :: dense_to_sym_triplet
    public :: dense_to_sym_triplet_z
    public :: gen_sparse_spd_problem
    public :: gen_sparse_hpd_problem_z

    ! Fortran reproducible-pseudo-random based on a seed array.
    integer, parameter :: rng_seed_size = 33

contains

    subroutine seed_rng(seed)
        integer, intent(in) :: seed
        integer :: i
        integer, allocatable :: seed_arr(:)
        integer :: n
        call random_seed(size=n)
        allocate(seed_arr(n))
        do i = 1, n
            seed_arr(i) = seed * 257 + i * 1009
        end do
        call random_seed(put=seed_arr)
        deallocate(seed_arr)
    end subroutine seed_rng

    ! Generate `A * x_true = b` with general (unsymmetric) A.
    ! A is diagonally-dominant so it's well-conditioned and pivot-stable
    ! for both LU factorization and MUMPS's default partial pivoting.
    subroutine gen_dense_problem(n, A, x_true, b, seed)
        integer,  intent(in)               :: n, seed
        real(ep), allocatable, intent(out) :: A(:,:), x_true(:), b(:)
        real(ep), allocatable              :: r(:,:), rx(:)
        integer :: i

        call seed_rng(seed)
        allocate(A(n, n), x_true(n), b(n), r(n, n), rx(n))
        call random_number(r);  A = 2.0_ep * real(r, ep) - 1.0_ep
        call random_number(rx); x_true = 2.0_ep * real(rx, ep) - 1.0_ep

        ! Diagonal dominance: |A(i,i)| > sum_{j!=i} |A(i,j)|.
        do i = 1, n
            A(i, i) = sign(real(n, ep) + 1.0_ep, A(i, i))
        end do

        b = matmul(A, x_true)
        deallocate(r, rx)
    end subroutine gen_dense_problem

    subroutine gen_dense_problem_z(n, A, x_true, b, seed)
        integer,     intent(in)               :: n, seed
        complex(ep), allocatable, intent(out) :: A(:,:), x_true(:), b(:)
        real(ep),    allocatable              :: rr(:,:), ri(:,:)
        real(ep),    allocatable              :: xr(:), xi(:)
        integer :: i

        call seed_rng(seed)
        allocate(A(n, n), x_true(n), b(n))
        allocate(rr(n, n), ri(n, n), xr(n), xi(n))
        call random_number(rr); call random_number(ri)
        call random_number(xr); call random_number(xi)
        A      = cmplx(2.0_ep * rr - 1.0_ep, 2.0_ep * ri - 1.0_ep, kind=ep)
        x_true = cmplx(2.0_ep * xr - 1.0_ep, 2.0_ep * xi - 1.0_ep, kind=ep)

        do i = 1, n
            A(i, i) = cmplx( &
                sign(real(n, ep) + 1.0_ep, real(A(i, i), ep)), &
                aimag(A(i, i)), kind=ep)
        end do
        b = matmul(A, x_true)
        deallocate(rr, ri, xr, xi)
    end subroutine gen_dense_problem_z

    ! Symmetric positive definite real: A = X*X^T + n*I.
    subroutine gen_spd_dense_problem(n, A, x_true, b, seed)
        integer,  intent(in)               :: n, seed
        real(ep), allocatable, intent(out) :: A(:,:), x_true(:), b(:)
        real(ep), allocatable              :: x(:,:), rx(:)
        integer :: i

        call seed_rng(seed)
        allocate(A(n, n), x_true(n), b(n), x(n, n), rx(n))
        call random_number(x);  x = 2.0_ep * x - 1.0_ep
        call random_number(rx); x_true = 2.0_ep * rx - 1.0_ep

        A = matmul(x, transpose(x))
        do i = 1, n
            A(i, i) = A(i, i) + real(n, ep)
        end do
        b = matmul(A, x_true)
        deallocate(x, rx)
    end subroutine gen_spd_dense_problem

    ! Hermitian positive definite complex: A = X*X^H + n*I.
    subroutine gen_hpd_dense_problem(n, A, x_true, b, seed)
        integer,     intent(in)               :: n, seed
        complex(ep), allocatable, intent(out) :: A(:,:), x_true(:), b(:)
        real(ep),    allocatable              :: rr(:,:), ri(:,:)
        real(ep),    allocatable              :: xr(:), xi(:)
        complex(ep), allocatable              :: x(:,:)
        integer :: i

        call seed_rng(seed)
        allocate(A(n, n), x_true(n), b(n))
        allocate(rr(n, n), ri(n, n), xr(n), xi(n), x(n, n))
        call random_number(rr); call random_number(ri)
        call random_number(xr); call random_number(xi)
        x      = cmplx(2.0_ep * rr - 1.0_ep, 2.0_ep * ri - 1.0_ep, kind=ep)
        x_true = cmplx(2.0_ep * xr - 1.0_ep, 2.0_ep * xi - 1.0_ep, kind=ep)

        A = matmul(x, transpose(conjg(x)))
        do i = 1, n
            A(i, i) = A(i, i) + cmplx(real(n, ep), 0.0_ep, kind=ep)
        end do
        b = matmul(A, x_true)
        deallocate(rr, ri, xr, xi, x)
    end subroutine gen_hpd_dense_problem

    ! General symmetric (not positive definite) real: A = (X + X^T)/2
    ! plus diagonal boost for indefinite-but-stable factorization.
    subroutine gen_general_sym_problem(n, A, x_true, b, seed)
        integer,  intent(in)               :: n, seed
        real(ep), allocatable, intent(out) :: A(:,:), x_true(:), b(:)
        real(ep), allocatable              :: r(:,:), rx(:)
        integer :: i

        call seed_rng(seed)
        allocate(A(n, n), x_true(n), b(n), r(n, n), rx(n))
        call random_number(r);  r = 2.0_ep * r - 1.0_ep
        call random_number(rx); x_true = 2.0_ep * rx - 1.0_ep

        A = 0.5_ep * (r + transpose(r))
        do i = 1, n
            A(i, i) = sign(real(n, ep) + 1.0_ep, A(i, i))
        end do
        b = matmul(A, x_true)
        deallocate(r, rx)
    end subroutine gen_general_sym_problem

    ! Sparse SPD problem in MUMPS triplet form: the 2D 5-point Laplacian
    ! on an nx-by-ny grid (n = nx*ny). Unlike the dense generators, this
    ! yields a genuinely SPARSE adjacency graph with real vertex
    ! separators — the structure that nested-dissection orderings (PORD,
    ! METIS, Scotch; ICNTL(7)=4/5/3) require. Fed a complete/dense graph
    ! those orderings abort ("no valid number of stages in multisector").
    !
    ! Each node couples to its N/S/E/W grid neighbours with -1; the
    ! diagonal is (degree + 1), so A is symmetric, strictly diagonally
    ! dominant and hence SPD and well-conditioned. The full matrix (both
    ! (i,j) and (j,i)) is emitted so it can be sent with SYM=0. x_true is
    ! random and b = A*x_true is formed by sparse mat-vec, so the caller
    ! has an exact reference solution without any dense factorization.
    subroutine gen_sparse_spd_problem(nx, ny, n, x_true, b, &
                                      irn, jcn, A_trip, nz, seed)
        integer,  intent(in)               :: nx, ny, seed
        integer,  intent(out)              :: n, nz
        real(ep), allocatable, intent(out) :: x_true(:), b(:), A_trip(:)
        integer,  allocatable, intent(out) :: irn(:), jcn(:)
        real(ep), allocatable              :: rx(:)
        integer :: ix, iy, node, deg, k, cap
        integer :: nbr(4), nnbr, m

        n = nx * ny
        ! Upper bound on nonzeros: n diagonal + 2 per grid edge.
        cap = n + 2 * ((nx - 1) * ny + nx * (ny - 1))
        allocate(irn(cap), jcn(cap), A_trip(cap))
        allocate(x_true(n), b(n), rx(n))

        call seed_rng(seed)
        call random_number(rx); x_true = 2.0_ep * rx - 1.0_ep

        k = 0
        do iy = 1, ny
            do ix = 1, nx
                node = (iy - 1) * nx + ix
                ! Collect in-grid neighbours (4-connectivity).
                nnbr = 0
                if (ix > 1)  then; nnbr = nnbr + 1; nbr(nnbr) = node - 1;  end if
                if (ix < nx) then; nnbr = nnbr + 1; nbr(nnbr) = node + 1;  end if
                if (iy > 1)  then; nnbr = nnbr + 1; nbr(nnbr) = node - nx; end if
                if (iy < ny) then; nnbr = nnbr + 1; nbr(nnbr) = node + nx; end if
                deg = nnbr
                ! Diagonal: degree + 1  -> strict diagonal dominance.
                k = k + 1
                irn(k) = node; jcn(k) = node
                A_trip(k) = real(deg, ep) + 1.0_ep
                ! Off-diagonals: -1 to each neighbour.
                do m = 1, nnbr
                    k = k + 1
                    irn(k) = node; jcn(k) = nbr(m)
                    A_trip(k) = -1.0_ep
                end do
            end do
        end do
        nz = k

        ! b = A * x_true by sparse mat-vec over the triplets.
        b = 0.0_ep
        do k = 1, nz
            b(irn(k)) = b(irn(k)) + A_trip(k) * x_true(jcn(k))
        end do
        deallocate(rx)
    end subroutine gen_sparse_spd_problem

    ! Complex Hermitian-positive-definite counterpart of
    ! gen_sparse_spd_problem: the same 2D grid graph, but each edge
    ! carries a genuinely complex coupling c on (p,q) and conj(c) on
    ! (q,p) so A is Hermitian; the real (deg+1) diagonal keeps it
    ! strictly diagonally dominant and hence HPD. Emits the full matrix
    ! (SYM=0). x_true is complex-random and b = A*x_true exactly.
    subroutine gen_sparse_hpd_problem_z(nx, ny, n, x_true, b, &
                                        irn, jcn, A_trip, nz, seed)
        integer,     intent(in)               :: nx, ny, seed
        integer,     intent(out)              :: n, nz
        complex(ep), allocatable, intent(out) :: x_true(:), b(:), A_trip(:)
        integer,     allocatable, intent(out) :: irn(:), jcn(:)
        real(ep),    allocatable              :: xr(:), xi(:)
        complex(ep), parameter :: cpl = (-1.0_ep, 0.25_ep)  ! edge coupling
        integer :: ix, iy, node, deg, k, cap
        integer :: nbr(4), nnbr, m

        n = nx * ny
        cap = n + 2 * ((nx - 1) * ny + nx * (ny - 1))
        allocate(irn(cap), jcn(cap), A_trip(cap))
        allocate(x_true(n), b(n), xr(n), xi(n))

        call seed_rng(seed)
        call random_number(xr); call random_number(xi)
        x_true = cmplx(2.0_ep * xr - 1.0_ep, 2.0_ep * xi - 1.0_ep, kind=ep)

        k = 0
        do iy = 1, ny
            do ix = 1, nx
                node = (iy - 1) * nx + ix
                nnbr = 0
                if (ix > 1)  then; nnbr = nnbr + 1; nbr(nnbr) = node - 1;  end if
                if (ix < nx) then; nnbr = nnbr + 1; nbr(nnbr) = node + 1;  end if
                if (iy > 1)  then; nnbr = nnbr + 1; nbr(nnbr) = node - nx; end if
                if (iy < ny) then; nnbr = nnbr + 1; nbr(nnbr) = node + nx; end if
                deg = nnbr
                k = k + 1
                irn(k) = node; jcn(k) = node
                A_trip(k) = cmplx(real(deg, ep) + 1.0_ep, 0.0_ep, kind=ep)
                do m = 1, nnbr
                    k = k + 1
                    irn(k) = node; jcn(k) = nbr(m)
                    ! Hermitian: c on the (lower->higher) edge, conj(c) back.
                    if (node < nbr(m)) then
                        A_trip(k) = cpl
                    else
                        A_trip(k) = conjg(cpl)
                    end if
                end do
            end do
        end do
        nz = k

        b = (0.0_ep, 0.0_ep)
        do k = 1, nz
            b(irn(k)) = b(irn(k)) + A_trip(k) * x_true(jcn(k))
        end do
        deallocate(xr, xi)
    end subroutine gen_sparse_hpd_problem_z

    ! Convert a dense real matrix to MUMPS triplet (IRN/JCN/A_trip).
    ! Every nonzero entry is included; n^2 entries total (small n only).
    subroutine dense_to_triplet(A, irn, jcn, A_trip, nz)
        real(ep), intent(in)               :: A(:,:)
        integer,  allocatable, intent(out) :: irn(:), jcn(:)
        real(ep), allocatable, intent(out) :: A_trip(:)
        integer,  intent(out)              :: nz
        integer :: i, j, n, k

        n = size(A, 1)
        nz = n * n
        allocate(irn(nz), jcn(nz), A_trip(nz))
        k = 0
        do j = 1, n
            do i = 1, n
                k = k + 1
                irn(k)    = i
                jcn(k)    = j
                A_trip(k) = A(i, j)
            end do
        end do
    end subroutine dense_to_triplet

    subroutine dense_to_triplet_z(A, irn, jcn, A_trip, nz)
        complex(ep), intent(in)               :: A(:,:)
        integer,     allocatable, intent(out) :: irn(:), jcn(:)
        complex(ep), allocatable, intent(out) :: A_trip(:)
        integer,     intent(out)              :: nz
        integer :: i, j, n, k

        n = size(A, 1)
        nz = n * n
        allocate(irn(nz), jcn(nz), A_trip(nz))
        k = 0
        do j = 1, n
            do i = 1, n
                k = k + 1
                irn(k)    = i
                jcn(k)    = j
                A_trip(k) = A(i, j)
            end do
        end do
    end subroutine dense_to_triplet_z

    ! Symmetric / Hermitian variants — only the upper triangle is sent
    ! to MUMPS (SYM=1 / SYM=2 contract; the user supplies one half of
    ! the matrix and MUMPS infers the other).
    subroutine dense_to_sym_triplet(A, irn, jcn, A_trip, nz)
        real(ep), intent(in)               :: A(:,:)
        integer,  allocatable, intent(out) :: irn(:), jcn(:)
        real(ep), allocatable, intent(out) :: A_trip(:)
        integer,  intent(out)              :: nz
        integer :: i, j, n, k

        n = size(A, 1)
        nz = n * (n + 1) / 2
        allocate(irn(nz), jcn(nz), A_trip(nz))
        k = 0
        do j = 1, n
            do i = 1, j
                k = k + 1
                irn(k)    = i
                jcn(k)    = j
                A_trip(k) = A(i, j)
            end do
        end do
    end subroutine dense_to_sym_triplet

    subroutine dense_to_sym_triplet_z(A, irn, jcn, A_trip, nz)
        complex(ep), intent(in)               :: A(:,:)
        integer,     allocatable, intent(out) :: irn(:), jcn(:)
        complex(ep), allocatable, intent(out) :: A_trip(:)
        integer,     intent(out)              :: nz
        integer :: i, j, n, k

        n = size(A, 1)
        nz = n * (n + 1) / 2
        allocate(irn(nz), jcn(nz), A_trip(nz))
        k = 0
        do j = 1, n
            do i = 1, j
                k = k + 1
                irn(k)    = i
                jcn(k)    = j
                A_trip(k) = A(i, j)
            end do
        end do
    end subroutine dense_to_sym_triplet_z

end module test_data_mumps
