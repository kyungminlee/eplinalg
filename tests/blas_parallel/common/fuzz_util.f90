! Shared utilities for parallel-BLAS overlay fuzz/consistency tests.
!
! Provides:
!   - read_seed() / read_cases() : env-driven knobs (BLAS_FUZZ_SEED,
!     BLAS_FUZZ_CASES). Default seed is time-based and printed so
!     failures are reproducible.
!   - seed_rng(seed)             : deterministic Fortran RNG seeding.
!   - rand_int_log(lo, hi)       : log-uniform integer in [lo, hi]
!     for size sweeps that emphasize small AND large.
!   - rand_trans()               : 'N' / 'T' / 'C' uniform.
!   - rand_real10_uniform        : scalar in [-1, 1] at REAL(KIND=10).
!   - fill_matrix_10             : column-major matrix [-1, 1].
!   - max_rel_err_10             : ||A - B||_max / max(|B|_max, tiny)
!     over a column-major slice — same precision both sides (kind10).

module fuzz_util
    implicit none
    private
    public :: read_seed, read_cases, seed_rng
    public :: rand_int_log, rand_trans, rand_trans_complex, rand_alpha_beta
    public :: fill_matrix_10, max_rel_err_10, check_sentinels_10
    public :: rand_alpha_beta_c10, fill_matrix_c10, max_rel_err_c10
    public :: check_sentinels_c10
    public :: rand_pad
    public :: rk10, eps10, sentinel_real, sentinel_imag

    integer, parameter :: rk10 = 10
    real(rk10), parameter :: eps10 = epsilon(1.0_rk10)

    ! Magic values planted in matrix padding rows (i > m). Any kernel
    ! that writes past row m gets caught by check_sentinels_*. Values
    ! are deliberately non-tiny / non-typical so even a wild stray
    ! write is obvious. NaN would be cleaner but trips FP exceptions
    ! on some flag combinations — sentinels stay as concrete numbers.
    real(rk10), parameter :: sentinel_real = -7.3737373737e+18_rk10
    real(rk10), parameter :: sentinel_imag =  4.2424242424e+17_rk10

contains

    function read_seed() result(s)
        integer :: s
        character(len=32) :: buf
        integer :: ios, t(8)
        call get_environment_variable('BLAS_FUZZ_SEED', buf, status=ios)
        if (ios == 0 .and. len_trim(buf) > 0) then
            read(buf, *, iostat=ios) s
            if (ios == 0) return
        end if
        call date_and_time(values=t)
        s = abs(t(7) * 60000 + t(8) * 1000 + t(6) * 1000000 + t(5) * 36000000)
        if (s == 0) s = 1
    end function read_seed

    function read_cases() result(n)
        integer :: n
        character(len=32) :: buf
        integer :: ios
        call get_environment_variable('BLAS_FUZZ_CASES', buf, status=ios)
        if (ios == 0 .and. len_trim(buf) > 0) then
            read(buf, *, iostat=ios) n
            if (ios == 0 .and. n > 0) return
        end if
        n = 200
    end function read_cases

    subroutine seed_rng(s)
        integer, intent(in) :: s
        integer :: n, i
        integer, allocatable :: arr(:)
        call random_seed(size=n)
        allocate(arr(n))
        do i = 1, n
            arr(i) = s + i * 1000003
        end do
        call random_seed(put=arr)
    end subroutine seed_rng

    ! Log-uniform integer in [lo, hi]. Emphasizes small sizes (edge
    ! cases) AND samples large ones — exactly what BLAS fuzz needs.
    ! Handles lo == 0 by shifting: sample in [1, hi-lo+1] then subtract.
    function rand_int_log(lo, hi) result(k)
        integer, intent(in) :: lo, hi
        integer :: k, lo_s, hi_s, off
        real(8) :: u, t
        if (hi <= lo) then
            k = lo; return
        end if
        off = 0
        if (lo < 1) off = 1 - lo
        lo_s = lo + off
        hi_s = hi + off
        call random_number(u)
        t = log(real(lo_s, 8)) + u * (log(real(hi_s, 8)) - log(real(lo_s, 8)))
        k = int(exp(t)) - off
        if (k < lo) k = lo
        if (k > hi) k = hi
    end function rand_int_log

    function rand_trans() result(c)
        character :: c
        real(8) :: u
        call random_number(u)
        if (u < 0.34_8) then
            c = 'N'
        else if (u < 0.67_8) then
            c = 'T'
        else
            c = 'C'  ! same as 'T' for real types — exercises code path
        end if
    end function rand_trans

    ! For complex routines 'C' is distinct (conjugate transpose).
    ! Same distribution as rand_trans — kept separate to make intent
    ! at call sites explicit.
    function rand_trans_complex() result(c)
        character :: c
        c = rand_trans()
    end function rand_trans_complex

    subroutine rand_alpha_beta(alpha, beta)
        real(rk10), intent(out) :: alpha, beta
        real(8) :: u
        ! Mix in {0, 1, -1} corner cases ~30% of the time each.
        call random_number(u)
        if (u < 0.10_8) then
            alpha = 0.0_rk10
        else if (u < 0.20_8) then
            alpha = 1.0_rk10
        else if (u < 0.30_8) then
            alpha = -1.0_rk10
        else
            call random_number(u)
            alpha = real(2.0_8 * u - 1.0_8, rk10)
        end if
        call random_number(u)
        if (u < 0.10_8) then
            beta = 0.0_rk10
        else if (u < 0.20_8) then
            beta = 1.0_rk10
        else if (u < 0.30_8) then
            beta = -1.0_rk10
        else
            call random_number(u)
            beta = real(2.0_8 * u - 1.0_8, rk10)
        end if
    end subroutine rand_alpha_beta

    subroutine fill_matrix_10(A, m, n, lda)
        real(rk10), intent(out) :: A(:)
        integer, intent(in) :: m, n, lda
        real(8) :: u
        integer :: i, j, idx
        do j = 1, n
            do i = 1, m
                idx = (j - 1) * lda + i
                call random_number(u)
                A(idx) = real(2.0_8 * u - 1.0_8, rk10)
            end do
            ! Sentinel pad rows [m+1, lda]. check_sentinels_10 will
            ! flag any kernel that writes here (catches stride-handling
            ! bugs where the kernel uses LDA but iterates to LDA rows).
            do i = m + 1, lda
                idx = (j - 1) * lda + i
                A(idx) = sentinel_real
            end do
        end do
    end subroutine fill_matrix_10

    ! Verify that padding rows [m+1, lda] of a freshly-filled (and
    ! possibly written-to) column-major matrix still hold the
    ! sentinel. Returns 0 if clean, otherwise the column index (1-based)
    ! of the first mismatch.
    function check_sentinels_10(A, m, n, lda) result(bad_col)
        real(rk10), intent(in) :: A(:)
        integer, intent(in) :: m, n, lda
        integer :: bad_col
        integer :: i, j, idx
        bad_col = 0
        if (lda <= m) return
        do j = 1, n
            do i = m + 1, lda
                idx = (j - 1) * lda + i
                if (A(idx) /= sentinel_real) then
                    bad_col = j
                    return
                end if
            end do
        end do
    end function check_sentinels_10

    ! Random padding amount with a fat tail. Mix of:
    !   60% small  : 0..4   (tight strides, common in practice)
    !   30% medium : 5..32  (cache-line / page boundary stress)
    !   10% large  : 33..128 (catches kernels that loop to LDA)
    function rand_pad() result(p)
        integer :: p
        real(8) :: u
        call random_number(u)
        if (u < 0.6_8) then
            p = rand_int_log(0, 4)
        else if (u < 0.9_8) then
            p = rand_int_log(5, 32)
        else
            p = rand_int_log(33, 128)
        end if
    end function rand_pad

    ! Compute max_i,j |Cout(i,j) - Cref(i,j)| / max(|Cref|_max, tiny)
    ! over the active m x n submatrix (skipping padding rows).
    function max_rel_err_10(C_out, C_ref, m, n, ldc) result(err)
        real(rk10), intent(in) :: C_out(:), C_ref(:)
        integer, intent(in) :: m, n, ldc
        real(rk10) :: err
        real(rk10) :: a, ref_max, diff_max
        integer :: i, j, idx
        diff_max = 0.0_rk10
        ref_max  = 0.0_rk10
        do j = 1, n
            do i = 1, m
                idx = (j - 1) * ldc + i
                a = abs(C_ref(idx))
                if (a > ref_max) ref_max = a
                a = abs(C_out(idx) - C_ref(idx))
                if (a > diff_max) diff_max = a
            end do
        end do
        if (ref_max < tiny(1.0_rk10)) then
            err = diff_max
        else
            err = diff_max / ref_max
        end if
    end function max_rel_err_10

    subroutine rand_alpha_beta_c10(alpha, beta)
        complex(rk10), intent(out) :: alpha, beta
        real(8) :: ur, ui, u
        ! Same corner-case shape as the real version: ~30% chance of
        ! {0, 1, -1} for either part, else uniform in unit disc.
        call random_number(u)
        if (u < 0.10_8) then
            alpha = (0.0_rk10, 0.0_rk10)
        else if (u < 0.20_8) then
            alpha = (1.0_rk10, 0.0_rk10)
        else if (u < 0.30_8) then
            alpha = (-1.0_rk10, 0.0_rk10)
        else
            call random_number(ur); call random_number(ui)
            alpha = cmplx(real(2.0_8*ur-1.0_8, rk10), &
                          real(2.0_8*ui-1.0_8, rk10), rk10)
        end if
        call random_number(u)
        if (u < 0.10_8) then
            beta = (0.0_rk10, 0.0_rk10)
        else if (u < 0.20_8) then
            beta = (1.0_rk10, 0.0_rk10)
        else if (u < 0.30_8) then
            beta = (-1.0_rk10, 0.0_rk10)
        else
            call random_number(ur); call random_number(ui)
            beta = cmplx(real(2.0_8*ur-1.0_8, rk10), &
                         real(2.0_8*ui-1.0_8, rk10), rk10)
        end if
    end subroutine rand_alpha_beta_c10

    subroutine fill_matrix_c10(A, m, n, lda)
        complex(rk10), intent(out) :: A(:)
        integer, intent(in) :: m, n, lda
        real(8) :: ur, ui
        integer :: i, j, idx
        do j = 1, n
            do i = 1, m
                idx = (j - 1) * lda + i
                call random_number(ur); call random_number(ui)
                A(idx) = cmplx(real(2.0_8*ur-1.0_8, rk10), &
                               real(2.0_8*ui-1.0_8, rk10), rk10)
            end do
            do i = m + 1, lda
                idx = (j - 1) * lda + i
                A(idx) = cmplx(sentinel_real, sentinel_imag, rk10)
            end do
        end do
    end subroutine fill_matrix_c10

    function check_sentinels_c10(A, m, n, lda) result(bad_col)
        complex(rk10), intent(in) :: A(:)
        integer, intent(in) :: m, n, lda
        integer :: bad_col
        complex(rk10) :: sentinel
        integer :: i, j, idx
        bad_col = 0
        if (lda <= m) return
        sentinel = cmplx(sentinel_real, sentinel_imag, rk10)
        do j = 1, n
            do i = m + 1, lda
                idx = (j - 1) * lda + i
                if (A(idx) /= sentinel) then
                    bad_col = j
                    return
                end if
            end do
        end do
    end function check_sentinels_c10

    function max_rel_err_c10(C_out, C_ref, m, n, ldc) result(err)
        complex(rk10), intent(in) :: C_out(:), C_ref(:)
        integer, intent(in) :: m, n, ldc
        real(rk10) :: err
        real(rk10) :: a, ref_max, diff_max
        integer :: i, j, idx
        diff_max = 0.0_rk10
        ref_max  = 0.0_rk10
        do j = 1, n
            do i = 1, m
                idx = (j - 1) * ldc + i
                a = abs(C_ref(idx))
                if (a > ref_max) ref_max = a
                a = abs(C_out(idx) - C_ref(idx))
                if (a > diff_max) diff_max = a
            end do
        end do
        if (ref_max < tiny(1.0_rk10)) then
            err = diff_max
        else
            err = diff_max / ref_max
        end if
    end function max_rel_err_c10

end module fuzz_util
