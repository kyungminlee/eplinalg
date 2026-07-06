module compare
    use prec_kinds, only: ep
    implicit none
    private
    public :: rel_err_scalar, rel_err_scalar_z
    public :: max_rel_err_vec, max_rel_err_vec_z
    public :: max_rel_err_mat, max_rel_err_mat_z
    public :: sort_eig_lex_z
contains

    ! Relative error of two scalars: |a-b| / max(|b|, tiny).
    ! Falls back to absolute error when |b| underflows.
    pure function rel_err_scalar(a, b) result(err)
        real(ep), intent(in) :: a, b
        real(ep) :: err, denom
        denom = abs(b)
        if (denom < tiny(1.0_ep)) then
            err = abs(a - b)
        else
            err = abs(a - b) / denom
        end if
    end function rel_err_scalar

    pure function rel_err_scalar_z(a, b) result(err)
        complex(ep), intent(in) :: a, b
        real(ep) :: err, denom
        denom = abs(b)
        if (denom < tiny(1.0_ep)) then
            err = abs(a - b)
        else
            err = abs(a - b) / denom
        end if
    end function rel_err_scalar_z

    ! Max relative error across a vector — same formula, denom is
    ! max|b| so a single near-zero element doesn't blow up the metric.
    pure function max_rel_err_vec(a, b) result(err)
        real(ep), intent(in) :: a(:), b(:)
        real(ep) :: err, denom
        denom = maxval(abs(b))
        if (denom < tiny(1.0_ep)) then
            err = maxval(abs(a - b))
        else
            err = maxval(abs(a - b)) / denom
        end if
    end function max_rel_err_vec

    pure function max_rel_err_vec_z(a, b) result(err)
        complex(ep), intent(in) :: a(:), b(:)
        real(ep) :: err, denom
        denom = maxval(abs(b))
        if (denom < tiny(1.0_ep)) then
            err = maxval(abs(a - b))
        else
            err = maxval(abs(a - b)) / denom
        end if
    end function max_rel_err_vec_z

    pure function max_rel_err_mat(a, b) result(err)
        real(ep), intent(in) :: a(:,:), b(:,:)
        real(ep) :: err, denom
        denom = maxval(abs(b))
        if (denom < tiny(1.0_ep)) then
            err = maxval(abs(a - b))
        else
            err = maxval(abs(a - b)) / denom
        end if
    end function max_rel_err_mat

    pure function max_rel_err_mat_z(a, b) result(err)
        complex(ep), intent(in) :: a(:,:), b(:,:)
        real(ep) :: err, denom
        denom = maxval(abs(b))
        if (denom < tiny(1.0_ep)) then
            err = maxval(abs(a - b))
        else
            err = maxval(abs(a - b)) / denom
        end if
    end function max_rel_err_mat_z

    ! In-place insertion sort on a complex array, ascending lexicographic
    ! by (real, imag). Used to canonicalize eigenvalue spectra before
    ! cross-implementation comparison.
    pure subroutine sort_eig_lex_z(w)
        complex(ep), intent(inout) :: w(:)
        complex(ep) :: key
        integer :: i, j

        do i = 2, size(w)
            key = w(i)
            j = i - 1
            do while (j >= 1)
                if (real(w(j), ep) < real(key, ep)) exit
                if (real(w(j), ep) == real(key, ep) .and. &
                    aimag(w(j)) <= aimag(key)) exit
                w(j + 1) = w(j)
                j = j - 1
            end do
            w(j + 1) = key
        end do
    end subroutine sort_eig_lex_z

end module compare
