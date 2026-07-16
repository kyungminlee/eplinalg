program test_ddisna
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_ddisna
    use ref_quad_lapack, only: ddisna
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info
    real(ep), allocatable :: D(:), Dsorted(:), sep_ref(:), sep_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ddisna', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n, D, seed = 300001 + 47 * i)
        ! ddisna requires sorted D (ascending or descending). Sort ascending.
        allocate(Dsorted(n))
        Dsorted = abs(D)
        call sort_asc(Dsorted, n)
        allocate(sep_ref(n), sep_got(n))
        call ddisna('E', n, n, Dsorted, sep_ref, info)
        call target_ddisna('E', n, n, Dsorted, sep_got, info)
        err = max_rel_err_vec(sep_got, sep_ref)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(D, Dsorted, sep_ref, sep_got)
    end do
    call report_finalize()
contains
    subroutine sort_asc(x, n)
        real(ep), intent(inout) :: x(:)
        integer,  intent(in)    :: n
        integer :: ii, jj
        real(ep) :: tt
        do ii = 1, n - 1
            do jj = ii + 1, n
                if (x(ii) > x(jj)) then
                    tt = x(ii); x(ii) = x(jj); x(jj) = tt
                end if
            end do
        end do
    end subroutine
end program test_ddisna
