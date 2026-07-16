! zlarnv: random number generator (complex).
program test_zlarnv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec_z
    use target_lapack,   only: target_name, target_eps, target_zlarnv
    use ref_quad_lapack, only: zlarnv
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer, parameter :: idists(*) = [1, 2, 3, 4, 5]
    integer :: i, k, n
    integer :: iseed_r(4), iseed_g(4)
    complex(ep), allocatable :: X_r(:), X_g(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zlarnv', target_name)
    do i = 1, size(ns)
        n = ns(i)
        do k = 1, size(idists)
            allocate(X_r(n), X_g(n))
            iseed_r = [13, 27, 41, 53]
            iseed_g = iseed_r
            call zlarnv(idists(k), iseed_r, n, X_r)
            call target_zlarnv(idists(k), iseed_g, n, X_g)
            err = max_rel_err_vec_z(X_g, X_r)
            if (idists(k) >= 3) then
                tol = 32.0_ep * sqrt(target_eps)
            else
                tol = 16.0_ep * target_eps
            end if
            write(label, '(a,i0,a,i0)') 'idist=', idists(k), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(X_r, X_g)
        end do
    end do
    call report_finalize()
end program test_zlarnv
