! zpteqr: eigenvalues of a positive-definite tridiagonal matrix
! (real D/E; Z is complex when eigenvectors are computed).
program test_zpteqr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use target_lapack,   only: target_name, target_eps, target_zpteqr
    use ref_quad_lapack, only: zpteqr
    implicit none

    integer, parameter :: ns(*) = [10, 32, 64]
    integer :: i, n, info, j
    real(ep),    allocatable :: D0(:), E0(:), D_ref(:), E_ref(:), D_got(:), E_got(:), work(:)
    complex(ep), allocatable :: Z(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zpteqr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        allocate(D0(n), E0(max(1, n-1)), D_ref(n), E_ref(max(1, n-1)), &
                 D_got(n), E_got(max(1, n-1)), Z(1, 1), work(max(1, 4*n)))
        do j = 1, n
            D0(j) = 2.0_ep + 0.01_ep * real(j, ep)
        end do
        do j = 1, n-1
            E0(j) = -1.0_ep + 0.001_ep * real(j, ep)
        end do
        D_ref = D0; E_ref = E0; D_got = D0; E_got = E0
        call zpteqr('N', n, D_ref, E_ref, Z, 1, work, info)
        call target_zpteqr('N', n, D_got, E_got, Z, 1, info)
        call sort_desc(D_ref); call sort_desc(D_got)
        err = max_rel_err_vec(D_got, D_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(D0, E0, D_ref, E_ref, D_got, E_got, Z, work)
    end do
    call report_finalize()
contains
    subroutine sort_desc(x)
        real(ep), intent(inout) :: x(:)
        integer :: a, b
        real(ep) :: t
        do a = 1, size(x)-1
            do b = a+1, size(x)
                if (x(b) > x(a)) then
                    t = x(a); x(a) = x(b); x(b) = t
                end if
            end do
        end do
    end subroutine
end program test_zpteqr
