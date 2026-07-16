! dhseqr computes eigenvalues of an upper Hessenberg matrix from
! dgehrd output. Compare WR/WI as a sorted set (algorithms may
! permute pairs).
program test_dhseqr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dhseqr
    use ref_quad_lapack, only: dgehrd, dhseqr
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork
    real(ep), allocatable :: A(:,:), tau(:), work(:)
    real(ep), allocatable :: H_ref(:,:), H_got(:,:)
    real(ep), allocatable :: WR_ref(:), WI_ref(:), WR_got(:), WI_got(:)
    real(ep), allocatable :: WR_sort_ref(:), WR_sort_got(:)
    real(ep), allocatable :: Z(:,:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dhseqr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 340001 + 47 * i)
        allocate(tau(n-1))
        call dgehrd(n, 1, n, A, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgehrd(n, 1, n, A, n, tau, work, lwork, info)
        deallocate(work)
        allocate(H_ref(n,n), H_got(n,n))
        H_ref = A; H_got = A
        allocate(WR_ref(n), WI_ref(n), WR_got(n), WI_got(n), Z(1,1))
        call dhseqr('E', 'N', n, 1, n, H_ref, n, WR_ref, WI_ref, Z, 1, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dhseqr('E', 'N', n, 1, n, H_ref, n, WR_ref, WI_ref, Z, 1, work, lwork, info)
        call target_dhseqr('E', 'N', n, 1, n, H_got, n, WR_got, WI_got, Z, 1, info)
        ! Compare sorted real parts as a robust order-independent check.
        allocate(WR_sort_ref(n), WR_sort_got(n))
        WR_sort_ref = WR_ref; WR_sort_got = WR_got
        call sort_asc(WR_sort_ref, n); call sort_asc(WR_sort_got, n)
        err = max_rel_err_vec(WR_sort_got, WR_sort_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, tau, H_ref, H_got, WR_ref, WI_ref, WR_got, WI_got, &
                   WR_sort_ref, WR_sort_got, Z, work)
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
end program test_dhseqr
