program test_zhseqr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zhseqr
    use ref_quad_lapack, only: zgehrd, zhseqr
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork
    complex(ep), allocatable :: A(:,:), tau(:), work(:)
    complex(ep), allocatable :: H_ref(:,:), H_got(:,:)
    complex(ep), allocatable :: W_ref(:), W_got(:), Z(:,:)
    real(ep), allocatable :: Wabs_ref(:), Wabs_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label
    integer :: j

    call report_init('zhseqr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 340011 + 47 * i)
        allocate(tau(n-1))
        call zgehrd(n, 1, n, A, n, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgehrd(n, 1, n, A, n, tau, work, lwork, info)
        deallocate(work)
        allocate(H_ref(n,n), H_got(n,n))
        H_ref = A; H_got = A
        allocate(W_ref(n), W_got(n), Z(1,1))
        call zhseqr('E', 'N', n, 1, n, H_ref, n, W_ref, Z, 1, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zhseqr('E', 'N', n, 1, n, H_ref, n, W_ref, Z, 1, work, lwork, info)
        call target_zhseqr('E', 'N', n, 1, n, H_got, n, W_got, Z, 1, info)
        allocate(Wabs_ref(n), Wabs_got(n))
        do j = 1, n; Wabs_ref(j) = abs(W_ref(j)); Wabs_got(j) = abs(W_got(j)); end do
        call sort_asc(Wabs_ref, n); call sort_asc(Wabs_got, n)
        err = max_rel_err_vec(Wabs_got, Wabs_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, tau, H_ref, H_got, W_ref, W_got, Wabs_ref, Wabs_got, Z, work)
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
end program test_zhseqr
