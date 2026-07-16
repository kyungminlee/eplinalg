program test_zsytri2x
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_complex_symmetric_quad
    use target_lapack,   only: target_name, target_eps, target_zsytrf, target_zsytri2x
    use ref_quad_lapack, only: zsytrf, zsytri2x
    implicit none
    integer, parameter :: ns(*) = [8, 16, 32]
    integer, parameter :: nb    = 2
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, info, ju, lwork, j
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:), work2(:,:)
    integer,     allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('zsytri2x', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_complex_symmetric_quad(n, A0, seed = 65101 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + cmplx(real(n, ep), 0.0_ep, ep); end do
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
            A_ref = A0; A_got = A0
            call zsytrf(uplos(ju), n, A_ref, n, ipiv_ref, wopt, -1, info)
            lwork = max(1, int(real(wopt(1), ep)))
            allocate(work(lwork))
            call zsytrf(uplos(ju), n, A_ref, n, ipiv_ref, work, lwork, info)
            deallocate(work)
            call target_zsytrf(uplos(ju), n, A_got, n, ipiv_got, info)
            allocate(work2(n+nb+1, nb+3))
            call zsytri2x(uplos(ju), n, A_ref, n, ipiv_ref, work2, nb, info)
            deallocate(work2)
            call target_zsytri2x(uplos(ju), n, A_got, n, ipiv_got, nb, info)
            err = max_rel_err_mat_z(A_got, A_ref)
            tol = 16.0_ep * real(n, ep)**3 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, ipiv_ref, ipiv_got)
        end do
    end do
    call report_finalize()
end program test_zsytri2x
