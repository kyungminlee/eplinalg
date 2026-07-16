program test_dsytri2x
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_symmetric_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsytrf, target_dsytri2x
    use ref_quad_lapack, only: dsytrf, dsytri2x
    implicit none
    integer, parameter :: ns(*) = [8, 16, 32]
    integer, parameter :: nb    = 2
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, info, ju, lwork, j
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:), work2(:,:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label
    call report_init('dsytri2x', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 25101 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + real(n, ep); end do
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
            A_ref = A0; A_got = A0
            call dsytrf(uplos(ju), n, A_ref, n, ipiv_ref, wopt, -1, info)
            lwork = max(1, int(wopt(1)))
            allocate(work(lwork))
            call dsytrf(uplos(ju), n, A_ref, n, ipiv_ref, work, lwork, info)
            deallocate(work)
            call target_dsytrf(uplos(ju), n, A_got, n, ipiv_got, info)
            allocate(work2(n+nb+1, nb+3))
            call dsytri2x(uplos(ju), n, A_ref, n, ipiv_ref, work2, nb, info)
            deallocate(work2)
            call target_dsytri2x(uplos(ju), n, A_got, n, ipiv_got, nb, info)
            err = max_rel_err_mat(A_got, A_ref)
            tol = 16.0_ep * real(n, ep)**3 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, ipiv_ref, ipiv_got)
        end do
    end do
    call report_finalize()
end program test_dsytri2x
