! dtpttr: triangular packed -> triangular full.
program test_dtpttr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtpttr
    use ref_quad_lapack, only: dtrttp, dtpttr
    implicit none

    integer, parameter :: ns(*) = [7, 16, 33]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, u, n, info_ref, info_got, nt
    real(ep), allocatable :: A0(:,:), AP(:), A_ref(:,:), A_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtpttr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        nt = n*(n+1)/2
        call gen_matrix_quad(n, n, A0, seed = 17101 + 79 * i)
        do u = 1, size(uplos)
            allocate(AP(nt), A_ref(n, n), A_got(n, n))
            call dtrttp(uplos(u), n, A0, n, AP, info_ref)
            A_ref = 0.0_ep; A_got = 0.0_ep
            call dtpttr(uplos(u), n, AP, A_ref, n, info_ref)
            call target_dtpttr(uplos(u), n, AP, A_got, n, info_got)
            err = max_rel_err_mat(A_got, A_ref)
            tol = 16.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(u), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(AP, A_ref, A_got)
        end do
        deallocate(A0)
    end do
    call report_finalize()
end program test_dtpttr
