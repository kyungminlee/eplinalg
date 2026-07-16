! ztpttr: triangular packed -> triangular full (complex).
program test_ztpttr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztpttr
    use ref_quad_lapack, only: ztrttp, ztpttr
    implicit none

    integer, parameter :: ns(*) = [7, 16, 33]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, u, n, info_ref, info_got, nt
    complex(ep), allocatable :: A0(:,:), AP(:), A_ref(:,:), A_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztpttr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        nt = n*(n+1)/2
        call gen_matrix_complex(n, n, A0, seed = 17151 + 83 * i)
        do u = 1, size(uplos)
            allocate(AP(nt), A_ref(n, n), A_got(n, n))
            call ztrttp(uplos(u), n, A0, n, AP, info_ref)
            A_ref = (0.0_ep, 0.0_ep); A_got = (0.0_ep, 0.0_ep)
            call ztpttr(uplos(u), n, AP, A_ref, n, info_ref)
            call target_ztpttr(uplos(u), n, AP, A_got, n, info_got)
            err = max_rel_err_mat_z(A_got, A_ref)
            tol = 16.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(u), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(AP, A_ref, A_got)
        end do
        deallocate(A0)
    end do
    call report_finalize()
end program test_ztpttr
