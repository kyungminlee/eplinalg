! ztrttp: triangular full -> triangular packed (complex).
program test_ztrttp
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztrttp
    use ref_quad_lapack, only: ztrttp
    implicit none

    integer, parameter :: ns(*) = [7, 16, 33]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, u, n, info_ref, info_got, nt
    complex(ep), allocatable :: A(:,:), AP_ref(:), AP_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztrttp', target_name)
    do i = 1, size(ns)
        n = ns(i)
        nt = n*(n+1)/2
        call gen_matrix_complex(n, n, A, seed = 17051 + 73 * i)
        do u = 1, size(uplos)
            allocate(AP_ref(nt), AP_got(nt))
            AP_ref = (0.0_ep, 0.0_ep); AP_got = (0.0_ep, 0.0_ep)
            call ztrttp(uplos(u), n, A, n, AP_ref, info_ref)
            call target_ztrttp(uplos(u), n, A, n, AP_got, info_got)
            err = max_rel_err_vec_z(AP_got, AP_ref)
            tol = 16.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(u), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(AP_ref, AP_got)
        end do
        deallocate(A)
    end do
    call report_finalize()
end program test_ztrttp
