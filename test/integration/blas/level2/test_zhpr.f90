program test_zhpr
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zhpr
    use ref_quad_blas, only: zhpr
    implicit none

    integer, parameter :: cases(*)            = [10, 50, 200]
    character(len=1), parameter :: uplos(*)  = ['U', 'L', 'U']
    integer :: i, n, aps
    complex(ep), allocatable :: ap0(:), x(:), ap_ref(:), ap_got(:)
    real(ep) :: alpha, err, tol
    character(len=64) :: label

    call report_init('zhpr', target_name)
    do i = 1, size(cases)
        n = cases(i)
        aps = n * (n + 1) / 2
        call gen_vector_complex(aps, ap0, seed = 1801 + 17 * i)
        call gen_vector_complex(n,    x,  seed = 1811 + 17 * i)
        alpha = 0.55_ep
        allocate(ap_ref(aps), ap_got(aps))
        ap_ref = ap0; ap_got = ap0
        call zhpr(uplos(i), n, alpha, x, 1, ap_ref)
        call target_zhpr(uplos(i), n, alpha, x, 1, ap_got)
        err = max_rel_err_vec_z(ap_got, ap_ref)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(ap_ref, ap_got)
    end do
    call report_finalize()
end program test_zhpr
