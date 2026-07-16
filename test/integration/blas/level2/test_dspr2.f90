program test_dspr2
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dspr2
    use ref_quad_blas, only: dspr2
    implicit none

    integer, parameter :: cases(*)            = [10, 50, 200]
    character(len=1), parameter :: uplos(*)  = ['U', 'L', 'U']
    integer :: i, n, aps
    real(ep), allocatable :: ap0(:), x(:), y(:), ap_ref(:), ap_got(:)
    real(ep) :: alpha, err, tol
    character(len=32) :: label

    call report_init('dspr2', target_name)
    do i = 1, size(cases)
        n = cases(i)
        aps = n * (n + 1) / 2
        call gen_vector_quad(aps, ap0, seed = 801 + 17 * i)
        call gen_vector_quad(n,    x,  seed = 811 + 17 * i)
        call gen_vector_quad(n,    y,  seed = 821 + 17 * i)
        alpha = 0.7_ep
        allocate(ap_ref(aps), ap_got(aps))
        ap_ref = ap0; ap_got = ap0
        call dspr2(uplos(i), n, alpha, x, 1, y, 1, ap_ref)
        call target_dspr2(uplos(i), n, alpha, x, 1, y, 1, ap_got)
        err = max_rel_err_vec(ap_got, ap_ref)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(ap_ref, ap_got)
    end do
    call report_finalize()
end program test_dspr2
