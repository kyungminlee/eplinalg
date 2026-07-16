program test_ztpmv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_ztpmv
    use ref_quad_blas, only: ztpmv
    implicit none

    integer, parameter :: cases(*)              = [10, 50, 200]
    character(len=1), parameter :: uplos(*)    = ['U', 'L', 'U']
    character(len=1), parameter :: transes(*)  = ['N', 'C', 'T']
    character(len=1), parameter :: diags(*)    = ['N', 'N', 'U']
    integer :: i, n, aps
    complex(ep), allocatable :: ap(:), x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztpmv', target_name)
    do i = 1, size(cases)
        n = cases(i)
        aps = n * (n + 1) / 2
        call gen_vector_complex(aps, ap, seed = 2301 + 17 * i)
        call gen_vector_complex(n,    x0, seed = 2311 + 17 * i)
        allocate(x_ref(n), x_got(n))
        x_ref = x0; x_got = x0
        call ztpmv(uplos(i), transes(i), diags(i), n, ap, x_ref, 1)
        call target_ztpmv(uplos(i), transes(i), diags(i), n, ap, x_got, 1)
        err = max_rel_err_vec_z(x_got, x_ref)
        tol = 32.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,a,a,a,a,i0)') 'uplo=', uplos(i), &
            ',trans=', transes(i), ',diag=', diags(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(x_ref, x_got)
    end do
    call report_finalize()
end program test_ztpmv
