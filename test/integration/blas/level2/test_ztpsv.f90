program test_ztpsv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_ztpsv
    use ref_quad_blas, only: ztpsv
    implicit none

    ! Sweep UPLO × TRANS = 6 combos. DIAG='U' is pinned below.
    integer, parameter :: cases(*)              = [10, 50, 200, 64, 100, 50]
    character(len=1), parameter :: uplos(*)    = ['U', 'L', 'U', 'L', 'U', 'L']
    character(len=1), parameter :: transes(*)  = ['N', 'N', 'T', 'T', 'C', 'C']
    integer :: i, n, aps
    complex(ep), allocatable :: ap(:), x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztpsv', target_name)
    do i = 1, size(cases)
        n = cases(i)
        aps = n * (n + 1) / 2
        call gen_vector_complex(aps, ap, seed = 2401 + 17 * i)
        ap = 0.1_ep * ap
        call gen_vector_complex(n,    x0, seed = 2411 + 17 * i)
        allocate(x_ref(n), x_got(n))
        x_ref = x0; x_got = x0
        call ztpsv(uplos(i), transes(i), 'U', n, ap, x_ref, 1)
        call target_ztpsv(uplos(i), transes(i), 'U', n, ap, x_got, 1)
        err = max_rel_err_vec_z(x_got, x_ref)
        tol = 64.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,a,a,i0)') 'uplo=', uplos(i), &
            ',trans=', transes(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(x_ref, x_got)
    end do
    call report_finalize()
end program test_ztpsv
