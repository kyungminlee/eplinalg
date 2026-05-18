program test_dtpsv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dtpsv
    use ref_quad_blas, only: dtpsv
    implicit none

    ! Sweep UPLO × TRANS = 4 combos. DIAG='U' is pinned below.
    integer, parameter :: cases(*)              = [10, 50, 200, 64]
    character(len=1), parameter :: uplos(*)    = ['U', 'L', 'U', 'L']
    character(len=1), parameter :: transes(*)  = ['N', 'N', 'T', 'T']
    integer :: i, n, aps
    real(ep), allocatable :: ap(:), x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtpsv', target_name)
    do i = 1, size(cases)
        n = cases(i)
        aps = n * (n + 1) / 2
        call gen_vector_quad(aps, ap, seed = 1201 + 17 * i)
        ap = 0.1_ep * ap
        call gen_vector_quad(n,    x0, seed = 1211 + 17 * i)
        allocate(x_ref(n), x_got(n))
        x_ref = x0; x_got = x0
        call dtpsv(uplos(i), transes(i), 'U', n, ap, x_ref, 1)
        call target_dtpsv(uplos(i), transes(i), 'U', n, ap, x_got, 1)
        err = max_rel_err_vec(x_got, x_ref)
        tol = 32.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,a,a,i0)') 'uplo=', uplos(i), &
            ',trans=', transes(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(x_ref, x_got)
    end do
    call report_finalize()
end program test_dtpsv
