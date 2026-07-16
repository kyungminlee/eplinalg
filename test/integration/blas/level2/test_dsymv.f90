program test_dsymv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_matrix_quad, gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dsymv
    use ref_quad_blas, only: dsymv
    implicit none

    integer, parameter :: cases(*)            = [10, 50, 200]
    character(len=1), parameter :: uplos(*)  = ['U', 'L', 'U']
    integer :: i, n
    real(ep), allocatable :: A(:,:), x(:), y0(:), y_ref(:), y_got(:)
    real(ep) :: alpha, beta, err, tol
    character(len=32) :: label

    call report_init('dsymv', target_name)
    do i = 1, size(cases)
        n = cases(i)
        ! Routine reads only one triangle, so any matrix is fine.
        call gen_matrix_quad(n, n, A, seed = 491 + 17 * i)
        call gen_vector_quad(n, x,  seed = 501 + 17 * i)
        call gen_vector_quad(n, y0, seed = 511 + 17 * i)
        alpha = real(0.6_ep, ep)
        beta  = real(0.4_ep, ep)
        allocate(y_ref(n), y_got(n))
        y_ref = y0
        y_got = y0
        call dsymv(uplos(i), n, alpha, A, n, x, 1, beta, y_ref, 1)
        call target_dsymv(uplos(i), n, alpha, A, n, x, 1, beta, y_got, 1)
        err = max_rel_err_vec(y_got, y_ref)
        tol = 16.0_ep * 2.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(y_ref, y_got)
    end do
    call report_finalize()
end program test_dsymv
