program test_dsyr2
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad, gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dsyr2
    use ref_quad_blas, only: dsyr2
    implicit none

    integer, parameter :: cases(*)            = [10, 50, 200]
    character(len=1), parameter :: uplos(*)  = ['U', 'L', 'U']
    integer :: i, n
    real(ep), allocatable :: A0(:,:), x(:), y(:), A_ref(:,:), A_got(:,:)
    real(ep) :: alpha, err, tol
    character(len=32) :: label

    call report_init('dsyr2', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_matrix_quad(n, n, A0, seed = 1001 + 17 * i)
        call gen_vector_quad(n,     x, seed = 1011 + 17 * i)
        call gen_vector_quad(n,     y, seed = 1021 + 17 * i)
        alpha = 0.65_ep
        allocate(A_ref(n,n), A_got(n,n))
        A_ref = A0; A_got = A0
        call dsyr2(uplos(i), n, alpha, x, 1, y, 1, A_ref, n)
        call target_dsyr2(uplos(i), n, alpha, x, 1, y, 1, A_got, n)
        err = max_rel_err_mat(A_got, A_ref)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got)
    end do
    call report_finalize()
end program test_dsyr2
