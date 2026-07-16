program test_dsyr
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad, gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dsyr
    use ref_quad_blas, only: dsyr
    implicit none

    integer, parameter :: cases(*)            = [10, 50, 200]
    character(len=1), parameter :: uplos(*)  = ['U', 'L', 'U']
    integer :: i, n
    real(ep), allocatable :: A0(:,:), x(:), A_ref(:,:), A_got(:,:)
    real(ep) :: alpha, err, tol
    character(len=32) :: label

    call report_init('dsyr', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_matrix_quad(n, n, A0, seed = 901 + 17 * i)
        call gen_vector_quad(n,     x, seed = 911 + 17 * i)
        alpha = 0.55_ep
        allocate(A_ref(n,n), A_got(n,n))
        A_ref = A0; A_got = A0
        call dsyr(uplos(i), n, alpha, x, 1, A_ref, n)
        call target_dsyr(uplos(i), n, alpha, x, 1, A_got, n)
        err = max_rel_err_mat(A_got, A_ref)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got)
    end do
    call report_finalize()
end program test_dsyr
