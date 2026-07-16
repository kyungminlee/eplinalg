! dlapll: smallest singular value of [X Y].
program test_dlapll
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_dlapll
    use ref_quad_lapack, only: dlapll
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n
    real(ep), allocatable :: X0(:), Y0(:), X_r(:), Y_r(:), X_g(:), Y_g(:)
    real(ep) :: ssmin_r, ssmin_g, err, tol
    character(len=48) :: label

    call report_init('dlapll', target_name)
    do i = 1, size(ns)
        n = ns(i)
        allocate(X0(n), Y0(n))
        call gen_vector_quad(n, X0, seed = 19401 + 47 * i)
        call gen_vector_quad(n, Y0, seed = 19411 + 47 * i)
        allocate(X_r(n), Y_r(n), X_g(n), Y_g(n))
        X_r = X0; Y_r = Y0; X_g = X0; Y_g = Y0
        call dlapll(n, X_r, 1, Y_r, 1, ssmin_r)
        call target_dlapll(n, X_g, 1, Y_g, 1, ssmin_g)
        err = rel_err_scalar(ssmin_g, ssmin_r)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(X0, Y0, X_r, Y_r, X_g, Y_g)
    end do
    call report_finalize()
end program test_dlapll
