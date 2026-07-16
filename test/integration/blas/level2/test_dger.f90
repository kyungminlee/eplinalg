program test_dger
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad, gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dger
    use ref_quad_blas, only: dger
    implicit none

    integer, parameter :: ms(*) = [5, 50, 200]
    integer, parameter :: ns(*) = [7, 60, 250]
    integer :: i, m, n
    real(ep), allocatable :: A0(:,:), x(:), y(:)
    real(ep), allocatable :: A_ref(:,:), A_got(:,:)
    real(ep) :: alpha, err, tol
    character(len=32) :: label

    call report_init('dger', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 461 + 17 * i)
        call gen_vector_quad(m, x,  seed = 471 + 17 * i)
        call gen_vector_quad(n, y,  seed = 481 + 17 * i)
        alpha = real(0.4_ep, ep)
        allocate(A_ref(m, n), A_got(m, n))
        A_ref = A0
        A_got = A0
        call dger(m, n, alpha, x, 1, y, 1, A_ref, m)
        call target_dger(m, n, alpha, x, 1, y, 1, A_got, m)
        err = max_rel_err_mat(A_got, A_ref)
        tol = 16.0_ep * 4.0_ep * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got)
    end do
    call report_finalize()
end program test_dger
