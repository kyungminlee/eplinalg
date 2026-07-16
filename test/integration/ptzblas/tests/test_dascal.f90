program test_dascal
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_vec
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_dascal
    use test_data,            only: gen_vector_quad
    use target_ptzblas,       only: target_name, target_eps, target_dascal
    implicit none

    integer, parameter :: cases(*) = [50, 500, 5000]
    integer :: i, n
    real(ep), allocatable :: x_got(:), x_ref(:)
    real(ep) :: alpha, err, tol
    character(len=32) :: label

    call report_init('dascal', target_name, 0)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_quad(n, x_got, seed = 200 + 5 * i)
        x_ref = x_got
        alpha = -1.3_ep + 0.05_ep * real(i, ep)
        call target_dascal(n, alpha, x_got, 1)
        call ref_dascal(n, alpha, x_ref)
        err = max_rel_err_vec(x_got, x_ref)
        tol = 8.0_ep * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x_got, x_ref)
    end do
    call report_finalize()
end program test_dascal
