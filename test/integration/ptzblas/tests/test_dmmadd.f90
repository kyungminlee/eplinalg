program test_dmmadd
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_dmmadd
    use test_data,            only: gen_matrix_quad
    use target_ptzblas,       only: target_name, target_eps, target_dmmadd
    implicit none

    integer, parameter :: m_cases(*) = [10,  64, 100]
    integer, parameter :: n_cases(*) = [12,  48, 120]
    integer :: i, m, n
    real(ep), allocatable :: A(:,:), B_got(:,:), B_ref(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=32) :: label

    call report_init('dmmadd', target_name, 0)
    do i = 1, size(m_cases)
        m = m_cases(i); n = n_cases(i)
        call gen_matrix_quad(m, n, A,     seed = 300 + 7 * i)
        call gen_matrix_quad(m, n, B_got, seed = 400 + 11 * i)
        B_ref = B_got
        alpha = 1.7_ep; beta = -0.5_ep
        call target_dmmadd(m, n, alpha, A, m, beta, B_got, m)
        call ref_dmmadd(m, n, alpha, A, beta, B_ref)
        err = max_rel_err_mat(B_got, B_ref)
        tol = 16.0_ep * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B_got, B_ref)
    end do
    call report_finalize()
end program test_dmmadd
