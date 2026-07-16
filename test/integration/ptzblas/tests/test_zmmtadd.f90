program test_zmmtadd
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat_z
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_zmmtadd
    use test_data,            only: gen_matrix_complex
    use target_ptzblas,       only: target_name, target_eps, target_zmmtadd
    implicit none

    integer, parameter :: m_cases(*) = [10, 64, 100]
    integer, parameter :: n_cases(*) = [12, 48, 120]
    integer :: i, m, n
    complex(ep), allocatable :: A(:,:), B_got(:,:), B_ref(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('zmmtadd', target_name, 0)
    do i = 1, size(m_cases)
        m = m_cases(i); n = n_cases(i)
        ! A is M×N, B is N×M.
        call gen_matrix_complex(m, n, A,     seed = 3100 + 7 * i)
        call gen_matrix_complex(n, m, B_got, seed = 3200 + 11 * i)
        B_ref = B_got
        alpha = cmplx(0.9_ep, -0.2_ep, ep)
        beta  = cmplx(0.3_ep,  0.4_ep, ep)
        call target_zmmtadd(m, n, alpha, A, m, beta, B_got, n)
        call ref_zmmtadd(m, n, alpha, A, beta, B_ref)
        err = max_rel_err_mat_z(B_got, B_ref)
        tol = 32.0_ep * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B_got, B_ref)
    end do
    call report_finalize()
end program test_zmmtadd
