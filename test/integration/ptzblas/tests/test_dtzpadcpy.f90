program test_dtzpadcpy
    ! Trapezoidal pad-copy. Cycle UPLO and DIAG; IOFFD = 0.
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_dtzpadcpy
    use test_data,            only: gen_matrix_quad
    use target_ptzblas,       only: target_name, target_eps, target_dtzpadcpy
    implicit none

    integer, parameter :: m_cases(*) = [10, 32, 80]
    integer, parameter :: n_cases(*) = [12, 32, 60]
    character(len=1), parameter :: uplos(*) = ['U', 'L']
    character(len=1), parameter :: diags(*) = ['N', 'U']
    integer :: i, ku, kd, m, n
    real(ep), allocatable :: A(:,:), B_got(:,:), B_ref(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtzpadcpy', target_name, 0)
    do i = 1, size(m_cases)
        m = m_cases(i); n = n_cases(i)
        do ku = 1, size(uplos)
            do kd = 1, size(diags)
                call gen_matrix_quad(m, n, A,     seed = 1700 + 7*i + 31*ku + 53*kd)
                call gen_matrix_quad(m, n, B_got, seed = 1800 + 11*i + 31*ku + 53*kd)
                B_ref = B_got
                call target_dtzpadcpy(uplos(ku), diags(kd), m, n, 0, A, m, B_got, m)
                call ref_dtzpadcpy(uplos(ku), diags(kd), m, n, 0, A, B_ref)
                err = max_rel_err_mat(B_got, B_ref)
                tol = 16.0_ep * target_eps
                write(label, '(a,a1,a,a1,a,i0,a,i0)') &
                    'uplo=', uplos(ku), ',diag=', diags(kd), &
                    ',m=', m, ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(A, B_got, B_ref)
            end do
        end do
    end do
    call report_finalize()
end program test_dtzpadcpy
