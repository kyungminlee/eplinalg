program test_ztzpadcpy
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat_z
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_ztzpadcpy
    use test_data,            only: gen_matrix_complex
    use target_ptzblas,       only: target_name, target_eps, target_ztzpadcpy
    implicit none

    integer, parameter :: m_cases(*) = [10, 32, 80]
    integer, parameter :: n_cases(*) = [12, 32, 60]
    character(len=1), parameter :: uplos(*) = ['U', 'L']
    character(len=1), parameter :: diags(*) = ['N', 'U']
    integer :: i, ku, kd, m, n
    complex(ep), allocatable :: A(:,:), B_got(:,:), B_ref(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztzpadcpy', target_name, 0)
    do i = 1, size(m_cases)
        m = m_cases(i); n = n_cases(i)
        do ku = 1, size(uplos)
            do kd = 1, size(diags)
                call gen_matrix_complex(m, n, A,     seed = 1900 + 7*i + 31*ku + 53*kd)
                call gen_matrix_complex(m, n, B_got, seed = 2000 + 11*i + 31*ku + 53*kd)
                B_ref = B_got
                call target_ztzpadcpy(uplos(ku), diags(kd), m, n, 0, A, m, B_got, m)
                call ref_ztzpadcpy(uplos(ku), diags(kd), m, n, 0, A, B_ref)
                err = max_rel_err_mat_z(B_got, B_ref)
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
end program test_ztzpadcpy
