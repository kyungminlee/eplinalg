program test_dcshft
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_dcshft
    use test_data,            only: gen_matrix_quad
    use target_ptzblas,       only: target_name, target_eps, target_dcshft
    implicit none

    integer, parameter :: m_cases(*)   = [10, 50, 128, 24,  60]
    integer, parameter :: n_cases(*)   = [12, 50,  64, 16,  40]
    integer, parameter :: off_cases(*) = [ 1,  5,  17, -2, -11]
    integer :: i, m, n, off, ncols
    real(ep), allocatable :: A_got(:,:), A_ref(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dcshft', target_name, 0)
    do i = 1, size(m_cases)
        m = m_cases(i); n = n_cases(i); off = off_cases(i)
        ! Upstream qcshft writes A(I, J+offset); column dim ≥ n + |offset|.
        ncols = n + abs(off)
        call gen_matrix_quad(m, ncols, A_got, seed = 800 + 7 * i)
        A_ref = A_got
        call target_dcshft(m, n, off, A_got, m)
        call ref_dcshft(m, n, off, A_ref)
        err = max_rel_err_mat(A_got(1:m, 1:ncols), A_ref(1:m, 1:ncols))
        tol = 4.0_ep * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',off=', off
        call report_case(trim(label), err, tol)
        deallocate(A_got, A_ref)
    end do
    call report_finalize()
end program test_dcshft
