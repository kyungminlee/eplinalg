program test_drshft
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_drshft
    use test_data,            only: gen_matrix_quad
    use target_ptzblas,       only: target_name, target_eps, target_drshft
    implicit none

    integer, parameter :: m_cases(*)   = [10, 50, 128, 24,  60]
    integer, parameter :: n_cases(*)   = [12, 50,  64, 16,  40]
    integer, parameter :: off_cases(*) = [ 1,  7,  33, -3, -19]
    integer :: i, m, n, off, lda
    real(ep), allocatable :: A_got(:,:), A_ref(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('drshft', target_name, 0)
    do i = 1, size(m_cases)
        m = m_cases(i); n = n_cases(i); off = off_cases(i)
        ! Upstream qrshft requires LDA ≥ M + |offset|.
        lda = m + abs(off)
        call gen_matrix_quad(lda, n, A_got, seed = 700 + 7 * i)
        A_ref = A_got
        call target_drshft(m, n, off, A_got, lda)
        call ref_drshft(m, n, off, A_ref)
        err = max_rel_err_mat(A_got(1:lda, 1:n), A_ref(1:lda, 1:n))
        tol = 4.0_ep * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',off=', off
        call report_case(trim(label), err, tol)
        deallocate(A_got, A_ref)
    end do
    call report_finalize()
end program test_drshft
