program test_datrmv
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_vec
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_datrmv
    use test_data,            only: gen_matrix_quad, gen_vector_quad
    use target_ptzblas,       only: target_name, target_eps, target_datrmv
    implicit none

    integer, parameter :: n = 24
    character(len=1), parameter :: uplo_set(*)  = ['L', 'U', 'L', 'U']
    character(len=1), parameter :: trans_set(*) = ['N', 'N', 'T', 'T']
    character(len=1), parameter :: diag_set(*)  = ['N', 'U', 'N', 'U']
    integer :: i
    character :: uplo, trans, diag
    real(ep), allocatable :: A(:,:), x(:), y_got(:), y_ref(:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call report_init('datrmv', target_name, 0)
    do i = 1, size(uplo_set)
        uplo  = uplo_set(i); trans = trans_set(i); diag = diag_set(i)
        call gen_matrix_quad(n, n, A,    seed = 1800 + 3 * i)
        call gen_vector_quad(n, x,       seed = 1900 + 5 * i)
        call gen_vector_quad(n, y_got,   seed = 2000 + 7 * i)
        allocate(y_ref(n)); y_ref = y_got
        alpha = 1.1_ep; beta = -0.2_ep
        call target_datrmv(uplo, trans, diag, n, alpha, A, n, x, 1, beta, y_got, 1)
        call ref_datrmv(uplo, trans, diag, n, alpha, A, x, beta, y_ref)
        err = max_rel_err_vec(y_got, y_ref)
        tol = 32.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,a,a,a)') 'uplo=', uplo, ',trans=', trans, ',diag=', diag
        call report_case(trim(label), err, tol)
        deallocate(A, x, y_got, y_ref)
    end do
    call report_finalize()
end program test_datrmv
