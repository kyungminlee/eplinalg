program test_dtzscal
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_dtzscal
    use test_data,            only: gen_matrix_quad
    use target_ptzblas,       only: target_name, target_eps, target_dtzscal
    implicit none

    integer, parameter :: m = 24, n = 24
    character(len=1), parameter :: uplo_set(*) = ['L','U','L','U','L','U','D']
    integer,          parameter :: ioffd_set(*) = [0, 0, 2, 2, -3, -3, 0]
    integer :: i, ioffd
    character :: uplo
    real(ep), allocatable :: A_got(:,:), A_ref(:,:)
    real(ep) :: alpha, err, tol
    character(len=48) :: label

    call report_init('dtzscal', target_name, 0)
    do i = 1, size(uplo_set)
        uplo = uplo_set(i); ioffd = ioffd_set(i)
        call gen_matrix_quad(m, n, A_got, seed = 1100 + 5 * i)
        A_ref = A_got
        alpha = 0.6_ep + 0.1_ep * real(i, ep)
        call target_dtzscal(uplo, m, n, ioffd, alpha, A_got, m)
        call ref_dtzscal(uplo, m, n, ioffd, alpha, A_ref)
        err = max_rel_err_mat(A_got, A_ref)
        tol = 8.0_ep * target_eps
        write(label, '(a,a,a,i0)') 'uplo=', uplo, ',ioffd=', ioffd
        call report_case(trim(label), err, tol)
        deallocate(A_got, A_ref)
    end do
    call report_finalize()
end program test_dtzscal
