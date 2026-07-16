program test_ztzcnjg
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat_z
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_ztzcnjg
    use test_data,            only: gen_matrix_complex
    use target_ptzblas,       only: target_name, target_eps, target_ztzcnjg
    implicit none

    integer, parameter :: m = 20, n = 20
    character(len=1), parameter :: uplo_set(*)  = ['L','U','L','U','L','U','D']
    integer,          parameter :: ioffd_set(*) = [0, 0, 2, 2, -3, -3, 0]
    integer :: i, ioffd
    character :: uplo
    complex(ep), allocatable :: A_got(:,:), A_ref(:,:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztzcnjg', target_name, 0)
    do i = 1, size(uplo_set)
        uplo = uplo_set(i); ioffd = ioffd_set(i)
        call gen_matrix_complex(m, n, A_got, seed = 2800 + 5 * i)
        A_ref = A_got
        alpha = cmplx(0.7_ep + 0.1_ep * real(i, ep), 0.3_ep, ep)
        call target_ztzcnjg(uplo, m, n, ioffd, alpha, A_got, m)
        call ref_ztzcnjg(uplo, m, n, ioffd, alpha, A_ref)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 8.0_ep * target_eps
        write(label, '(a,a,a,i0)') 'uplo=', uplo, ',ioffd=', ioffd
        call report_case(trim(label), err, tol)
        deallocate(A_got, A_ref)
    end do
    call report_finalize()
end program test_ztzcnjg
