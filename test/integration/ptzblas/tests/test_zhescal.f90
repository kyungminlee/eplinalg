program test_zhescal
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat_z
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_zhescal
    use test_data,            only: gen_matrix_complex
    use target_ptzblas,       only: target_name, target_eps, target_zhescal
    implicit none

    integer, parameter :: m = 24, n = 24
    ! Cover the three upstream code paths (ALPHA==1 zero-imag-only,
    ! ALPHA==0 delegate-to-tzpad, general scale) across UPLO ∈ {L,U,D}
    ! and IOFFD ∈ {0, +2, -3}.
    character(len=1), parameter :: uplo_set(*) = &
        ['L','U','L','U','L','U', &
         'L','U','D','L','U','D', &
         'L','U','D','L','U','D']
    real(ep), parameter :: alpha_set(*) = &
        [1.5_ep, 1.6_ep, 1.0_ep, 1.0_ep, 0.0_ep, 0.0_ep, &
         1.7_ep, 1.8_ep, 1.9_ep, 1.0_ep, 1.0_ep, 1.0_ep, &
         0.0_ep, 0.0_ep, 0.0_ep, 1.3_ep, 1.4_ep, 1.5_ep]
    integer, parameter :: ioffd_set(*) = &
        [ 0,  0,  0,  0,  0,  0, &
          2,  2,  2,  2,  2,  2, &
         -3, -3, -3, -3, -3, -3]
    integer :: i, ioffd
    character :: uplo
    complex(ep), allocatable :: A_got(:,:), A_ref(:,:)
    real(ep) :: alpha, err, tol
    character(len=64) :: label

    call report_init('zhescal', target_name, 0)
    do i = 1, size(uplo_set)
        uplo = uplo_set(i); alpha = alpha_set(i); ioffd = ioffd_set(i)
        call gen_matrix_complex(m, n, A_got, seed = 2700 + 5 * i)
        A_ref = A_got
        call target_zhescal(uplo, m, n, ioffd, alpha, A_got, m)
        call ref_zhescal(uplo, m, n, ioffd, alpha, A_ref)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 8.0_ep * target_eps
        write(label, '(a,a,a,f4.1,a,i0)') 'uplo=', uplo, ',alpha=', alpha, ',ioffd=', ioffd
        call report_case(trim(label), err, tol)
        deallocate(A_got, A_ref)
    end do
    call report_finalize()
end program test_zhescal
