program test_ztzpad
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat_z
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_ztzpad
    use test_data,            only: gen_matrix_complex
    use target_ptzblas,       only: target_name, target_eps, target_ztzpad
    implicit none

    integer, parameter :: m = 16, n = 16
    ! Cover UPLO ∈ {L,U,D} × HERM ∈ {N,Z} × IOFFD ∈ {0, +2, -1}.  The
    ! UPLO='D' branch is upstream-supported only for ztzpad (real has no
    ! imag to zero); HERM='Z' exercises the zero-imag-of-diagonal path.
    character(len=1), parameter :: uplo_set(*)  = &
        ['L','U','L','U','D','D', &
         'L','U','L','U','D','D', &
         'L','U','L','U','D','D']
    character(len=1), parameter :: herm_set(*)  = &
        ['N','N','Z','Z','N','Z', &
         'N','N','Z','Z','N','Z', &
         'N','N','Z','Z','N','Z']
    integer,          parameter :: ioffd_set(*) = &
        [ 0,  0,  0,  0,  0,  0, &
          2,  2,  2,  2,  2,  2, &
         -1, -1, -1, -1, -1, -1]
    integer :: i, ioffd
    character :: uplo, herm
    complex(ep), allocatable :: A_got(:,:), A_ref(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztzpad', target_name, 0)
    do i = 1, size(uplo_set)
        uplo = uplo_set(i); herm = herm_set(i); ioffd = ioffd_set(i)
        call gen_matrix_complex(m, n, A_got, seed = 3300 + 7 * i)
        A_ref = A_got
        alpha = cmplx(0.25_ep + 0.1_ep * i, 0.3_ep, ep)
        beta  = cmplx(-1.0_ep, 0.5_ep, ep)
        call target_ztzpad(uplo, herm, m, n, ioffd, alpha, beta, A_got, m)
        call ref_ztzpad(uplo, herm, m, n, ioffd, alpha, beta, A_ref)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 4.0_ep * target_eps
        write(label, '(a,a,a,a,a,i0)') 'uplo=', uplo, ',herm=', herm, ',ioffd=', ioffd
        call report_case(trim(label), err, tol)
        deallocate(A_got, A_ref)
    end do
    call report_finalize()
end program test_ztzpad
