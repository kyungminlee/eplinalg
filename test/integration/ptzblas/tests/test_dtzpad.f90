program test_dtzpad
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_mat
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_dtzpad
    use test_data,            only: gen_matrix_quad
    use target_ptzblas,       only: target_name, target_eps, target_dtzpad
    implicit none

    integer, parameter :: m = 16, n = 16
    ! UPLO ∈ {L,U} × IOFFD ∈ {0, +2, -1} — exercises the upstream two-loop
    ! "full-fill / diagonal+strict" structure that triggers off-by-one
    ! bugs when IOFFD≠0.
    integer,          parameter :: ioffd_set(*) = [ 0,  0,  2,  2, -1, -1,  0]
    character(len=1), parameter :: uplo_set(*)  = ['L','U','L','U','L','U','D']
    character(len=1), parameter :: herm_set(*)  = ['N','N','N','N','N','N','N']
    integer :: i, ioffd
    character :: uplo, herm
    real(ep), allocatable :: A_got(:,:), A_ref(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=64) :: label

    call report_init('dtzpad', target_name, 0)
    do i = 1, size(uplo_set)
        uplo  = uplo_set(i); herm = herm_set(i); ioffd = ioffd_set(i)
        call gen_matrix_quad(m, n, A_got, seed = 900 + 7 * i)
        A_ref = A_got
        alpha = 0.25_ep + 0.1_ep * real(i, ep)
        beta  = -1.0_ep - 0.05_ep * real(i, ep)
        call target_dtzpad(uplo, herm, m, n, ioffd, alpha, beta, A_got, m)
        call ref_dtzpad(uplo, herm, m, n, ioffd, alpha, beta, A_ref)
        err = max_rel_err_mat(A_got, A_ref)
        tol = 4.0_ep * target_eps
        write(label, '(a,a,a,a,a,i0)') 'uplo=', uplo, ',herm=', herm, ',ioffd=', ioffd
        call report_case(trim(label), err, tol)
        deallocate(A_got, A_ref)
    end do
    call report_finalize()
end program test_dtzpad
