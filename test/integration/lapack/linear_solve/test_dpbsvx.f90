program test_dpbsvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, rel_err_scalar
    use test_data,       only: gen_spd_matrix_quad, gen_matrix_quad, pack_sym_band_quad
    use target_lapack,   only: target_name, target_eps, target_dpbsvx
    use ref_quad_lapack, only: dpbsvx
    implicit none

    integer, parameter :: ns(*)  = [16, 32, 48]
    integer, parameter :: kds(*) = [2, 4, 6]
    integer, parameter :: nrhs    = 2
    integer :: i, n, kd, ldab, info
    real(ep), allocatable :: A(:,:), B0(:,:)
    real(ep), allocatable :: AB_ref(:,:), AFB_ref(:,:), B_ref(:,:), X_ref(:,:), S_ref(:)
    real(ep), allocatable :: AB_got(:,:), AFB_got(:,:), B_got(:,:), X_got(:,:), S_got(:)
    real(ep), allocatable :: ferr_ref(:), berr_ref(:), ferr_got(:), berr_got(:), work(:)
    integer,  allocatable :: iwork(:)
    real(ep) :: rcond_ref, rcond_got, err, tol
    character :: equed_ref, equed_got
    character(len=48) :: label

    call report_init('dpbsvx', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_spd_matrix_quad(n, A, seed = 410001 + 47 * i)
        call gen_matrix_quad(n, nrhs, B0, seed = 410011 + 47 * i)
        allocate(AB_ref(ldab, n), AFB_ref(ldab, n), B_ref(n, nrhs), X_ref(n, nrhs), S_ref(n))
        allocate(AB_got(ldab, n), AFB_got(ldab, n), B_got(n, nrhs), X_got(n, nrhs), S_got(n))
        allocate(ferr_ref(nrhs), berr_ref(nrhs), ferr_got(nrhs), berr_got(nrhs))
        allocate(work(3*n), iwork(n))
        call pack_sym_band_quad('U', n, kd, A, AB_ref); AB_got = AB_ref
        B_ref = B0; B_got = B0
        equed_ref = 'N'; equed_got = 'N'
        call dpbsvx('E', 'U', n, kd, nrhs, AB_ref, ldab, AFB_ref, ldab, equed_ref, &
                    S_ref, B_ref, n, X_ref, n, rcond_ref, ferr_ref, berr_ref, work, iwork, info)
        call target_dpbsvx('E', 'U', n, kd, nrhs, AB_got, ldab, AFB_got, ldab, equed_got, &
                           S_got, B_got, n, X_got, n, rcond_got, ferr_got, berr_got, info)
        err = max(max_rel_err_mat(X_got, X_ref), rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'n=', n, ',kd=', kd
        call report_case(trim(label), err, tol)
        deallocate(A, B0, AB_ref, AFB_ref, B_ref, X_ref, S_ref, &
                   AB_got, AFB_got, B_got, X_got, S_got, &
                   ferr_ref, berr_ref, ferr_got, berr_got, work, iwork)
    end do
    call report_finalize()
end program test_dpbsvx
