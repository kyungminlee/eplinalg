program test_zppsvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_hpd_matrix_quad, gen_matrix_complex, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zppsvx
    use ref_quad_lapack, only: zppsvx
    implicit none

    integer, parameter :: ns(*)  = [16, 32, 48]
    integer, parameter :: nrhs   = 2
    integer :: i, n, np, info
    complex(ep), allocatable :: A(:,:), B0(:,:)
    complex(ep), allocatable :: AP_ref(:), AFP_ref(:), B_ref(:,:), X_ref(:,:)
    complex(ep), allocatable :: AP_got(:), AFP_got(:), B_got(:,:), X_got(:,:), work(:)
    real(ep), allocatable :: S_ref(:), S_got(:)
    real(ep), allocatable :: ferr_ref(:), berr_ref(:), ferr_got(:), berr_got(:), rwork(:)
    real(ep) :: rcond_ref, rcond_got, err, tol
    character :: equed_ref, equed_got
    character(len=48) :: label

    call report_init('zppsvx', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_hpd_matrix_quad(n, A, seed = 410061 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 410071 + 47 * i)
        allocate(AP_ref(np), AFP_ref(np), B_ref(n, nrhs), X_ref(n, nrhs), S_ref(n))
        allocate(AP_got(np), AFP_got(np), B_got(n, nrhs), X_got(n, nrhs), S_got(n))
        allocate(ferr_ref(nrhs), berr_ref(nrhs), ferr_got(nrhs), berr_got(nrhs))
        allocate(work(2*n), rwork(n))
        call pack_herm_packed_quad('U', n, A, AP_ref); AP_got = AP_ref
        B_ref = B0; B_got = B0
        equed_ref = 'N'; equed_got = 'N'
        call zppsvx('E', 'U', n, nrhs, AP_ref, AFP_ref, equed_ref, S_ref, &
                    B_ref, n, X_ref, n, rcond_ref, ferr_ref, berr_ref, work, rwork, info)
        call target_zppsvx('E', 'U', n, nrhs, AP_got, AFP_got, equed_got, S_got, &
                           B_got, n, X_got, n, rcond_got, ferr_got, berr_got, info)
        err = max(max_rel_err_mat_z(X_got, X_ref), rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B0, AP_ref, AFP_ref, B_ref, X_ref, S_ref, &
                   AP_got, AFP_got, B_got, X_got, S_got, &
                   ferr_ref, berr_ref, ferr_got, berr_got, work, rwork)
    end do
    call report_finalize()
end program test_zppsvx
