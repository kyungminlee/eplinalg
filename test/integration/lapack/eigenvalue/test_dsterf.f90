! Eigenvalues of a symmetric tridiagonal matrix via QR. Build (D, E)
! from a random symmetric matrix's dsytrd reduction so the test data
! has well-conditioned spectrum.
program test_dsterf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsterf
    use ref_quad_lapack, only: dsytrd, dsterf
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, lwork
    real(ep), allocatable :: A(:,:), D0(:), E0(:), tau(:), work(:)
    real(ep), allocatable :: D_ref(:), E_ref(:), D_got(:), E_got(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dsterf', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A, seed = 280001 + 47 * i)
        allocate(D0(n), E0(n-1), tau(n-1))
        call dsytrd('U', n, A, n, D0, E0, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dsytrd('U', n, A, n, D0, E0, tau, work, lwork, info)
        deallocate(work)
        allocate(D_ref(n), E_ref(n-1), D_got(n), E_got(n-1))
        D_ref = D0; E_ref = E0; D_got = D0; E_got = E0
        call dsterf(n, D_ref, E_ref, info)
        call target_dsterf(n, D_got, E_got, info)
        err = max_rel_err_vec(D_got, D_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, D0, E0, tau, D_ref, E_ref, D_got, E_got)
    end do
    call report_finalize()
end program test_dsterf
