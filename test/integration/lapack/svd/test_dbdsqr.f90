! Singular values of a bidiagonal matrix via QR. Build (D, E) from a
! random matrix's dgebrd reduction.
program test_dbdsqr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dbdsqr
    use ref_quad_lapack, only: dgebrd, dbdsqr
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, lwork
    real(ep), allocatable :: A(:,:), D0(:), E0(:), tauq(:), taup(:), work(:)
    real(ep), allocatable :: D_ref(:), E_ref(:), D_got(:), E_got(:)
    real(ep), allocatable :: VT(:,:), U(:,:), C(:,:), wbds(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dbdsqr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 290001 + 47 * i)
        allocate(D0(n), E0(n-1), tauq(n), taup(n))
        call dgebrd(n, n, A, n, D0, E0, tauq, taup, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgebrd(n, n, A, n, D0, E0, tauq, taup, work, lwork, info)
        deallocate(work)
        allocate(D_ref(n), E_ref(n-1), D_got(n), E_got(n-1))
        allocate(VT(1,1), U(1,1), C(1,1), wbds(4*n))
        D_ref = D0; E_ref = E0; D_got = D0; E_got = E0
        call dbdsqr('U', n, 0, 0, 0, D_ref, E_ref, VT, 1, U, 1, C, 1, wbds, info)
        call target_dbdsqr('U', n, 0, 0, 0, D_got, E_got, VT, 1, U, 1, C, 1, info)
        err = max_rel_err_vec(D_got, D_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, D0, E0, tauq, taup, D_ref, E_ref, D_got, E_got, VT, U, C, wbds)
    end do
    call report_finalize()
end program test_dbdsqr
