program test_zbdsqr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zbdsqr
    use ref_quad_lapack, only: zgebrd, zbdsqr
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, lwork
    complex(ep), allocatable :: A(:,:), tauq(:), taup(:), work(:), VT(:,:), U(:,:), C(:,:)
    real(ep), allocatable :: D0(:), E0(:), rwork(:)
    real(ep), allocatable :: D_ref(:), E_ref(:), D_got(:), E_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zbdsqr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 290011 + 47 * i)
        allocate(D0(n), E0(n-1), tauq(n), taup(n))
        call zgebrd(n, n, A, n, D0, E0, tauq, taup, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgebrd(n, n, A, n, D0, E0, tauq, taup, work, lwork, info)
        deallocate(work)
        allocate(D_ref(n), E_ref(n-1), D_got(n), E_got(n-1))
        allocate(VT(1,1), U(1,1), C(1,1), rwork(4*n))
        D_ref = D0; E_ref = E0; D_got = D0; E_got = E0
        call zbdsqr('U', n, 0, 0, 0, D_ref, E_ref, VT, 1, U, 1, C, 1, rwork, info)
        call target_zbdsqr('U', n, 0, 0, 0, D_got, E_got, VT, 1, U, 1, C, 1, info)
        err = max_rel_err_vec(D_got, D_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, D0, E0, tauq, taup, D_ref, E_ref, D_got, E_got, VT, U, C, rwork)
    end do
    call report_finalize()
end program test_zbdsqr
