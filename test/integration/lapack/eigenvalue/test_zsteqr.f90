program test_zsteqr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zsteqr
    use ref_quad_lapack, only: zhetrd, zsteqr
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, lwork
    complex(ep), allocatable :: A(:,:), tau(:), work(:), Z(:,:)
    real(ep), allocatable :: D0(:), E0(:), wsteqr(:)
    real(ep), allocatable :: D_ref(:), E_ref(:), D_got(:), E_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zsteqr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A, seed = 280021 + 47 * i)
        allocate(D0(n), E0(n-1), tau(n-1))
        call zhetrd('U', n, A, n, D0, E0, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zhetrd('U', n, A, n, D0, E0, tau, work, lwork, info)
        deallocate(work)
        allocate(D_ref(n), E_ref(n-1), D_got(n), E_got(n-1), Z(1, 1), wsteqr(1))
        D_ref = D0; E_ref = E0; D_got = D0; E_got = E0
        call zsteqr('N', n, D_ref, E_ref, Z, 1, wsteqr, info)
        call target_zsteqr('N', n, D_got, E_got, Z, 1, info)
        err = max_rel_err_vec(D_got, D_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, D0, E0, tau, D_ref, E_ref, D_got, E_got, Z, wsteqr)
    end do
    call report_finalize()
end program test_zsteqr
