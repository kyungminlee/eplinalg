program test_dsytrd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsytrd
    use ref_quad_lapack, only: dsytrd
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: D_ref(:), E_ref(:), tau_ref(:), work(:)
    real(ep), allocatable :: D_got(:), E_got(:), tau_got(:)
    real(ep) :: wopt(1), err, tol
    integer :: lwork
    character(len=48) :: label

    call report_init('dsytrd', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 220001 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n))
        allocate(D_ref(n), E_ref(n-1), tau_ref(n-1))
        allocate(D_got(n), E_got(n-1), tau_got(n-1))
        A_ref = A0; A_got = A0
        call dsytrd('U', n, A_ref, n, D_ref, E_ref, tau_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dsytrd('U', n, A_ref, n, D_ref, E_ref, tau_ref, work, lwork, info)
        call target_dsytrd('U', n, A_got, n, D_got, E_got, tau_got, info)
        err = max(max_rel_err_vec(D_got, D_ref), max_rel_err_vec(E_got, E_ref))
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, D_ref, E_ref, tau_ref, D_got, E_got, tau_got, work)
    end do
    call report_finalize()
end program test_dsytrd
