program test_dsytrd_2stage
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsytrd_2stage
    use ref_quad_lapack, only: dsytrd_2stage
    implicit none
    integer, parameter :: ns(*) = [16, 32, 64]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, info, ju, lhous2, lwork, j
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), D_r(:), E_r(:), tau_r(:)
    real(ep), allocatable :: D_g(:), E_g(:), tau_g(:), hous2(:), work(:)
    real(ep) :: hopt(1), wopt(1), err, tol
    character(len=48) :: label
    call report_init('dsytrd_2stage', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 110001 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + real(n, ep); end do
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n), D_r(n), E_r(max(1,n-1)), tau_r(max(1,n-1)))
            allocate(D_g(n), E_g(max(1,n-1)), tau_g(max(1,n-1)))
            A_ref = A0; A_got = A0
            call dsytrd_2stage('N', uplos(ju), n, A_ref, n, D_r, E_r, tau_r, hopt, -1, wopt, -1, info)
            lhous2 = max(1, int(hopt(1))); lwork = max(1, int(wopt(1)))
            allocate(hous2(lhous2), work(lwork))
            call dsytrd_2stage('N', uplos(ju), n, A_ref, n, D_r, E_r, tau_r, hous2, lhous2, work, lwork, info)
            deallocate(hous2, work)
            call target_dsytrd_2stage('N', uplos(ju), n, A_got, n, D_g, E_g, tau_g, info)
            err = max(max_rel_err_vec(D_g, D_r), max_rel_err_vec(E_g, E_r))
            tol = 256.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, D_r, E_r, tau_r, D_g, E_g, tau_g)
        end do
    end do
    call report_finalize()
end program test_dsytrd_2stage
