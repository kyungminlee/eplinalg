program test_dsytrd_sy2sb
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_symmetric_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsytrd_sy2sb
    use ref_quad_lapack, only: dsytrd_sy2sb
    implicit none
    integer, parameter :: ns(*) = [32, 64]
    integer, parameter :: kd    = 4
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, info, ju, lwork, ldab, j
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), AB_r(:,:), AB_g(:,:), tau_r(:), tau_g(:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label
    call report_init('dsytrd_sy2sb', target_name)
    do i = 1, size(ns)
        n = ns(i)
        ldab = kd + 1
        call gen_symmetric_matrix_quad(n, A0, seed = 111001 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + real(n, ep); end do
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n), AB_r(ldab, n), AB_g(ldab, n), tau_r(max(1,n-kd)), tau_g(max(1,n-kd)))
            A_ref = A0; A_got = A0
            call dsytrd_sy2sb(uplos(ju), n, kd, A_ref, n, AB_r, ldab, tau_r, wopt, -1, info)
            lwork = max(1, int(wopt(1)))
            allocate(work(lwork))
            call dsytrd_sy2sb(uplos(ju), n, kd, A_ref, n, AB_r, ldab, tau_r, work, lwork, info)
            deallocate(work)
            call target_dsytrd_sy2sb(uplos(ju), n, kd, A_got, n, AB_g, ldab, tau_g, info)
            err = real(abs(info), ep)
            tol = 0.5_ep
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, AB_r, AB_g, tau_r, tau_g)
        end do
    end do
    call report_finalize()
end program test_dsytrd_sy2sb
