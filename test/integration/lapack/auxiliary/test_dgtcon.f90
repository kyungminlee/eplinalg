program test_dgtcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_dgttrf, target_dgtcon
    use ref_quad_lapack, only: dgttrf, dgtcon
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer :: i, n, info, j
    real(ep), allocatable :: dl0(:), d0(:), du0(:)
    real(ep), allocatable :: dl_ref(:), d_ref(:), du_ref(:), du2_ref(:)
    real(ep), allocatable :: dl_got(:), d_got(:), du_got(:), du2_got(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:), iwork(:)
    real(ep), allocatable :: work(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('dgtcon', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n-1, dl0, seed = 116001 + 47 * i)
        call gen_vector_quad(n,   d0,  seed = 116011 + 47 * i)
        call gen_vector_quad(n-1, du0, seed = 116021 + 47 * i)
        do j = 1, n; d0(j) = d0(j) + real(4, ep); end do
        ! 1-norm: max over col j of |dl(j-1)| + |d(j)| + |du(j)|
        anorm = 0.0_ep
        do j = 1, n
            anorm = max(anorm, &
                merge(abs(dl0(j-1)), 0.0_ep, j > 1) + abs(d0(j)) + &
                merge(abs(du0(j)),   0.0_ep, j < n))
        end do
        allocate(dl_ref(n-1), d_ref(n), du_ref(n-1), du2_ref(n-2))
        allocate(dl_got(n-1), d_got(n), du_got(n-1), du2_got(n-2))
        allocate(ipiv_ref(n), ipiv_got(n), work(2*n), iwork(n))
        dl_ref = dl0; d_ref = d0; du_ref = du0
        dl_got = dl0; d_got = d0; du_got = du0
        call dgttrf(n, dl_ref, d_ref, du_ref, du2_ref, ipiv_ref, info)
        call target_dgttrf(n, dl_got, d_got, du_got, du2_got, ipiv_got, info)
        call dgtcon('1', n, dl_ref, d_ref, du_ref, du2_ref, ipiv_ref, &
                    anorm, rcond_ref, work, iwork, info)
        call target_dgtcon('1', n, dl_got, d_got, du_got, du2_got, ipiv_got, &
                           anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 100.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(dl0, d0, du0, dl_ref, d_ref, du_ref, du2_ref)
        deallocate(dl_got, d_got, du_got, du2_got, ipiv_ref, ipiv_got, work, iwork)
    end do
    call report_finalize()
end program test_dgtcon
