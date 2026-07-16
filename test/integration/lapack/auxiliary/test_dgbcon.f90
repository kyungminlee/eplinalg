program test_dgbcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgbtrf, target_dgbcon
    use ref_quad_lapack, only: dgbtrf, dgbcon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kl = 2, ku = 3
    integer :: i, n, ldab, info, j, k
    real(ep), allocatable :: Adense(:,:), AB0(:,:), AB_ref(:,:), AB_got(:,:), work(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:), iwork(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('dgbcon', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = 2*kl + ku + 1
        call gen_matrix_quad(n, n, Adense, seed = 114001 + 47 * i)
        do j = 1, n; Adense(j, j) = Adense(j, j) + real(2*n, ep); end do
        allocate(AB0(ldab, n)); AB0 = 0.0_ep
        do j = 1, n
            do k = max(1, j-ku), min(n, j+kl)
                AB0(kl + ku + 1 + k - j, j) = Adense(k, j)
            end do
        end do
        anorm = maxval(sum(abs(AB0(kl+1:ldab, :)), dim=1))
        allocate(AB_ref(ldab, n), AB_got(ldab, n), ipiv_ref(n), ipiv_got(n))
        allocate(work(3*n), iwork(n))
        AB_ref = AB0; AB_got = AB0
        call dgbtrf(n, n, kl, ku, AB_ref, ldab, ipiv_ref, info)
        call target_dgbtrf(n, n, kl, ku, AB_got, ldab, ipiv_got, info)
        call dgbcon('1', n, kl, ku, AB_ref, ldab, ipiv_ref, anorm, rcond_ref, &
                    work, iwork, info)
        call target_dgbcon('1', n, kl, ku, AB_got, ldab, ipiv_got, anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 100.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(Adense, AB0, AB_ref, AB_got, ipiv_ref, ipiv_got, work, iwork)
    end do
    call report_finalize()
end program test_dgbcon
