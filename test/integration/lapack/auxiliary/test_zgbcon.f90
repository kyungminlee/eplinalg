program test_zgbcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgbtrf, target_zgbcon
    use ref_quad_lapack, only: zgbtrf, zgbcon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kl = 2, ku = 3
    integer :: i, n, ldab, info, j, k
    complex(ep), allocatable :: Adense(:,:), AB0(:,:), AB_ref(:,:), AB_got(:,:), work(:)
    real(ep),    allocatable :: rwork(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('zgbcon', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = 2*kl + ku + 1
        call gen_matrix_complex(n, n, Adense, seed = 115001 + 47 * i)
        do j = 1, n; Adense(j, j) = Adense(j, j) + cmplx(real(2*n, ep), 0.0_ep, ep); end do
        allocate(AB0(ldab, n)); AB0 = (0.0_ep, 0.0_ep)
        do j = 1, n
            do k = max(1, j-ku), min(n, j+kl)
                AB0(kl + ku + 1 + k - j, j) = Adense(k, j)
            end do
        end do
        anorm = maxval(sum(abs(AB0(kl+1:ldab, :)), dim=1))
        allocate(AB_ref(ldab, n), AB_got(ldab, n), ipiv_ref(n), ipiv_got(n))
        allocate(work(2*n), rwork(n))
        AB_ref = AB0; AB_got = AB0
        call zgbtrf(n, n, kl, ku, AB_ref, ldab, ipiv_ref, info)
        call target_zgbtrf(n, n, kl, ku, AB_got, ldab, ipiv_got, info)
        call zgbcon('1', n, kl, ku, AB_ref, ldab, ipiv_ref, anorm, rcond_ref, &
                    work, rwork, info)
        call target_zgbcon('1', n, kl, ku, AB_got, ldab, ipiv_got, anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 100.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(Adense, AB0, AB_ref, AB_got, ipiv_ref, ipiv_got, work, rwork)
    end do
    call report_finalize()
end program test_zgbcon
