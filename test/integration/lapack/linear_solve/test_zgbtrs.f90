program test_zgbtrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgbtrf, target_zgbtrs
    use ref_quad_lapack, only: zgbtrf, zgbtrs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kl = 2, ku = 3, nrhs = 3
    integer :: i, n, ldab, info, j, k
    complex(ep), allocatable :: Adense(:,:), AB0(:,:), AB_ref(:,:), AB_got(:,:)
    complex(ep), allocatable :: B0(:,:), B_ref(:,:), B_got(:,:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgbtrs', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = 2*kl + ku + 1
        call gen_matrix_complex(n, n, Adense, seed = 93001 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 93011 + 47 * i)
        do j = 1, n
            Adense(j, j) = Adense(j, j) + cmplx(real(2*n, ep), 0.0_ep, ep)
        end do
        allocate(AB0(ldab, n)); AB0 = (0.0_ep, 0.0_ep)
        do j = 1, n
            do k = max(1, j-ku), min(n, j+kl)
                AB0(kl + ku + 1 + k - j, j) = Adense(k, j)
            end do
        end do
        allocate(AB_ref(ldab, n), AB_got(ldab, n), ipiv_ref(n), ipiv_got(n))
        allocate(B_ref(n, nrhs), B_got(n, nrhs))
        AB_ref = AB0; AB_got = AB0
        call zgbtrf(n, n, kl, ku, AB_ref, ldab, ipiv_ref, info)
        call target_zgbtrf(n, n, kl, ku, AB_got, ldab, ipiv_got, info)
        B_ref = B0; B_got = B0
        call zgbtrs('N', n, kl, ku, nrhs, AB_ref, ldab, ipiv_ref, B_ref, n, info)
        call target_zgbtrs('N', n, kl, ku, nrhs, AB_got, ldab, ipiv_got, B_got, n, info)
        err = max_rel_err_mat_z(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(Adense, AB0, AB_ref, AB_got, ipiv_ref, ipiv_got, B_ref, B_got)
    end do
    call report_finalize()
end program test_zgbtrs
