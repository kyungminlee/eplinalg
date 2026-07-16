! dpbtrf: banded Cholesky factorization. For UPLO='U', AB(kd+1+i-j, j)
! holds A(i,j) for max(1, j-kd) <= i <= j. ldab = kd+1.
program test_dpbtrf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_spd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dpbtrf
    use ref_quad_lapack, only: dpbtrf
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kd    = 3
    integer :: i, n, ldab, info, j, k
    real(ep), allocatable :: Adense(:,:), AB0(:,:), AB_ref(:,:), AB_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dpbtrf', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = kd + 1
        call gen_spd_matrix_quad(n, Adense, seed = 63001 + 47 * i)
        ! Drop entries outside the band so the test matrix is genuinely
        ! banded.
        do j = 1, n
            do k = 1, n
                if (abs(j - k) > kd) Adense(k, j) = 0.0_ep
            end do
        end do
        do j = 1, n
            Adense(j, j) = Adense(j, j) + real(2*kd, ep)
        end do
        allocate(AB0(ldab, n))
        AB0 = 0.0_ep
        do j = 1, n
            do k = max(1, j-kd), j
                AB0(kd + 1 + k - j, j) = Adense(k, j)
            end do
        end do
        allocate(AB_ref(ldab, n), AB_got(ldab, n))
        AB_ref = AB0; AB_got = AB0
        call dpbtrf('U', n, kd, AB_ref, ldab, info)
        call target_dpbtrf('U', n, kd, AB_got, ldab, info)
        err = max_rel_err_mat(AB_got, AB_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(Adense, AB0, AB_ref, AB_got)
    end do
    call report_finalize()
end program test_dpbtrf
