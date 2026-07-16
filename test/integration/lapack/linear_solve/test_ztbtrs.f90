program test_ztbtrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztbtrs
    use ref_quad_lapack, only: ztbtrs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kd = 3, nrhs = 3
    integer :: i, n, ldab, info, j, k
    complex(ep), allocatable :: A(:,:), AB(:,:), B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztbtrs', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = kd + 1
        call gen_matrix_complex(n, n, A, seed = 105001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + cmplx(real(2*n, ep), 0.0_ep, ep); end do
        allocate(AB(ldab, n)); AB = (0.0_ep, 0.0_ep)
        do j = 1, n
            do k = max(1, j-kd), j
                AB(kd + 1 + k - j, j) = A(k, j)
            end do
        end do
        call gen_matrix_complex(n, nrhs, B0, seed = 105011 + 47 * i)
        allocate(B_ref(n, nrhs), B_got(n, nrhs))
        B_ref = B0; B_got = B0
        call ztbtrs('U', 'N', 'N', n, kd, nrhs, AB, ldab, B_ref, n, info)
        call target_ztbtrs('U', 'N', 'N', n, kd, nrhs, AB, ldab, B_got, n, info)
        err = max_rel_err_mat_z(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AB, B0, B_ref, B_got)
    end do
    call report_finalize()
end program test_ztbtrs
