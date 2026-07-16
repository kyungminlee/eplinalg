program test_ztptri
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec_z
    use test_data,       only: gen_matrix_complex, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_ztptri
    use ref_quad_lapack, only: ztptri
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, np, info, j
    complex(ep), allocatable :: A(:,:), AP_ref(:), AP_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztptri', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_matrix_complex(n, n, A, seed = 109001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + cmplx(real(2*n, ep), 0.0_ep, ep); end do
        allocate(AP_ref(np), AP_got(np))
        call pack_herm_packed_quad('U', n, A, AP_ref); AP_got = AP_ref
        call ztptri('U', 'N', n, AP_ref, info)
        call target_ztptri('U', 'N', n, AP_got, info)
        err = max_rel_err_vec_z(AP_got, AP_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AP_ref, AP_got)
    end do
    call report_finalize()
end program test_ztptri
