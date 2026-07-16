program test_dtptri
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad, pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, target_dtptri
    use ref_quad_lapack, only: dtptri
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, np, info, j
    real(ep), allocatable :: A(:,:), AP_ref(:), AP_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtptri', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_matrix_quad(n, n, A, seed = 108001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + real(2*n, ep); end do
        allocate(AP_ref(np), AP_got(np))
        call pack_sym_packed_quad('U', n, A, AP_ref); AP_got = AP_ref
        call dtptri('U', 'N', n, AP_ref, info)
        call target_dtptri('U', 'N', n, AP_got, info)
        err = max_rel_err_vec(AP_got, AP_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AP_ref, AP_got)
    end do
    call report_finalize()
end program test_dtptri
