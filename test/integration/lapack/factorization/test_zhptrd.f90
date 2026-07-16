program test_zhptrd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zhptrd
    use ref_quad_lapack, only: zhptrd
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, np, info
    complex(ep), allocatable :: A(:,:), AP_ref(:), AP_got(:), tau_ref(:), tau_got(:)
    real(ep), allocatable :: D_ref(:), E_ref(:), D_got(:), E_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zhptrd', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_hermitian_matrix_quad(n, A, seed = 220031 + 47 * i)
        allocate(AP_ref(np), AP_got(np))
        call pack_herm_packed_quad('U', n, A, AP_ref); AP_got = AP_ref
        allocate(D_ref(n), E_ref(n-1), tau_ref(n-1))
        allocate(D_got(n), E_got(n-1), tau_got(n-1))
        call zhptrd('U', n, AP_ref, D_ref, E_ref, tau_ref, info)
        call target_zhptrd('U', n, AP_got, D_got, E_got, tau_got, info)
        err = max(max_rel_err_vec(D_got, D_ref), max_rel_err_vec(E_got, E_ref))
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AP_ref, AP_got, D_ref, E_ref, tau_ref, D_got, E_got, tau_got)
    end do
    call report_finalize()
end program test_zhptrd
