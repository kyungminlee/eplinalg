program test_zpptri
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec_z
    use test_data,       only: gen_hpd_matrix_quad, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zpptrf, target_zpptri
    use ref_quad_lapack, only: zpptrf, zpptri
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info, np
    complex(ep), allocatable :: A0(:,:), AP_ref(:), AP_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zpptri', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_hpd_matrix_quad(n, A0, seed = 91001 + 47 * i)
        allocate(AP_ref(np), AP_got(np))
        call pack_herm_packed_quad('U', n, A0, AP_ref); AP_got = AP_ref
        call zpptrf('U', n, AP_ref, info)
        call target_zpptrf('U', n, AP_got, info)
        call zpptri('U', n, AP_ref, info)
        call target_zpptri('U', n, AP_got, info)
        err = max_rel_err_vec_z(AP_got, AP_ref)
        tol = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(AP_ref, AP_got)
    end do
    call report_finalize()
end program test_zpptri
