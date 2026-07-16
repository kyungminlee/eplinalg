program test_dsptri
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad, pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, target_dsptrf, target_dsptri
    use ref_quad_lapack, only: dsptrf, dsptri
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info, np
    integer, allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep), allocatable :: A0(:,:), AP_ref(:), AP_got(:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dsptri', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_symmetric_matrix_quad(n, A0, seed = 75001 + 47 * i)
        allocate(AP_ref(np), AP_got(np), ipiv_ref(n), ipiv_got(n), work(n))
        call pack_sym_packed_quad('U', n, A0, AP_ref); AP_got = AP_ref
        call dsptrf('U', n, AP_ref, ipiv_ref, info)
        call target_dsptrf('U', n, AP_got, ipiv_got, info)
        call dsptri('U', n, AP_ref, ipiv_ref, work, info)
        call target_dsptri('U', n, AP_got, ipiv_got, info)
        err = max_rel_err_vec(AP_got, AP_ref)
        tol = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(AP_ref, AP_got, ipiv_ref, ipiv_got, work)
    end do
    call report_finalize()
end program test_dsptri
