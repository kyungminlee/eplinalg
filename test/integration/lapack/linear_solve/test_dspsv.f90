program test_dspsv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_symmetric_matrix_quad, gen_matrix_quad, &
                                pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, target_dspsv
    use ref_quad_lapack, only: dspsv
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs = 3
    integer :: i, n, info, np
    integer, allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep), allocatable :: A0(:,:), AP_ref(:), AP_got(:), B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dspsv', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_symmetric_matrix_quad(n, A0, seed = 74001 + 47 * i)
        call gen_matrix_quad(n, nrhs, B0, seed = 74011 + 47 * i)
        allocate(AP_ref(np), AP_got(np), ipiv_ref(n), ipiv_got(n))
        allocate(B_ref(n, nrhs), B_got(n, nrhs))
        call pack_sym_packed_quad('U', n, A0, AP_ref); AP_got = AP_ref
        B_ref = B0; B_got = B0
        call dspsv('U', n, nrhs, AP_ref, ipiv_ref, B_ref, n, info)
        call target_dspsv('U', n, nrhs, AP_got, ipiv_got, B_got, n, info)
        err = max_rel_err_mat(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(AP_ref, AP_got, ipiv_ref, ipiv_got, B_ref, B_got)
    end do
    call report_finalize()
end program test_dspsv
