program test_zppsv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hpd_matrix_quad, gen_matrix_complex, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zppsv
    use ref_quad_lapack, only: zppsv
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs = 3
    integer :: i, n, info, np
    complex(ep), allocatable :: A0(:,:), AP_ref(:), AP_got(:), B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zppsv', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_hpd_matrix_quad(n, A0, seed = 90001 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 90011 + 47 * i)
        allocate(AP_ref(np), AP_got(np), B_ref(n, nrhs), B_got(n, nrhs))
        call pack_herm_packed_quad('U', n, A0, AP_ref); AP_got = AP_ref
        B_ref = B0; B_got = B0
        call zppsv('U', n, nrhs, AP_ref, B_ref, n, info)
        call target_zppsv('U', n, nrhs, AP_got, B_got, n, info)
        err = max_rel_err_mat_z(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(AP_ref, AP_got, B_ref, B_got)
    end do
    call report_finalize()
end program test_zppsv
