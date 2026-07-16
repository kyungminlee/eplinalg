! dspgst: reduce packed real symmetric gen-eig to standard form.
program test_dspgst
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_spd_matrix_quad, pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, target_dspgst
    use ref_quad_lapack, only: dspgst, dpptrf
    implicit none

    integer, parameter :: ns(*) = [12, 24, 48]
    integer :: i, n, info
    real(ep), allocatable :: A0(:,:), B0(:,:), AP0(:), BP0(:)
    real(ep), allocatable :: AP_r(:), AP_g(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dspgst', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_spd_matrix_quad(n, A0, seed = 21201 + 47 * i)
        call gen_spd_matrix_quad(n, B0, seed = 21211 + 47 * i)
        allocate(AP0(n*(n+1)/2), BP0(n*(n+1)/2))
        call pack_sym_packed_quad('U', n, A0, AP0)
        call pack_sym_packed_quad('U', n, B0, BP0)
        call dpptrf('U', n, BP0, info)
        allocate(AP_r(n*(n+1)/2), AP_g(n*(n+1)/2))
        AP_r = AP0; AP_g = AP0
        call dspgst(1, 'U', n, AP_r, BP0, info)
        call target_dspgst(1, 'U', n, AP_g, BP0, info)
        err = max_rel_err_vec(AP_g, AP_r)
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, AP0, BP0, AP_r, AP_g)
    end do
    call report_finalize()
end program test_dspgst
