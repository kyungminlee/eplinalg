program test_dspgv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad, gen_spd_matrix_quad, &
                                 pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, target_dspgv
    use ref_quad_lapack, only: dspgv
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, np, info
    real(ep), allocatable :: A(:,:), B(:,:), work(:)
    real(ep), allocatable :: AP_ref(:), BP_ref(:), AP_got(:), BP_got(:)
    real(ep), allocatable :: W_ref(:), W_got(:), Z(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dspgv', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_symmetric_matrix_quad(n, A, seed = 250001 + 47 * i)
        call gen_spd_matrix_quad(n, B, seed = 250011 + 47 * i)
        allocate(AP_ref(np), BP_ref(np), AP_got(np), BP_got(np), work(3*n))
        call pack_sym_packed_quad('U', n, A, AP_ref); AP_got = AP_ref
        call pack_sym_packed_quad('U', n, B, BP_ref); BP_got = BP_ref
        allocate(W_ref(n), W_got(n), Z(1, 1))
        call dspgv(1, 'N', 'U', n, AP_ref, BP_ref, W_ref, Z, 1, work, info)
        call target_dspgv(1, 'N', 'U', n, AP_got, BP_got, W_got, Z, 1, info)
        err = max_rel_err_vec(W_got, W_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, work, AP_ref, BP_ref, AP_got, BP_got, W_ref, W_got, Z)
    end do
    call report_finalize()
end program test_dspgv
