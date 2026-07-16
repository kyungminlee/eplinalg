program test_dsbgv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad, gen_spd_matrix_quad, &
                                 pack_sym_band_quad
    use target_lapack,   only: target_name, target_eps, target_dsbgv
    use ref_quad_lapack, only: dsbgv
    implicit none

    integer, parameter :: ns(*)  = [16, 32, 48]
    integer, parameter :: kas(*) = [2, 4, 6]
    integer, parameter :: kbs(*) = [1, 2, 3]
    integer :: i, n, ka, kb, ldab, ldbb, info
    real(ep), allocatable :: A(:,:), B(:,:), work(:)
    real(ep), allocatable :: AB_ref(:,:), BB_ref(:,:), AB_got(:,:), BB_got(:,:)
    real(ep), allocatable :: W_ref(:), W_got(:), Z(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dsbgv', target_name)
    do i = 1, size(ns)
        n = ns(i); ka = kas(i); kb = kbs(i)
        ldab = ka + 1; ldbb = kb + 1
        call gen_symmetric_matrix_quad(n, A, seed = 250041 + 47 * i)
        call gen_spd_matrix_quad(n, B, seed = 250051 + 47 * i)
        allocate(AB_ref(ldab, n), BB_ref(ldbb, n), AB_got(ldab, n), BB_got(ldbb, n))
        call pack_sym_band_quad('U', n, ka, A, AB_ref); AB_got = AB_ref
        call pack_sym_band_quad('U', n, kb, B, BB_ref); BB_got = BB_ref
        allocate(W_ref(n), W_got(n), Z(1, 1), work(3*n))
        call dsbgv('N', 'U', n, ka, kb, AB_ref, ldab, BB_ref, ldbb, W_ref, Z, 1, work, info)
        call target_dsbgv('N', 'U', n, ka, kb, AB_got, ldab, BB_got, ldbb, W_got, Z, 1, info)
        err = max_rel_err_vec(W_got, W_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'n=', n, ',ka=', ka, ',kb=', kb
        call report_case(trim(label), err, tol)
        deallocate(A, B, AB_ref, BB_ref, AB_got, BB_got, W_ref, W_got, Z, work)
    end do
    call report_finalize()
end program test_dsbgv
