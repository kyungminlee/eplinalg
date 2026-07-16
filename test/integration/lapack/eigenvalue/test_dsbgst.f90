! dsbgst: reduce real symmetric banded generalized eigenproblem
! A*x = lambda*B*x to standard form. Smoke test.
program test_dsbgst
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_spd_matrix_quad, pack_sym_band_quad
    use target_lapack,   only: target_name, target_eps, target_dsbgst
    use ref_quad_lapack, only: dsbgst, dpbtrf
    implicit none

    integer, parameter :: ns(*)  = [12, 24]
    integer, parameter :: kas(*) = [3, 5]
    integer, parameter :: kbs(*) = [3, 5]
    integer :: i, n, ka, kb, ldab, ldbb, info
    real(ep), allocatable :: A0(:,:), B0(:,:), AB0(:,:), BB0(:,:)
    real(ep), allocatable :: AB_r(:,:), BB_r(:,:), X_r(:,:)
    real(ep), allocatable :: AB_g(:,:), BB_g(:,:), X_g(:,:)
    real(ep), allocatable :: work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dsbgst', target_name)
    do i = 1, size(ns)
        n = ns(i); ka = kas(i); kb = kbs(i)
        ldab = ka + 1; ldbb = kb + 1
        call gen_spd_matrix_quad(n, A0, seed = 21001 + 47 * i)
        call gen_spd_matrix_quad(n, B0, seed = 21011 + 47 * i)
        allocate(AB0(ldab,n), BB0(ldbb,n))
        call pack_sym_band_quad('U', n, ka, A0, AB0)
        call pack_sym_band_quad('U', n, kb, B0, BB0)
        call dpbtrf('U', n, kb, BB0, ldbb, info)
        allocate(AB_r(ldab,n), BB_r(ldbb,n), X_r(n,n))
        allocate(AB_g(ldab,n), BB_g(ldbb,n), X_g(n,n))
        AB_r = AB0; BB_r = BB0; AB_g = AB0; BB_g = BB0
        allocate(work(2*n))
        call dsbgst('N', 'U', n, ka, kb, AB_r, ldab, BB_r, ldbb, X_r, n, work, info)
        call target_dsbgst('N', 'U', n, ka, kb, AB_g, ldab, BB_g, ldbb, X_g, n, info)
        err = max_rel_err_mat(AB_g, AB_r)
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'n=', n, ',ka=', ka, ',kb=', kb
        call report_case(trim(label), err, tol)
        deallocate(work, A0, B0, AB0, BB0, AB_r, BB_r, X_r, AB_g, BB_g, X_g)
    end do
    call report_finalize()
end program test_dsbgst
