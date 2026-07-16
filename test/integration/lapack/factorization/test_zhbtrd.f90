program test_zhbtrd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad, pack_herm_band_quad
    use target_lapack,   only: target_name, target_eps, target_zhbtrd
    use ref_quad_lapack, only: zhbtrd
    implicit none

    integer, parameter :: ns(*)  = [16, 32, 64]
    integer, parameter :: kds(*) = [2, 4, 8]
    integer :: i, n, kd, ldab, info
    complex(ep), allocatable :: A(:,:), AB0(:,:), AB_ref(:,:), AB_got(:,:), Q(:,:), work(:)
    real(ep), allocatable :: D_ref(:), E_ref(:), D_got(:), E_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zhbtrd', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_hermitian_matrix_quad(n, A, seed = 220051 + 47 * i)
        allocate(AB0(ldab, n), AB_ref(ldab, n), AB_got(ldab, n), Q(1,1), work(n))
        call pack_herm_band_quad('U', n, kd, A, AB0)
        AB_ref = AB0; AB_got = AB0
        allocate(D_ref(n), E_ref(n-1), D_got(n), E_got(n-1))
        call zhbtrd('N', 'U', n, kd, AB_ref, ldab, D_ref, E_ref, Q, 1, work, info)
        call target_zhbtrd('N', 'U', n, kd, AB_got, ldab, D_got, E_got, Q, 1, info)
        err = max(max_rel_err_vec(D_got, D_ref), max_rel_err_vec(E_got, E_ref))
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'n=', n, ',kd=', kd
        call report_case(trim(label), err, tol)
        deallocate(A, AB0, AB_ref, AB_got, Q, work, D_ref, E_ref, D_got, E_got)
    end do
    call report_finalize()
end program test_zhbtrd
