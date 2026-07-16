program test_dstebz
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dstebz
    use ref_quad_lapack, only: dsytrd, dstebz
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, lwork, m_ref, m_got, ns_ref, ns_got
    real(ep), allocatable :: A(:,:), D(:), E(:), tau(:), work(:)
    real(ep), allocatable :: W_ref(:), W_got(:), wstebz(:)
    integer,  allocatable :: iblock_ref(:), isplit_ref(:), iblock_got(:), isplit_got(:), iworks(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dstebz', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A, seed = 280051 + 47 * i)
        allocate(D(n), E(n-1), tau(n-1))
        call dsytrd('U', n, A, n, D, E, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dsytrd('U', n, A, n, D, E, tau, work, lwork, info)
        deallocate(work)
        allocate(W_ref(n), W_got(n), iblock_ref(n), isplit_ref(n), &
                 iblock_got(n), isplit_got(n), wstebz(4*n), iworks(3*n))
        call dstebz('A', 'E', n, 0.0_ep, 0.0_ep, 0, 0, 0.0_ep, D, E, &
                    m_ref, ns_ref, W_ref, iblock_ref, isplit_ref, wstebz, iworks, info)
        call target_dstebz('A', 'E', n, 0.0_ep, 0.0_ep, 0, 0, 0.0_ep, D, E, &
                           m_got, ns_got, W_got, iblock_got, isplit_got, info)
        err = max_rel_err_vec(W_got(1:m_got), W_ref(1:m_ref))
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, D, E, tau, W_ref, W_got, iblock_ref, isplit_ref, &
                   iblock_got, isplit_got, wstebz, iworks)
    end do
    call report_finalize()
end program test_dstebz
