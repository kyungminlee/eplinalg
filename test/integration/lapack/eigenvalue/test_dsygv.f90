! dsygv: generalized symmetric eigenvalue problem A*x = lambda*B*x
! with B SPD. Compare eigenvalues only (eigenvectors have B-norm
! sign convention).
program test_dsygv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad, gen_spd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsygv
    use ref_quad_lapack, only: dsygv
    implicit none

    integer, parameter :: ns(*) = [8, 24, 48]
    integer :: i, n, info, lwork
    real(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:), w_ref(:), w_got(:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dsygv', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 34001 + 47 * i)
        call gen_spd_matrix_quad(n, B0, seed = 34011 + 47 * i)
        allocate(A_ref(n,n), B_ref(n,n), A_got(n,n), B_got(n,n), w_ref(n), w_got(n))
        A_ref = A0; A_got = A0; B_ref = B0; B_got = B0
        call dsygv(1, 'N', 'U', n, A_ref, n, B_ref, n, w_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dsygv(1, 'N', 'U', n, A_ref, n, B_ref, n, w_ref, work, lwork, info)
        deallocate(work)
        call target_dsygv(1, 'N', 'U', n, A_got, n, B_got, n, w_got, info)
        err = max_rel_err_vec(w_got, w_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, B_ref, A_got, B_got, w_ref, w_got)
    end do
    call report_finalize()
end program test_dsygv
