! dgesvj: one-sided Jacobi SVD. Smoke test.
! Historical note: previously failed on multifloats (Phase L5 — tgesvj
! returned Infinity); resolved per tests/lapack/CHANGELOG.md.
program test_dgesvj
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgesvj
    use ref_quad_lapack, only: dgesvj
    implicit none

    integer, parameter :: ms(*) = [16, 24]
    integer, parameter :: ns(*) = [8, 16]
    integer :: i, m, n, info_r, info_g, lwork
    real(ep), allocatable :: A0(:,:), A_r(:,:), A_g(:,:), V_r(:,:), V_g(:,:)
    real(ep), allocatable :: sva_r(:), sva_g(:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgesvj', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 155001 + 47 * i)
        allocate(A_r(m,n), A_g(m,n), V_r(n,n), V_g(n,n))
        allocate(sva_r(n), sva_g(n))
        A_r = A0; A_g = A0
        lwork = max(6, m + n)
        allocate(work(lwork))
        work(1) = 0.0_ep   ! ctol = default
        call dgesvj('G', 'U', 'V', m, n, A_r, m, sva_r, 0, V_r, n, &
                    work, lwork, info_r)
        deallocate(work)
        call target_dgesvj('G', 'U', 'V', m, n, A_g, m, sva_g, 0, V_g, n, &
                           0.0_ep, info_g)
        err = max_rel_err_vec(sva_g, sva_r)
        tol = 256.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_r, A_g, V_r, V_g, sva_r, sva_g)
    end do
    call report_finalize()
end program test_dgesvj
