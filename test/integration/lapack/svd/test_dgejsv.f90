! dgejsv: Jacobi SVD with optional preconditioning. Smoke test —
! compare singular values for joba='G' (general full pivoting),
! jobu='U', jobv='V', no extra preconditioning.
program test_dgejsv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgejsv
    use ref_quad_lapack, only: dgejsv
    implicit none

    integer, parameter :: ms(*) = [16, 24]
    integer, parameter :: ns(*) = [8, 16]
    integer :: i, m, n, info_r, info_g, lwork
    real(ep), allocatable :: A0(:,:), A_r(:,:), A_g(:,:), U_r(:,:), V_r(:,:), U_g(:,:), V_g(:,:)
    real(ep), allocatable :: sva_r(:), sva_g(:), work(:)
    integer, allocatable :: iwork(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgejsv', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 154001 + 47 * i)
        allocate(A_r(m,n), A_g(m,n), U_r(m,n), V_r(n,n), U_g(m,n), V_g(n,n))
        allocate(sva_r(n), sva_g(n))
        A_r = A0; A_g = A0
        lwork = max(7, 7*n + 4*n*n + 6 + 2*max(m,n))
        allocate(work(lwork), iwork(max(1, m+3*n)))
        call dgejsv('G', 'U', 'V', 'N', 'N', 'N', m, n, A_r, m, sva_r, &
                    U_r, m, V_r, n, work, lwork, iwork, info_r)
        deallocate(work, iwork)
        call target_dgejsv('G', 'U', 'V', 'N', 'N', 'N', m, n, A_g, m, sva_g, &
                           U_g, m, V_g, n, info_g)
        err = max_rel_err_vec(sva_g, sva_r)
        tol = 256.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_r, A_g, U_r, V_r, U_g, V_g, sva_r, sva_g)
    end do
    call report_finalize()
end program test_dgejsv
