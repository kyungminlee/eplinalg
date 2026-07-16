program test_dptrfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_vector_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, &
                                target_dpttrf, target_dpttrs, target_dptrfs
    use ref_quad_lapack, only: dpttrf, dpttrs, dptrfs
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer, parameter :: nrhs = 2
    integer :: i, n, info, j
    real(ep), allocatable :: d0(:), e0(:), B0(:,:)
    real(ep), allocatable :: df_ref(:), ef_ref(:), df_got(:), ef_got(:)
    real(ep), allocatable :: X_ref(:,:), X_got(:,:), work(:), ferr(:), berr(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dptrfs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n,   d0, seed = 144001 + 47 * i)
        call gen_vector_quad(n-1, e0, seed = 144011 + 47 * i)
        call gen_matrix_quad(n, nrhs, B0, seed = 144021 + 47 * i)
        do j = 1, n; d0(j) = abs(d0(j)) + real(4, ep); end do
        allocate(df_ref(n), ef_ref(n-1), df_got(n), ef_got(n-1))
        allocate(X_ref(n, nrhs), X_got(n, nrhs), work(2*n), ferr(nrhs), berr(nrhs))
        df_ref = d0; ef_ref = e0; df_got = d0; ef_got = e0
        call dpttrf(n, df_ref, ef_ref, info)
        call target_dpttrf(n, df_got, ef_got, info)
        X_ref = B0; X_got = B0
        call dpttrs(n, nrhs, df_ref, ef_ref, X_ref, n, info)
        call target_dpttrs(n, nrhs, df_got, ef_got, X_got, n, info)
        call dptrfs(n, nrhs, d0, e0, df_ref, ef_ref, B0, n, X_ref, n, &
                    ferr, berr, work, info)
        call target_dptrfs(n, nrhs, d0, e0, df_got, ef_got, B0, n, X_got, n, &
                           ferr, berr, info)
        err = max_rel_err_mat(X_got, X_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(d0, e0, B0, df_ref, ef_ref, df_got, ef_got, X_ref, X_got, work, ferr, berr)
    end do
    call report_finalize()
end program test_dptrfs
