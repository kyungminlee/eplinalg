program test_dptsv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_vector_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dptsv
    use ref_quad_lapack, only: dptsv
    implicit none
    integer, parameter :: ns(*) = [10, 32, 64]
    integer, parameter :: nrhs  = 2
    integer :: i, n, info, j
    real(ep), allocatable :: d0(:), e0(:), B0(:,:)
    real(ep), allocatable :: d_r(:), e_r(:), B_r(:,:)
    real(ep), allocatable :: d_g(:), e_g(:), B_g(:,:)
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('dptsv', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n,   d0, seed = 81001 + 19 * i)
        call gen_vector_quad(n-1, e0, seed = 81011 + 19 * i)
        do j = 1, n; d0(j) = abs(d0(j)) + 4.0_ep; end do
        call gen_matrix_quad(n, nrhs, B0, seed = 81021 + 19 * i)
        allocate(d_r(n), e_r(max(1,n-1)), B_r(n, nrhs))
        allocate(d_g(n), e_g(max(1,n-1)), B_g(n, nrhs))
        d_r = d0; e_r = e0; B_r = B0
        d_g = d0; e_g = e0; B_g = B0
        call dptsv(n, nrhs, d_r, e_r, B_r, n, info)
        call target_dptsv(n, nrhs, d_g, e_g, B_g, n, info)
        err = max_rel_err_mat(B_g, B_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(d_r, e_r, B_r, d_g, e_g, B_g)
    end do
    call report_finalize()
end program test_dptsv
