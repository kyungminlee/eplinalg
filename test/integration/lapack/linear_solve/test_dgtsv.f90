program test_dgtsv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_vector_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgtsv
    use ref_quad_lapack, only: dgtsv
    implicit none
    integer, parameter :: ns(*) = [10, 32, 64]
    integer, parameter :: nrhs  = 2
    integer :: i, n, info, j
    real(ep), allocatable :: dl0(:), d0(:), du0(:), B0(:,:)
    real(ep), allocatable :: dl_r(:), d_r(:), du_r(:), B_r(:,:)
    real(ep), allocatable :: dl_g(:), d_g(:), du_g(:), B_g(:,:)
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('dgtsv', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n,   d0, seed = 80001 + 19 * i)
        call gen_vector_quad(n-1, dl0, seed = 80011 + 19 * i)
        call gen_vector_quad(n-1, du0, seed = 80021 + 19 * i)
        do j = 1, n; d0(j) = d0(j) + 4.0_ep; end do
        call gen_matrix_quad(n, nrhs, B0, seed = 80031 + 19 * i)
        allocate(dl_r(max(1,n-1)), d_r(n), du_r(max(1,n-1)), B_r(n, nrhs))
        allocate(dl_g(max(1,n-1)), d_g(n), du_g(max(1,n-1)), B_g(n, nrhs))
        dl_r = dl0; d_r = d0; du_r = du0; B_r = B0
        dl_g = dl0; d_g = d0; du_g = du0; B_g = B0
        call dgtsv(n, nrhs, dl_r, d_r, du_r, B_r, n, info)
        call target_dgtsv(n, nrhs, dl_g, d_g, du_g, B_g, n, info)
        err = max_rel_err_mat(B_g, B_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(dl_r, d_r, du_r, B_r, dl_g, d_g, du_g, B_g)
    end do
    call report_finalize()
end program test_dgtsv
