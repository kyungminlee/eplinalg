program test_dopmtr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_symmetric_matrix_quad, gen_matrix_quad, pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, target_dopmtr
    use ref_quad_lapack, only: dsptrd, dopmtr
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: ncols = 5
    integer :: i, n, np, info
    real(ep), allocatable :: A(:,:), AP(:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep), allocatable :: D(:), E(:), tau(:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dopmtr', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_symmetric_matrix_quad(n, A, seed = 230081 + 47 * i)
        call gen_matrix_quad(n, ncols, C0, seed = 230091 + 47 * i)
        allocate(AP(np), D(n), E(n-1), tau(n-1), work(ncols))
        call pack_sym_packed_quad('U', n, A, AP)
        call dsptrd('U', n, AP, D, E, tau, info)
        allocate(C_ref(n, ncols), C_got(n, ncols))
        C_ref = C0; C_got = C0
        call dopmtr('L', 'U', 'N', n, ncols, AP, tau, C_ref, n, work, info)
        call target_dopmtr('L', 'U', 'N', n, ncols, AP, tau, C_got, n, info)
        err = max_rel_err_mat(C_got, C_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AP, C0, D, E, tau, work, C_ref, C_got)
    end do
    call report_finalize()
end program test_dopmtr
