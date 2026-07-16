program test_zupmtr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hermitian_matrix_quad, gen_matrix_complex, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zupmtr
    use ref_quad_lapack, only: zhptrd, zupmtr
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: ncols = 5
    integer :: i, n, np, info
    complex(ep), allocatable :: A(:,:), AP(:), C0(:,:), C_ref(:,:), C_got(:,:), tau(:), work(:)
    real(ep), allocatable :: D(:), E(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zupmtr', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_hermitian_matrix_quad(n, A, seed = 230101 + 47 * i)
        call gen_matrix_complex(n, ncols, C0, seed = 230111 + 47 * i)
        allocate(AP(np), D(n), E(n-1), tau(n-1), work(ncols))
        call pack_herm_packed_quad('U', n, A, AP)
        call zhptrd('U', n, AP, D, E, tau, info)
        allocate(C_ref(n, ncols), C_got(n, ncols))
        C_ref = C0; C_got = C0
        call zupmtr('L', 'U', 'N', n, ncols, AP, tau, C_ref, n, work, info)
        call target_zupmtr('L', 'U', 'N', n, ncols, AP, tau, C_got, n, info)
        err = max_rel_err_mat_z(C_got, C_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AP, C0, D, E, tau, work, C_ref, C_got)
    end do
    call report_finalize()
end program test_zupmtr
