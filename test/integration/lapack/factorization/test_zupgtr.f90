program test_zupgtr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hermitian_matrix_quad, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zupgtr
    use ref_quad_lapack, only: zhptrd, zupgtr
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, np, info
    complex(ep), allocatable :: A(:,:), AP(:), tau(:), work(:), Q_ref(:,:), Q_got(:,:)
    real(ep), allocatable :: D(:), E(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zupgtr', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_hermitian_matrix_quad(n, A, seed = 230071 + 47 * i)
        allocate(AP(np), D(n), E(n-1), tau(n-1), work(n-1))
        call pack_herm_packed_quad('U', n, A, AP)
        call zhptrd('U', n, AP, D, E, tau, info)
        allocate(Q_ref(n,n), Q_got(n,n))
        call zupgtr('U', n, AP, tau, Q_ref, n, work, info)
        call target_zupgtr('U', n, AP, tau, Q_got, n, info)
        err = max_rel_err_mat_z(Q_got, Q_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AP, D, E, tau, work, Q_ref, Q_got)
    end do
    call report_finalize()
end program test_zupgtr
