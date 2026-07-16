program test_zungtr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zungtr
    use ref_quad_lapack, only: zhetrd, zungtr
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info, lwork
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), tau(:), work(:)
    real(ep), allocatable :: D(:), E(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zungtr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 230011 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), D(n), E(n-1), tau(n-1))
        A_ref = A0
        call zhetrd('U', n, A_ref, n, D, E, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zhetrd('U', n, A_ref, n, D, E, tau, work, lwork, info)
        deallocate(work)
        A_got = A_ref
        call zungtr('U', n, A_ref, n, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zungtr('U', n, A_ref, n, tau, work, lwork, info)
        call target_zungtr('U', n, A_got, n, tau, info)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, D, E, tau, work)
    end do
    call report_finalize()
end program test_zungtr
