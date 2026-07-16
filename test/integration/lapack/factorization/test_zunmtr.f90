program test_zunmtr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hermitian_matrix_quad, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zunmtr
    use ref_quad_lapack, only: zhetrd, zunmtr
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: ncols = 5
    integer :: i, n, info, lwork
    complex(ep), allocatable :: A(:,:), C0(:,:), C_ref(:,:), C_got(:,:), tau(:), work(:)
    real(ep), allocatable :: D(:), E(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zunmtr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A, seed = 230041 + 47 * i)
        call gen_matrix_complex(n, ncols, C0, seed = 230051 + 47 * i)
        allocate(D(n), E(n-1), tau(n-1))
        call zhetrd('U', n, A, n, D, E, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zhetrd('U', n, A, n, D, E, tau, work, lwork, info)
        deallocate(work)
        allocate(C_ref(n, ncols), C_got(n, ncols))
        C_ref = C0; C_got = C0
        call zunmtr('L', 'U', 'N', n, ncols, A, n, tau, C_ref, n, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zunmtr('L', 'U', 'N', n, ncols, A, n, tau, C_ref, n, work, lwork, info)
        call target_zunmtr('L', 'U', 'N', n, ncols, A, n, tau, C_got, n, info)
        err = max_rel_err_mat_z(C_got, C_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, C0, D, E, tau, work, C_ref, C_got)
    end do
    call report_finalize()
end program test_zunmtr
