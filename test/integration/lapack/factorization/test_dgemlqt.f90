program test_dgemlqt
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgemlqt
    use ref_quad_lapack, only: dgelqt, dgemlqt
    implicit none

    integer, parameter :: ms(*) = [12, 20, 28]
    integer, parameter :: ns(*) = [16, 32, 48]
    integer, parameter :: ncols = 5
    integer, parameter :: mb    = 4
    integer :: i, m, n, info
    real(ep), allocatable :: V(:,:), T(:,:), C0(:,:), C_ref(:,:), C_got(:,:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgemlqt', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, V, seed = 380081 + 47 * i)
        allocate(T(mb, min(m, n)), work(mb*max(m, n)))
        call dgelqt(m, n, mb, V, m, T, mb, work, info)
        deallocate(work)
        call gen_matrix_quad(n, ncols, C0, seed = 380091 + 47 * i)
        allocate(C_ref(n, ncols), C_got(n, ncols), work(ncols*mb))
        C_ref = C0; C_got = C0
        call dgemlqt('L', 'N', n, ncols, m, mb, V, m, T, mb, C_ref, n, work, info)
        call target_dgemlqt('L', 'N', n, ncols, m, mb, V, m, T, mb, C_got, n, info)
        err = max_rel_err_mat(C_got, C_ref)
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(V, T, C0, C_ref, C_got, work)
    end do
    call report_finalize()
end program test_dgemlqt
