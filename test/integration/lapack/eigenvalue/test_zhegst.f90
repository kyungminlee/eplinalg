! zhegst: reduce complex Hermitian gen-eig to standard form.
program test_zhegst
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hpd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zhegst
    use ref_quad_lapack, only: zhegst, zpotrf
    implicit none

    integer, parameter :: ns(*) = [12, 24, 48]
    integer :: i, n, info
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_r(:,:), A_g(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zhegst', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hpd_matrix_quad(n, A0, seed = 21501 + 47 * i)
        call gen_hpd_matrix_quad(n, B0, seed = 21511 + 47 * i)
        call zpotrf('U', n, B0, n, info)
        allocate(A_r(n,n), A_g(n,n))
        A_r = A0; A_g = A0
        call zhegst(1, 'U', n, A_r, n, B0, n, info)
        call target_zhegst(1, 'U', n, A_g, n, B0, n, info)
        err = max_rel_err_mat_z(A_g, A_r)
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_r, A_g)
    end do
    call report_finalize()
end program test_zhegst
