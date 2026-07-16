! dsygst: reduce real symmetric gen-eig to standard form.
program test_dsygst
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_spd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsygst
    use ref_quad_lapack, only: dsygst, dpotrf
    implicit none

    integer, parameter :: ns(*) = [12, 24, 48]
    integer :: i, n, info
    real(ep), allocatable :: A0(:,:), B0(:,:), A_r(:,:), A_g(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dsygst', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_spd_matrix_quad(n, A0, seed = 21401 + 47 * i)
        call gen_spd_matrix_quad(n, B0, seed = 21411 + 47 * i)
        call dpotrf('U', n, B0, n, info)
        allocate(A_r(n,n), A_g(n,n))
        A_r = A0; A_g = A0
        call dsygst(1, 'U', n, A_r, n, B0, n, info)
        call target_dsygst(1, 'U', n, A_g, n, B0, n, info)
        err = max_rel_err_mat(A_g, A_r)
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_r, A_g)
    end do
    call report_finalize()
end program test_dsygst
