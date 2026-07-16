! zggglm: complex generalized linear model minimization. Smoke test.
program test_zggglm
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec_z
    use test_data,       only: gen_matrix_complex, gen_vector_complex
    use target_lapack,   only: target_name, target_eps, target_zggglm
    use ref_quad_lapack, only: zggglm
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, m, p, info, lwork
    complex(ep), allocatable :: A(:,:), B(:,:), d(:)
    complex(ep), allocatable :: A_r(:,:), B_r(:,:), d_r(:), x_r(:), y_r(:), work(:)
    complex(ep), allocatable :: A_g(:,:), B_g(:,:), d_g(:), x_g(:), y_g(:)
    complex(ep) :: wopt(1)
    real(ep) :: err_x, err_y, err, tol
    character(len=48) :: label

    call report_init('zggglm', target_name)
    do i = 1, size(ns)
        n = ns(i); m = n / 2; p = n
        call gen_matrix_complex(n, m, A, seed = 142101 + 47 * i)
        call gen_matrix_complex(n, p, B, seed = 142111 + 47 * i)
        allocate(d(n)); call gen_vector_complex(n, d, seed = 142121 + 47 * i)
        allocate(A_r(n,m), B_r(n,p), d_r(n), x_r(m), y_r(p))
        allocate(A_g(n,m), B_g(n,p), d_g(n), x_g(m), y_g(p))
        A_r = A; B_r = B; d_r = d
        A_g = A; B_g = B; d_g = d
        call zggglm(n, m, p, A_r, n, B_r, n, d_r, x_r, y_r, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zggglm(n, m, p, A_r, n, B_r, n, d_r, x_r, y_r, work, lwork, info)
        deallocate(work)
        call target_zggglm(n, m, p, A_g, n, B_g, n, d_g, x_g, y_g, info)
        err_x = max_rel_err_vec_z(x_g, x_r)
        err_y = max_rel_err_vec_z(y_g, y_r)
        err = max(err_x, err_y)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, d, A_r, B_r, d_r, x_r, y_r, A_g, B_g, d_g, x_g, y_g)
    end do
    call report_finalize()
end program test_zggglm
