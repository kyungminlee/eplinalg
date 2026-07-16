program test_zgbmv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_matrix_complex, gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zgbmv
    use ref_quad_blas, only: zgbmv
    implicit none

    integer, parameter :: cases(*)             = [20, 100, 64]
    character(len=1), parameter :: transes(*) = ['N', 'C', 'T']
    integer,          parameter :: kls(*)     = [2, 3, 1]
    integer,          parameter :: kus(*)     = [3, 2, 4]
    integer :: i, n, kl, ku, lda
    complex(ep), allocatable :: A(:,:), x(:), y0(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('zgbmv', target_name)
    do i = 1, size(cases)
        n   = cases(i)
        kl  = kls(i); ku = kus(i)
        lda = kl + ku + 1
        call gen_matrix_complex(lda, n, A,  seed = 1301 + 17 * i)
        call gen_vector_complex(n,      x,  seed = 1311 + 17 * i)
        call gen_vector_complex(n,      y0, seed = 1321 + 17 * i)
        alpha = cmplx(0.6_ep, 0.2_ep, ep)
        beta  = cmplx(0.4_ep, -0.1_ep, ep)
        allocate(y_ref(n), y_got(n))
        y_ref = y0; y_got = y0
        call zgbmv(transes(i), n, n, kl, ku, alpha, A, lda, x, 1, beta, y_ref, 1)
        call target_zgbmv(transes(i), n, n, kl, ku, alpha, A, lda, x, 1, beta, y_got, 1)
        err = max_rel_err_vec_z(y_got, y_ref)
        tol = 32.0_ep * real(kl + ku + 1, ep) * target_eps
        write(label, '(a,a,a,i0,a,i0,a,i0)') 'trans=', transes(i), ',n=', n, &
            ',kl=', kl, ',ku=', ku
        call report_case(trim(label), err, tol)
        deallocate(y_ref, y_got)
    end do
    call report_finalize()
end program test_zgbmv
