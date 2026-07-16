! zggsvd3: generalized SVD (complex, 3-stage). Compare alpha/beta.
program test_zggsvd3
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zggsvd3
    use ref_quad_lapack, only: zggsvd3
    implicit none

    integer, parameter :: ms(*) = [10, 14]
    integer, parameter :: ns(*) = [8,  10]
    integer, parameter :: ps(*) = [6,  8]
    integer :: i, m, n, p, k_r, l_r, k_g, l_g, info
    complex(ep), allocatable :: A0(:,:), B0(:,:), Ar(:,:), Br(:,:), Ag(:,:), Bg(:,:)
    real(ep),    allocatable :: alpha_r(:), beta_r(:), alpha_g(:), beta_g(:)
    complex(ep), allocatable :: U(:,:), V(:,:), Q(:,:), work(:)
    real(ep),    allocatable :: rwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    integer, allocatable :: iwork(:)
    integer :: lwork
    character(len=64) :: label

    call report_init('zggsvd3', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); p = ps(i)
        call gen_matrix_complex(m, n, A0, seed = 26451 + 83 * i)
        call gen_matrix_complex(p, n, B0, seed = 26461 + 83 * i)
        allocate(Ar(m, n), Br(p, n), Ag(m, n), Bg(p, n))
        Ar = A0; Br = B0; Ag = A0; Bg = B0
        allocate(alpha_r(n), beta_r(n), alpha_g(n), beta_g(n), &
                 U(m, m), V(p, p), Q(n, n), iwork(n), rwork(2*n))
        call zggsvd3('U', 'V', 'Q', m, n, p, k_r, l_r, Ar, m, Br, p, &
                     alpha_r, beta_r, U, m, V, p, Q, n, wopt, -1, rwork, iwork, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zggsvd3('U', 'V', 'Q', m, n, p, k_r, l_r, Ar, m, Br, p, &
                     alpha_r, beta_r, U, m, V, p, Q, n, work, lwork, rwork, iwork, info)
        deallocate(work)
        call target_zggsvd3('U', 'V', 'Q', m, n, p, k_g, l_g, Ag, m, Bg, p, &
                            alpha_g, beta_g, U, m, V, p, Q, n, iwork, info)
        if (k_r /= k_g .or. l_r /= l_g) then
            err = huge(1.0_ep)
        else
            err = max(max_rel_err_vec(alpha_g(1:k_r+l_r), alpha_r(1:k_r+l_r)), &
                      max_rel_err_vec(beta_g(1:k_r+l_r), beta_r(1:k_r+l_r)))
        end if
        tol = 16.0_ep * real(max(m, n, p), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',p=', p
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, Ar, Br, Ag, Bg, alpha_r, beta_r, alpha_g, beta_g, &
                   U, V, Q, iwork, rwork)
    end do
    call report_finalize()
end program test_zggsvd3
