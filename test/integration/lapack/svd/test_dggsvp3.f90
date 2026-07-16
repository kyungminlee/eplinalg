! dggsvp3: GSVD preprocessing (3-stage variant). Compare K, L counts +
! Frobenius norm of the leading (k+l) diagonal of the post-preprocessing
! A — invariant under sign-flips of the orthogonal factors.
program test_dggsvp3
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dggsvp3
    use ref_quad_lapack, only: dggsvp3, dlange
    implicit none

    integer, parameter :: ms(*) = [10, 14]
    integer, parameter :: ns(*) = [8,  10]
    integer, parameter :: ps(*) = [6,  8]
    integer :: i, j, m, n, p, k_r, l_r, k_g, l_g, info, lwork
    real(ep), allocatable :: A0(:,:), B0(:,:), Ar(:,:), Br(:,:), Ag(:,:), Bg(:,:)
    real(ep), allocatable :: U(:,:), V(:,:), Q(:,:), tau(:), work(:), rwork(:)
    real(ep) :: tola, tolb, normA, normB, anormR, anormG, err, tol, wopt(1)
    integer, allocatable :: iwork(:)
    character(len=64) :: label

    call report_init('dggsvp3', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); p = ps(i)
        call gen_matrix_quad(m, n, A0, seed = 26501 + 83 * i)
        call gen_matrix_quad(p, n, B0, seed = 26511 + 83 * i)
        allocate(rwork(2*max(m,p)), Ar(m, n), Br(p, n), Ag(m, n), Bg(p, n))
        Ar = A0; Br = B0; Ag = A0; Bg = B0
        normA = dlange('1', m, n, A0, m, rwork)
        normB = dlange('1', p, n, B0, p, rwork)
        tola = real(max(m, n), ep) * normA * epsilon(1.0_ep)
        tolb = real(max(p, n), ep) * normB * epsilon(1.0_ep)
        allocate(U(m, m), V(p, p), Q(n, n), tau(n), iwork(n))
        call dggsvp3('U', 'V', 'Q', m, p, n, Ar, m, Br, p, tola, tolb, &
                     k_r, l_r, U, m, V, p, Q, n, iwork, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dggsvp3('U', 'V', 'Q', m, p, n, Ar, m, Br, p, tola, tolb, &
                     k_r, l_r, U, m, V, p, Q, n, iwork, tau, work, lwork, info)
        deallocate(work)
        call target_dggsvp3('U', 'V', 'Q', m, p, n, Ag, m, Bg, p, tola, tolb, &
                            k_g, l_g, U, m, V, p, Q, n, iwork, tau, info)
        if (k_r /= k_g .or. l_r /= l_g) then
            err = huge(1.0_ep)
        else
            anormR = 0.0_ep; anormG = 0.0_ep
            do j = 1, k_r + l_r
                anormR = anormR + Ar(j, n - (k_r + l_r) + j)**2
                anormG = anormG + Ag(j, n - (k_g + l_g) + j)**2
            end do
            anormR = sqrt(anormR); anormG = sqrt(anormG)
            err = rel_err_scalar(anormG, anormR)
        end if
        tol = 16.0_ep * real(max(m, n, p), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',p=', p
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, Ar, Br, Ag, Bg, U, V, Q, tau, iwork, rwork)
    end do
    call report_finalize()
end program test_dggsvp3
