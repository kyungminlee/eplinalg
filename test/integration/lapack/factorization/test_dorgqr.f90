! Tests target_dorgqr by first factoring A via dgeqrf (same algorithm
! on both sides), then calling dorgqr to form Q explicitly.
program test_dorgqr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, &
                                target_dgeqrf, target_dorgqr
    use ref_quad_lapack, only: dgeqrf, dorgqr
    implicit none

    integer, parameter :: ms(*) = [16, 48, 96]
    integer, parameter :: ns(*) = [8,  24, 64]
    integer :: i, m, n, info, lwork
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: tau_ref(:), tau_got(:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dorgqr', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 8001 + 61 * i)

        allocate(A_ref(m, n), A_got(m, n), tau_ref(min(m, n)), tau_got(min(m, n)))
        A_ref = A0;  A_got = A0

        ! Factor
        call dgeqrf(m, n, A_ref, m, tau_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgeqrf(m, n, A_ref, m, tau_ref, work, lwork, info)
        deallocate(work)
        call target_dgeqrf(m, n, A_got, m, tau_got, info)

        ! Form Q (n columns, since n <= m)
        call dorgqr(m, n, n, A_ref, m, tau_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dorgqr(m, n, n, A_ref, m, tau_ref, work, lwork, info)
        deallocate(work)
        call target_dorgqr(m, n, n, A_got, m, tau_got, info)

        err = max_rel_err_mat(A_got, A_ref)
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)

        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_dorgqr
