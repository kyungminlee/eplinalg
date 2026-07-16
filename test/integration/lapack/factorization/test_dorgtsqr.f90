! dorgtsqr: form Q from a TSQR factorization (output of dlatsqr).
program test_dorgtsqr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dorgtsqr
    use ref_quad_lapack, only: dlatsqr, dorgtsqr
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [4,  8]
    integer, parameter :: mbs(*) = [8, 16]
    integer, parameter :: nbs(*) = [2, 4]
    integer :: i, m, n, mb, nb, info, lwork, ntcols
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), T_arr(:,:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dorgtsqr', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); mb = mbs(i); nb = nbs(i)
        call gen_matrix_quad(m, n, A0, seed = 23201 + 67 * i)
        ntcols = n * max(1, (m - n + (mb - n) - 1) / max(1, mb - n))
        allocate(A_ref(m, n), A_got(m, n), T_arr(nb, ntcols))
        A_ref = A0; A_got = A0
        ! TSQR factorization (only A_ref needed; same A0 is the input)
        call dlatsqr(m, n, mb, nb, A_ref, m, T_arr, nb, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dlatsqr(m, n, mb, nb, A_ref, m, T_arr, nb, work, lwork, info)
        deallocate(work)
        A_got = A_ref       ! same TSQR factor for both
        ! Now form Q
        call dorgtsqr(m, n, mb, nb, A_ref, m, T_arr, nb, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dorgtsqr(m, n, mb, nb, A_ref, m, T_arr, nb, work, lwork, info)
        deallocate(work)
        call target_dorgtsqr(m, n, mb, nb, A_got, m, T_arr, nb, info)
        err = max_rel_err_mat(A_got, A_ref)
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',mb=', mb, ',nb=', nb
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, T_arr, A0)
    end do
    call report_finalize()
end program test_dorgtsqr
