! zungtsqr_row: form Q from a TSQR factorization (row variant, complex).
program test_zungtsqr_row
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zungtsqr_row
    use ref_quad_lapack, only: zlatsqr, zungtsqr_row
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [4,  8]
    integer, parameter :: mbs(*) = [8, 16]
    integer, parameter :: nbs(*) = [2, 4]
    integer :: i, m, n, mb, nb, info, lwork, ntcols
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), T_arr(:,:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zungtsqr_row', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); mb = mbs(i); nb = nbs(i)
        call gen_matrix_complex(m, n, A0, seed = 23351 + 83 * i)
        ntcols = n * max(1, (m - n + (mb - n) - 1) / max(1, mb - n))
        allocate(A_ref(m, n), A_got(m, n), T_arr(nb, ntcols))
        A_ref = A0; A_got = A0
        call zlatsqr(m, n, mb, nb, A_ref, m, T_arr, nb, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zlatsqr(m, n, mb, nb, A_ref, m, T_arr, nb, work, lwork, info)
        deallocate(work)
        A_got = A_ref
        call zungtsqr_row(m, n, mb, nb, A_ref, m, T_arr, nb, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zungtsqr_row(m, n, mb, nb, A_ref, m, T_arr, nb, work, lwork, info)
        deallocate(work)
        call target_zungtsqr_row(m, n, mb, nb, A_got, m, T_arr, nb, info)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',mb=', mb, ',nb=', nb
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, T_arr, A0)
    end do
    call report_finalize()
end program test_zungtsqr_row
