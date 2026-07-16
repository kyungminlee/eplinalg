! zunhr_col: Householder reconstruction (complex).
program test_zunhr_col
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, max_rel_err_vec_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zunhr_col
    use ref_quad_lapack, only: zlatsqr, zungtsqr, zunhr_col
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [4,  8]
    integer, parameter :: mbs(*) = [8, 16]
    integer, parameter :: nbs(*) = [2, 4]
    integer :: i, m, n, mb, nb, info, lwork, ntcols
    complex(ep), allocatable :: A0(:,:), Q(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: T_arr(:,:), T_ref(:,:), T_got(:,:)
    complex(ep), allocatable :: D_ref(:), D_got(:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zunhr_col', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); mb = mbs(i); nb = nbs(i)
        call gen_matrix_complex(m, n, A0, seed = 23451 + 97 * i)
        ntcols = n * max(1, (m - n + (mb - n) - 1) / max(1, mb - n))
        allocate(Q(m, n), T_arr(nb, ntcols))
        Q = A0
        call zlatsqr(m, n, mb, nb, Q, m, T_arr, nb, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zlatsqr(m, n, mb, nb, Q, m, T_arr, nb, work, lwork, info)
        deallocate(work)
        call zungtsqr(m, n, mb, nb, Q, m, T_arr, nb, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zungtsqr(m, n, mb, nb, Q, m, T_arr, nb, work, lwork, info)
        deallocate(work)
        allocate(A_ref(m, n), A_got(m, n), T_ref(nb, n), T_got(nb, n), D_ref(n), D_got(n))
        A_ref = Q; A_got = Q
        T_ref = (0.0_ep, 0.0_ep); T_got = (0.0_ep, 0.0_ep)
        call zunhr_col(m, n, nb, A_ref, m, T_ref, nb, D_ref, info)
        call target_zunhr_col(m, n, nb, A_got, m, T_got, nb, D_got, info)
        err = max(max_rel_err_mat_z(A_got, A_ref), max_rel_err_mat_z(T_got, T_ref), &
                  max_rel_err_vec_z(D_got, D_ref))
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',nb=', nb
        call report_case(trim(label), err, tol)
        deallocate(A0, Q, T_arr, A_ref, A_got, T_ref, T_got, D_ref, D_got)
    end do
    call report_finalize()
end program test_zunhr_col
