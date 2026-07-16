! dgetsqrhrt: TSQR + Householder reconstruction.
program test_dgetsqrhrt
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgetsqrhrt
    use ref_quad_lapack, only: dgetsqrhrt
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [4,  8]
    integer :: i, m, n, mb1, nb1, nb2, info, lwork
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), T_ref(:,:), T_got(:,:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgetsqrhrt', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); mb1 = 2*n; nb1 = 2; nb2 = 2  ! MB1 > N is required
        call gen_matrix_quad(m, n, A0, seed = 23501 + 67 * i)
        allocate(A_ref(m, n), A_got(m, n), T_ref(nb2, n), T_got(nb2, n))
        A_ref = A0; A_got = A0
        T_ref = 0.0_ep; T_got = 0.0_ep
        call dgetsqrhrt(m, n, mb1, nb1, nb2, A_ref, m, T_ref, nb2, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgetsqrhrt(m, n, mb1, nb1, nb2, A_ref, m, T_ref, nb2, work, lwork, info)
        deallocate(work)
        call target_dgetsqrhrt(m, n, mb1, nb1, nb2, A_got, m, T_got, nb2, info)
        err = max(max_rel_err_mat(A_got, A_ref), max_rel_err_mat(T_got, T_ref))
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, T_ref, T_got)
    end do
    call report_finalize()
end program test_dgetsqrhrt
