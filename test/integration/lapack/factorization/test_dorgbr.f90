! Square (m=n=k) case: VECT='Q' generates the M-by-N orthogonal Q
! from dgebrd's overwritten A and tauq.
program test_dorgbr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dorgbr
    use ref_quad_lapack, only: dgebrd, dorgbr
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, lwork
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:)
    real(ep), allocatable :: D(:), E(:), tauq(:), taup(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dorgbr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A0, seed = 240021 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), D(n), E(n-1), tauq(n), taup(n))
        A_ref = A0
        call dgebrd(n, n, A_ref, n, D, E, tauq, taup, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgebrd(n, n, A_ref, n, D, E, tauq, taup, work, lwork, info)
        deallocate(work)
        A_got = A_ref
        call dorgbr('Q', n, n, n, A_ref, n, tauq, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dorgbr('Q', n, n, n, A_ref, n, tauq, work, lwork, info)
        call target_dorgbr('Q', n, n, n, A_got, n, tauq, info)
        err = max_rel_err_mat(A_got, A_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, D, E, tauq, taup, work)
    end do
    call report_finalize()
end program test_dorgbr
