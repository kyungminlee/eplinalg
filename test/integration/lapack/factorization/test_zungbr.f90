program test_zungbr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zungbr
    use ref_quad_lapack, only: zgebrd, zungbr
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, lwork
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:), tauq(:), taup(:)
    real(ep), allocatable :: D(:), E(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zungbr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A0, seed = 240031 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), D(n), E(n-1), tauq(n), taup(n))
        A_ref = A0
        call zgebrd(n, n, A_ref, n, D, E, tauq, taup, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgebrd(n, n, A_ref, n, D, E, tauq, taup, work, lwork, info)
        deallocate(work)
        A_got = A_ref
        call zungbr('Q', n, n, n, A_ref, n, tauq, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zungbr('Q', n, n, n, A_ref, n, tauq, work, lwork, info)
        call target_zungbr('Q', n, n, n, A_got, n, tauq, info)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, D, E, tauq, taup, work)
    end do
    call report_finalize()
end program test_zungbr
