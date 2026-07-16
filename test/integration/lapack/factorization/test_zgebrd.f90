program test_zgebrd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgebrd
    use ref_quad_lapack, only: zgebrd
    implicit none

    integer, parameter :: ms(*) = [16, 32, 48]
    integer, parameter :: ns(*) = [16, 24, 64]
    integer :: i, m, n, mn, info, lwork
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:)
    complex(ep), allocatable :: tauq_ref(:), taup_ref(:), tauq_got(:), taup_got(:)
    real(ep), allocatable :: D_ref(:), E_ref(:), D_got(:), E_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgebrd', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); mn = min(m, n)
        call gen_matrix_complex(m, n, A0, seed = 240011 + 47 * i)
        allocate(A_ref(m,n), A_got(m,n))
        allocate(D_ref(mn), E_ref(mn-1), tauq_ref(mn), taup_ref(mn))
        allocate(D_got(mn), E_got(mn-1), tauq_got(mn), taup_got(mn))
        A_ref = A0; A_got = A0
        call zgebrd(m, n, A_ref, m, D_ref, E_ref, tauq_ref, taup_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgebrd(m, n, A_ref, m, D_ref, E_ref, tauq_ref, taup_ref, work, lwork, info)
        call target_zgebrd(m, n, A_got, m, D_got, E_got, tauq_got, taup_got, info)
        err = max(max_rel_err_vec(D_got, D_ref), max_rel_err_vec(E_got, E_ref))
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, D_ref, E_ref, tauq_ref, taup_ref, &
                   D_got, E_got, tauq_got, taup_got, work)
    end do
    call report_finalize()
end program test_zgebrd
