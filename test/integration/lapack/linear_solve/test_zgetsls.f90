program test_zgetsls
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgetsls
    use ref_quad_lapack, only: zgetsls
    implicit none

    integer, parameter :: ms(*)   = [16, 32, 48]
    integer, parameter :: ns(*)   = [12, 20, 28]
    integer, parameter :: nrhs    = 3
    integer :: i, m, n, info, lwork
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgetsls', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A0, seed = 390021 + 47 * i)
        call gen_matrix_complex(m, nrhs, B0, seed = 390031 + 47 * i)
        allocate(A_ref(m, n), B_ref(m, nrhs), A_got(m, n), B_got(m, nrhs))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        call zgetsls('N', m, n, nrhs, A_ref, m, B_ref, m, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgetsls('N', m, n, nrhs, A_ref, m, B_ref, m, work, lwork, info)
        call target_zgetsls('N', m, n, nrhs, A_got, m, B_got, m, info)
        err = max_rel_err_mat_z(B_got(1:n, 1:nrhs), B_ref(1:n, 1:nrhs))
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, A_got, B_got, work)
    end do
    call report_finalize()
end program test_zgetsls
