! zgelst: complex T-factor variant of zgels (LAPACK 3.10+).
program test_zgelst
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgelst
    use ref_quad_lapack, only: zgelst
    implicit none

    integer, parameter :: ms(*) = [16, 24]
    integer, parameter :: ns(*) = [8, 16]
    integer, parameter :: nrhs  = 2
    integer :: i, m, n, ldb, info, lwork
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_r(:,:), B_r(:,:), A_g(:,:), B_g(:,:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgelst', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); ldb = max(m, n)
        call gen_matrix_complex(m, n,    A0, seed = 153101 + 47 * i)
        call gen_matrix_complex(ldb, nrhs, B0, seed = 153111 + 47 * i)
        allocate(A_r(m,n), B_r(ldb,nrhs), A_g(m,n), B_g(ldb,nrhs))
        A_r = A0; B_r = B0; A_g = A0; B_g = B0
        call zgelst('N', m, n, nrhs, A_r, m, B_r, ldb, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgelst('N', m, n, nrhs, A_r, m, B_r, ldb, work, lwork, info)
        deallocate(work)
        call target_zgelst('N', m, n, nrhs, A_g, m, B_g, ldb, info)
        err = max_rel_err_mat_z(B_g(1:n,:), B_r(1:n,:))
        tol = 16.0_ep * real(max(m,n), ep)**3 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_r, B_r, A_g, B_g)
    end do
    call report_finalize()
end program test_zgelst
