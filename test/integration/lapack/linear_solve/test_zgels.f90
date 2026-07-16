program test_zgels
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgels
    use ref_quad_lapack, only: zgels
    implicit none

    integer, parameter :: ms(*) = [16, 24]
    integer, parameter :: ns(*) = [8, 16]
    integer, parameter :: nrhs  = 2
    integer :: i, m, n, ldb, info, lwork, krows
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: B_ref(:,:), B_got(:,:), work(:)
    complex(ep), allocatable :: Brblock(:,:), Bgblock(:,:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgels', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); ldb = max(m, n)
        call gen_matrix_complex(m, n, A0, seed = 31001 + 47 * i)
        call gen_matrix_complex(ldb, nrhs, B0, seed = 31011 + 47 * i)
        allocate(A_ref(m,n), A_got(m,n), B_ref(ldb,nrhs), B_got(ldb,nrhs))
        A_ref = A0; A_got = A0; B_ref = B0; B_got = B0
        call zgels('N', m, n, nrhs, A_ref, m, B_ref, ldb, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgels('N', m, n, nrhs, A_ref, m, B_ref, ldb, work, lwork, info)
        deallocate(work)
        call target_zgels('N', m, n, nrhs, A_got, m, B_got, ldb, info)
        krows = n
        allocate(Brblock(krows, nrhs), Bgblock(krows, nrhs))
        Brblock = B_ref(1:krows, :); Bgblock = B_got(1:krows, :)
        err = max_rel_err_mat_z(Bgblock, Brblock)
        tol = 16.0_ep * real(max(m,n), ep)**3 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, B_ref, B_got, Brblock, Bgblock)
    end do
    call report_finalize()
end program test_zgels
