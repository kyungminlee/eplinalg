! dgels: minimum-norm least-squares solution. The solution lives in
! the leading n-by-nrhs (overdetermined) or m-by-nrhs (underdetermined)
! block of B. Compare only that valid block.
program test_dgels
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgels
    use ref_quad_lapack, only: dgels
    implicit none

    integer, parameter :: ms(*) = [16, 24]   ! overdetermined
    integer, parameter :: ns(*) = [8, 16]
    integer, parameter :: nrhs  = 2
    integer :: i, m, n, ldb, info, lwork, krows
    real(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), A_got(:,:), B_ref(:,:), B_got(:,:), work(:)
    real(ep), allocatable :: Brblock(:,:), Bgblock(:,:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgels', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); ldb = max(m, n)
        call gen_matrix_quad(m, n,    A0, seed = 18001 + 47 * i)
        call gen_matrix_quad(ldb, nrhs, B0, seed = 18011 + 47 * i)
        allocate(A_ref(m,n), A_got(m,n), B_ref(ldb,nrhs), B_got(ldb,nrhs))
        A_ref = A0; A_got = A0; B_ref = B0; B_got = B0
        call dgels('N', m, n, nrhs, A_ref, m, B_ref, ldb, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgels('N', m, n, nrhs, A_ref, m, B_ref, ldb, work, lwork, info)
        deallocate(work)
        call target_dgels('N', m, n, nrhs, A_got, m, B_got, ldb, info)
        krows = n   ! overdetermined: solution in first n rows
        allocate(Brblock(krows, nrhs), Bgblock(krows, nrhs))
        Brblock = B_ref(1:krows, :); Bgblock = B_got(1:krows, :)
        err = max_rel_err_mat(Bgblock, Brblock)
        tol = 16.0_ep * real(max(m,n), ep)**3 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, B_ref, B_got, Brblock, Bgblock)
    end do
    call report_finalize()
end program test_dgels
