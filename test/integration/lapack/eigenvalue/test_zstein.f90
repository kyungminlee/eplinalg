! zstein: tridiagonal eigenvectors via inverse iteration (complex Z).
program test_zstein
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use target_lapack,   only: target_name, target_eps, target_zstein
    use ref_quad_lapack, only: dstebz, zstein
    implicit none

    integer, parameter :: ns(*) = [10, 32]
    integer :: i, n, m, nsplit, info, j
    real(ep), allocatable :: D(:), E(:), W(:), work(:)
    complex(ep), allocatable :: Z_r(:,:), Z_g(:,:)
    integer, allocatable :: iblock(:), isplit(:), iwork(:), ifail_r(:), ifail_g(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zstein', target_name)
    do i = 1, size(ns)
        n = ns(i)
        allocate(D(n), E(max(1, n-1)), W(n), iblock(n), isplit(n), iwork(3*n))
        do j = 1, n; D(j) = real(j, ep); end do
        do j = 1, n-1; E(j) = 0.5_ep; end do
        allocate(work(4*n))
        call dstebz('A', 'B', n, 0.0_ep, 0.0_ep, 0, 0, 0.0_ep, D, E, m, nsplit, &
                    W, iblock, isplit, work, iwork, info)
        deallocate(work)
        allocate(Z_r(n, m), Z_g(n, m), ifail_r(m), ifail_g(m))
        allocate(work(5*n))
        call zstein(n, D, E, m, W, iblock, isplit, Z_r, n, work, iwork, ifail_r, info)
        deallocate(work)
        call target_zstein(n, D, E, m, W, iblock, isplit, Z_g, n, ifail_g, info)
        err = max_rel_err_mat_z(Z_g, Z_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(D, E, W, iblock, isplit, iwork, Z_r, Z_g, ifail_r, ifail_g)
    end do
    call report_finalize()
end program test_zstein
