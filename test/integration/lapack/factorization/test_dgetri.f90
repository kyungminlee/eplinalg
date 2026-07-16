! dgetri: inverse from an LU factor. Factor first via dgetrf
! (same algorithm on both sides), then invert.
program test_dgetri
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgetrf, target_dgetri
    use ref_quad_lapack, only: dgetrf, dgetri
    implicit none

    integer, parameter :: ns(*) = [8, 32, 64]
    integer :: i, n, info, lwork
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgetri', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A0, seed = 23001 + 47 * i)
        ! Diagonal shift to ensure invertibility for random matrices.
        block
            integer :: j
            do j = 1, n
                A0(j, j) = A0(j, j) + real(n, ep)
            end do
        end block
        allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
        A_ref = A0; A_got = A0
        call dgetrf(n, n, A_ref, n, ipiv_ref, info)
        call target_dgetrf(n, n, A_got, n, ipiv_got, info)
        call dgetri(n, A_ref, n, ipiv_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgetri(n, A_ref, n, ipiv_ref, work, lwork, info)
        deallocate(work)
        call target_dgetri(n, A_got, n, ipiv_got, info)
        err = max_rel_err_mat(A_got, A_ref)
        tol = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, ipiv_ref, ipiv_got)
    end do
    call report_finalize()
end program test_dgetri
