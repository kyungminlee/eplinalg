program test_dgesv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgesv
    use ref_quad_lapack, only: dgesv
    implicit none

    integer, parameter :: ns(*)   = [8, 32, 96]
    integer, parameter :: nrhs    = 3
    integer :: i, n, info_ref, info_got
    real(ep), allocatable :: A0(:,:), B0(:,:)
    real(ep), allocatable :: A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgesv', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n,    A0, seed = 1001 + 31 * i)
        call gen_matrix_quad(n, nrhs, B0, seed = 1011 + 31 * i)

        allocate(A_ref(n, n), B_ref(n, nrhs), A_got(n, n), B_got(n, nrhs))
        allocate(ipiv_ref(n), ipiv_got(n))
        A_ref = A0;  B_ref = B0;  A_got = A0;  B_got = B0

        call dgesv(n, nrhs, A_ref, n, ipiv_ref, B_ref, n, info_ref)
        call target_dgesv(n, nrhs, A_got, n, ipiv_got, B_got, n, info_got)

        err = max_rel_err_mat(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs
        call report_case(trim(label), err, tol)

        deallocate(A_ref, B_ref, A_got, B_got, ipiv_ref, ipiv_got)
    end do
    call report_finalize()
end program test_dgesv
