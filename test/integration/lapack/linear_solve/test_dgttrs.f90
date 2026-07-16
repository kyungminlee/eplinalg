program test_dgttrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_vector_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgttrf, target_dgttrs
    use ref_quad_lapack, only: dgttrf, dgttrs
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer, parameter :: nrhs  = 2
    integer :: i, n, info, j
    real(ep), allocatable :: dl0(:), d0(:), du0(:), B0(:,:)
    real(ep), allocatable :: dl_ref(:), d_ref(:), du_ref(:), du2_ref(:)
    real(ep), allocatable :: dl_got(:), d_got(:), du_got(:), du2_got(:)
    real(ep), allocatable :: B_ref(:,:), B_got(:,:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgttrs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n-1, dl0, seed = 67001 + 47 * i)
        call gen_vector_quad(n,   d0,  seed = 67011 + 47 * i)
        call gen_vector_quad(n-1, du0, seed = 67021 + 47 * i)
        do j = 1, n; d0(j) = d0(j) + real(4, ep); end do
        call gen_matrix_quad(n, nrhs, B0, seed = 67031 + 47 * i)
        allocate(dl_ref(n-1), d_ref(n), du_ref(n-1), du2_ref(n-2))
        allocate(dl_got(n-1), d_got(n), du_got(n-1), du2_got(n-2))
        allocate(B_ref(n,nrhs), B_got(n,nrhs), ipiv_ref(n), ipiv_got(n))
        dl_ref = dl0; d_ref = d0; du_ref = du0
        dl_got = dl0; d_got = d0; du_got = du0
        call dgttrf(n, dl_ref, d_ref, du_ref, du2_ref, ipiv_ref, info)
        call target_dgttrf(n, dl_got, d_got, du_got, du2_got, ipiv_got, info)
        B_ref = B0; B_got = B0
        call dgttrs('N', n, nrhs, dl_ref, d_ref, du_ref, du2_ref, ipiv_ref, B_ref, n, info)
        call target_dgttrs('N', n, nrhs, dl_got, d_got, du_got, du2_got, ipiv_got, B_got, n, info)
        err = max_rel_err_mat(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(dl0, d0, du0, B0, dl_ref, d_ref, du_ref, du2_ref)
        deallocate(dl_got, d_got, du_got, du2_got, B_ref, B_got, ipiv_ref, ipiv_got)
    end do
    call report_finalize()
end program test_dgttrs
