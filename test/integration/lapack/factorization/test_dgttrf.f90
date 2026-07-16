! Tridiagonal LU: 3 vectors (subdiag dl, diag d, superdiag du) plus
! a fill-in vector du2 of length n-2 written by the factorization.
program test_dgttrf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_dgttrf
    use ref_quad_lapack, only: dgttrf
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer :: i, n, info, j
    real(ep), allocatable :: dl0(:), d0(:), du0(:)
    real(ep), allocatable :: dl_ref(:), d_ref(:), du_ref(:), du2_ref(:)
    real(ep), allocatable :: dl_got(:), d_got(:), du_got(:), du2_got(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgttrf', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n-1, dl0, seed = 66001 + 47 * i)
        call gen_vector_quad(n,   d0,  seed = 66011 + 47 * i)
        call gen_vector_quad(n-1, du0, seed = 66021 + 47 * i)
        do j = 1, n
            d0(j) = d0(j) + real(4, ep)
        end do
        allocate(dl_ref(n-1), d_ref(n), du_ref(n-1), du2_ref(n-2))
        allocate(dl_got(n-1), d_got(n), du_got(n-1), du2_got(n-2))
        allocate(ipiv_ref(n), ipiv_got(n))
        dl_ref = dl0; d_ref = d0; du_ref = du0
        dl_got = dl0; d_got = d0; du_got = du0
        call dgttrf(n, dl_ref, d_ref, du_ref, du2_ref, ipiv_ref, info)
        call target_dgttrf(n, dl_got, d_got, du_got, du2_got, ipiv_got, info)
        err = max(max_rel_err_vec(d_got, d_ref), &
                  max_rel_err_vec(du_got, du_ref))
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(dl0, d0, du0, dl_ref, d_ref, du_ref, du2_ref)
        deallocate(dl_got, d_got, du_got, du2_got, ipiv_ref, ipiv_got)
    end do
    call report_finalize()
end program test_dgttrf
