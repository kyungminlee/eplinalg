program test_zsytf2_rk
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, max_rel_err_vec_z
    use test_data,       only: gen_complex_symmetric_quad
    use target_lapack,   only: target_name, target_eps, target_zsytf2_rk
    use ref_quad_lapack, only: zsytf2_rk
    implicit none

    integer, parameter :: ns(*) = [8, 16, 32]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, info, ju, j
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), e_ref(:), e_got(:)
    integer,     allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zsytf2_rk', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_complex_symmetric_quad(n, A0, seed = 59101 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + cmplx(real(n, ep), 0.0_ep, ep); end do
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n), e_ref(n), e_got(n), ipiv_ref(n), ipiv_got(n))
            A_ref = A0; A_got = A0
            call zsytf2_rk(uplos(ju), n, A_ref, n, e_ref, ipiv_ref, info)
            call target_zsytf2_rk(uplos(ju), n, A_got, n, e_got, ipiv_got, info)
            err = max(max_rel_err_mat_z(A_got, A_ref), max_rel_err_vec_z(e_got, e_ref))
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, e_ref, e_got, ipiv_ref, ipiv_got)
        end do
    end do
    call report_finalize()
end program test_zsytf2_rk
