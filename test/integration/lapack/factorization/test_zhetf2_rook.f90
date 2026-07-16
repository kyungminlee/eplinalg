program test_zhetf2_rook
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zhetf2_rook
    use ref_quad_lapack, only: zhetf2_rook
    implicit none

    integer, parameter :: ns(*) = [8, 16, 32]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, info, ju, j
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    integer,     allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zhetf2_rook', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 44201 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + cmplx(real(n, ep), 0.0_ep, ep); end do
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
            A_ref = A0; A_got = A0
            call zhetf2_rook(uplos(ju), n, A_ref, n, ipiv_ref, info)
            call target_zhetf2_rook(uplos(ju), n, A_got, n, ipiv_got, info)
            err = max_rel_err_mat_z(A_got, A_ref)
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, ipiv_ref, ipiv_got)
        end do
    end do
    call report_finalize()
end program test_zhetf2_rook
