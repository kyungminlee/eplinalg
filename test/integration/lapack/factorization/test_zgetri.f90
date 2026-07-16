program test_zgetri
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgetrf, target_zgetri
    use ref_quad_lapack, only: zgetrf, zgetri
    implicit none

    integer, parameter :: ns(*) = [8, 32, 64]
    integer :: i, n, info, lwork, j
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), work(:)
    integer,     allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgetri', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A0, seed = 42001 + 47 * i)
        do j = 1, n
            A0(j, j) = A0(j, j) + cmplx(real(n, ep), 0.0_ep, ep)
        end do
        allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
        A_ref = A0; A_got = A0
        call zgetrf(n, n, A_ref, n, ipiv_ref, info)
        call target_zgetrf(n, n, A_got, n, ipiv_got, info)
        call zgetri(n, A_ref, n, ipiv_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgetri(n, A_ref, n, ipiv_ref, work, lwork, info)
        deallocate(work)
        call target_zgetri(n, A_got, n, ipiv_got, info)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, ipiv_ref, ipiv_got)
    end do
    call report_finalize()
end program test_zgetri
