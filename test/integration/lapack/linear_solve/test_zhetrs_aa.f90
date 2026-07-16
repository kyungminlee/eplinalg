program test_zhetrs_aa
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hermitian_matrix_quad, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zhetrf_aa, target_zhetrs_aa
    use ref_quad_lapack, only: zhetrf_aa, zhetrs_aa
    implicit none
    integer, parameter :: ns(*) = [8, 32, 64]
    integer, parameter :: nrhs  = 2
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, info, ju, lwork, j
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), A_got(:,:), B_ref(:,:), B_got(:,:), work(:)
    integer,     allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('zhetrs_aa', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 51101 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + cmplx(real(n, ep), 0.0_ep, ep); end do
        call gen_matrix_complex(n, nrhs, B0, seed = 51111 + 47 * i)
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n), B_ref(n,nrhs), B_got(n,nrhs))
            allocate(ipiv_ref(n), ipiv_got(n))
            A_ref = A0; A_got = A0
            call zhetrf_aa(uplos(ju), n, A_ref, n, ipiv_ref, wopt, -1, info)
            lwork = max(1, int(real(wopt(1), ep)))
            allocate(work(lwork))
            call zhetrf_aa(uplos(ju), n, A_ref, n, ipiv_ref, work, lwork, info)
            deallocate(work)
            call target_zhetrf_aa(uplos(ju), n, A_got, n, ipiv_got, info)
            B_ref = B0; B_got = B0
            call zhetrs_aa(uplos(ju), n, nrhs, A_ref, n, ipiv_ref, B_ref, n, wopt, -1, info)
            lwork = max(1, int(real(wopt(1), ep)))
            allocate(work(lwork))
            call zhetrs_aa(uplos(ju), n, nrhs, A_ref, n, ipiv_ref, B_ref, n, work, lwork, info)
            deallocate(work)
            call target_zhetrs_aa(uplos(ju), n, nrhs, A_got, n, ipiv_got, B_got, n, info)
            err = max_rel_err_mat_z(B_got, B_ref)
            tol = 16.0_ep * real(n, ep)**3 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, B_ref, B_got, ipiv_ref, ipiv_got)
        end do
    end do
    call report_finalize()
end program test_zhetrs_aa
