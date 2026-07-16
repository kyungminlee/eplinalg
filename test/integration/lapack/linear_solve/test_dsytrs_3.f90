! dsytrs_3: solve A*X=B given the bounded BK factorization from sytrf_rk.
program test_dsytrs_3
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_symmetric_matrix_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsytrf_rk, target_dsytrs_3
    use ref_quad_lapack, only: dsytrf_rk, dsytrs_3
    implicit none

    integer, parameter :: ns(*)    = [8, 16, 32]
    integer, parameter :: nrhss(*) = [1, 3]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, ir, n, nrhs, info, ju, lwork, j
    real(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: B_ref(:,:), B_got(:,:), e_ref(:), e_got(:), work(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: wopt(1), err, tol
    character(len=64) :: label

    call report_init('dsytrs_3', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 26201 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + real(n, ep); end do
        do ir = 1, size(nrhss)
            nrhs = nrhss(ir)
            call gen_matrix_quad(n, nrhs, B0, seed = 26211 + 47 * i + ir)
            do ju = 1, size(uplos)
                allocate(A_ref(n,n), A_got(n,n), B_ref(n,nrhs), B_got(n,nrhs), &
                         e_ref(n), e_got(n), ipiv_ref(n), ipiv_got(n))
                A_ref = A0; A_got = A0; B_ref = B0; B_got = B0
                call dsytrf_rk(uplos(ju), n, A_ref, n, e_ref, ipiv_ref, wopt, -1, info)
                lwork = max(1, int(wopt(1)))
                allocate(work(lwork))
                call dsytrf_rk(uplos(ju), n, A_ref, n, e_ref, ipiv_ref, work, lwork, info)
                deallocate(work)
                call target_dsytrf_rk(uplos(ju), n, A_got, n, e_got, ipiv_got, info)
                call dsytrs_3(uplos(ju), n, nrhs, A_ref, n, e_ref, ipiv_ref, B_ref, n, info)
                call target_dsytrs_3(uplos(ju), n, nrhs, A_got, n, e_got, ipiv_got, B_got, n, info)
                err = max_rel_err_mat(B_got, B_ref)
                tol = 16.0_ep * real(n, ep)**2 * target_eps
                write(label, '(a,a,a,i0,a,i0)') 'uplo=', uplos(ju), ',n=', n, ',nrhs=', nrhs
                call report_case(trim(label), err, tol)
                deallocate(A_ref, A_got, B_ref, B_got, e_ref, e_got, ipiv_ref, ipiv_got)
            end do
            deallocate(B0)
        end do
        deallocate(A0)
    end do
    call report_finalize()
end program test_dsytrs_3
