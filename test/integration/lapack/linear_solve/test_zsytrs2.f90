! zsytrs2: blocked solve A*X=B given the zsytrf factorization (complex sym).
program test_zsytrs2
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_complex_symmetric_quad, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zsytrf, target_zsytrs2
    use ref_quad_lapack, only: zsytrf, zsytrs2
    implicit none

    integer, parameter :: ns(*)    = [8, 16, 32]
    integer, parameter :: nrhss(*) = [1, 3]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, ir, n, nrhs, info, ju, lwork, j
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: B_ref(:,:), B_got(:,:), work(:)
    integer,     allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('zsytrs2', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_complex_symmetric_quad(n, A0, seed = 26151 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + cmplx(real(n, ep), 0.0_ep, ep); end do
        do ir = 1, size(nrhss)
            nrhs = nrhss(ir)
            call gen_matrix_complex(n, nrhs, B0, seed = 26161 + 47 * i + ir)
            do ju = 1, size(uplos)
                allocate(A_ref(n,n), A_got(n,n), B_ref(n,nrhs), B_got(n,nrhs), &
                         ipiv_ref(n), ipiv_got(n))
                A_ref = A0; A_got = A0; B_ref = B0; B_got = B0
                call zsytrf(uplos(ju), n, A_ref, n, ipiv_ref, wopt, -1, info)
                lwork = max(1, int(real(wopt(1), ep)))
                allocate(work(lwork))
                call zsytrf(uplos(ju), n, A_ref, n, ipiv_ref, work, lwork, info)
                deallocate(work)
                call target_zsytrf(uplos(ju), n, A_got, n, ipiv_got, info)
                allocate(work(n))
                call zsytrs2(uplos(ju), n, nrhs, A_ref, n, ipiv_ref, B_ref, n, work, info)
                deallocate(work)
                call target_zsytrs2(uplos(ju), n, nrhs, A_got, n, ipiv_got, B_got, n, info)
                err = max_rel_err_mat_z(B_got, B_ref)
                tol = 16.0_ep * real(n, ep)**2 * target_eps
                write(label, '(a,a,a,i0,a,i0)') 'uplo=', uplos(ju), ',n=', n, ',nrhs=', nrhs
                call report_case(trim(label), err, tol)
                deallocate(A_ref, A_got, B_ref, B_got, ipiv_ref, ipiv_got)
            end do
            deallocate(B0)
        end do
        deallocate(A0)
    end do
    call report_finalize()
end program test_zsytrs2
