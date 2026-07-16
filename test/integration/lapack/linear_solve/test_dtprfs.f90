program test_dtprfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad, pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, target_dtptrs, target_dtprfs
    use ref_quad_lapack, only: dtptrs, dtprfs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs = 3
    integer :: i, n, info, np, j
    real(ep), allocatable :: A(:,:), AP(:), B0(:,:), X_ref(:,:), X_got(:,:), work(:)
    real(ep), allocatable :: ferr_ref(:), berr_ref(:), ferr_got(:), berr_got(:)
    integer, allocatable :: iwork(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtprfs', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_matrix_quad(n, n, A, seed = 154001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + real(2*n, ep); end do
        call gen_matrix_quad(n, nrhs, B0, seed = 154011 + 47 * i)
        allocate(AP(np), X_ref(n, nrhs), X_got(n, nrhs), work(3*n), iwork(n))
        allocate(ferr_ref(nrhs), berr_ref(nrhs), ferr_got(nrhs), berr_got(nrhs))
        call pack_sym_packed_quad('U', n, A, AP)
        X_ref = B0; X_got = B0
        call dtptrs('U', 'N', 'N', n, nrhs, AP, X_ref, n, info)
        call target_dtptrs('U', 'N', 'N', n, nrhs, AP, X_got, n, info)
        call dtprfs('U', 'N', 'N', n, nrhs, AP, B0, n, X_ref, n, &
                    ferr_ref, berr_ref, work, iwork, info)
        call target_dtprfs('U', 'N', 'N', n, nrhs, AP, B0, n, X_ref, n, &
                           ferr_got, berr_got, info)
        err = maxval(abs(berr_got - berr_ref))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AP, B0, X_ref, X_got, work, iwork, ferr_ref, berr_ref, ferr_got, berr_got)
    end do
    call report_finalize()
end program test_dtprfs
