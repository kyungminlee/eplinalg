program test_dsprfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_symmetric_matrix_quad, gen_matrix_quad, pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, &
                                target_dsptrf, target_dsptrs, target_dsprfs
    use ref_quad_lapack, only: dsptrf, dsptrs, dsprfs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs = 3
    integer :: i, n, info, np
    integer, allocatable :: ipiv_ref(:), ipiv_got(:), iwork(:)
    real(ep), allocatable :: A0(:,:), AP0(:), AFP_ref(:), AFP_got(:)
    real(ep), allocatable :: B0(:,:), X_ref(:,:), X_got(:,:), work(:), ferr(:), berr(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dsprfs', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_symmetric_matrix_quad(n, A0, seed = 148001 + 47 * i)
        call gen_matrix_quad(n, nrhs, B0, seed = 148011 + 47 * i)
        allocate(AP0(np), AFP_ref(np), AFP_got(np), ipiv_ref(n), ipiv_got(n))
        allocate(X_ref(n, nrhs), X_got(n, nrhs), work(3*n), iwork(n), ferr(nrhs), berr(nrhs))
        call pack_sym_packed_quad('U', n, A0, AP0)
        AFP_ref = AP0; AFP_got = AP0
        call dsptrf('U', n, AFP_ref, ipiv_ref, info)
        call target_dsptrf('U', n, AFP_got, ipiv_got, info)
        X_ref = B0; X_got = B0
        call dsptrs('U', n, nrhs, AFP_ref, ipiv_ref, X_ref, n, info)
        call target_dsptrs('U', n, nrhs, AFP_got, ipiv_got, X_got, n, info)
        call dsprfs('U', n, nrhs, AP0, AFP_ref, ipiv_ref, B0, n, X_ref, n, &
                    ferr, berr, work, iwork, info)
        call target_dsprfs('U', n, nrhs, AP0, AFP_got, ipiv_got, B0, n, X_got, n, &
                           ferr, berr, info)
        err = max_rel_err_mat(X_got, X_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(AP0, AFP_ref, AFP_got, ipiv_ref, ipiv_got, X_ref, X_got, work, iwork, ferr, berr)
    end do
    call report_finalize()
end program test_dsprfs
