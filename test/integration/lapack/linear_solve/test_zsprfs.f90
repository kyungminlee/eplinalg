program test_zsprfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_complex_symmetric_quad, gen_matrix_complex, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, &
                                target_zsptrf, target_zsptrs, target_zsprfs
    use ref_quad_lapack, only: zsptrf, zsptrs, zsprfs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs = 3
    integer :: i, n, info, np
    integer, allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep), allocatable :: A0(:,:), AP0(:), AFP_ref(:), AFP_got(:)
    complex(ep), allocatable :: B0(:,:), X_ref(:,:), X_got(:,:), work(:)
    real(ep), allocatable :: rwork(:), ferr(:), berr(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zsprfs', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_complex_symmetric_quad(n, A0, seed = 149001 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 149011 + 47 * i)
        allocate(AP0(np), AFP_ref(np), AFP_got(np), ipiv_ref(n), ipiv_got(n))
        allocate(X_ref(n, nrhs), X_got(n, nrhs), work(2*n), rwork(n), ferr(nrhs), berr(nrhs))
        call pack_herm_packed_quad('U', n, A0, AP0)
        AFP_ref = AP0; AFP_got = AP0
        call zsptrf('U', n, AFP_ref, ipiv_ref, info)
        call target_zsptrf('U', n, AFP_got, ipiv_got, info)
        X_ref = B0; X_got = B0
        call zsptrs('U', n, nrhs, AFP_ref, ipiv_ref, X_ref, n, info)
        call target_zsptrs('U', n, nrhs, AFP_got, ipiv_got, X_got, n, info)
        call zsprfs('U', n, nrhs, AP0, AFP_ref, ipiv_ref, B0, n, X_ref, n, &
                    ferr, berr, work, rwork, info)
        call target_zsprfs('U', n, nrhs, AP0, AFP_got, ipiv_got, B0, n, X_got, n, &
                           ferr, berr, info)
        err = max_rel_err_mat_z(X_got, X_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(AP0, AFP_ref, AFP_got, ipiv_ref, ipiv_got, X_ref, X_got, work, rwork, ferr, berr)
    end do
    call report_finalize()
end program test_zsprfs
