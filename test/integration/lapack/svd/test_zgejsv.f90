! zgejsv: complex Jacobi SVD with optional preconditioning. Smoke test.
program test_zgejsv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgejsv
    use ref_quad_lapack, only: zgejsv
    implicit none

    integer, parameter :: ms(*) = [16, 24]
    integer, parameter :: ns(*) = [8, 16]
    integer :: i, m, n, info_r, info_g, lwork, lrwork
    complex(ep), allocatable :: A0(:,:), A_r(:,:), A_g(:,:)
    complex(ep), allocatable :: U_r(:,:), V_r(:,:), U_g(:,:), V_g(:,:), cwork(:)
    real(ep), allocatable :: sva_r(:), sva_g(:), rwork(:)
    integer, allocatable :: iwork(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgejsv', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A0, seed = 154101 + 47 * i)
        allocate(A_r(m,n), A_g(m,n), U_r(m,n), V_r(n,n), U_g(m,n), V_g(n,n))
        allocate(sva_r(n), sva_g(n))
        lwork  = max(2, 5*n + 2*n*n)
        lrwork = max(7, n + 2*m)
        allocate(cwork(lwork), rwork(lrwork), iwork(max(1, m+3*n)))
        A_r = A0; A_g = A0
        call zgejsv('G', 'U', 'V', 'N', 'N', 'N', m, n, A_r, m, sva_r, &
                    U_r, m, V_r, n, cwork, lwork, rwork, lrwork, iwork, info_r)
        deallocate(cwork, rwork, iwork)
        call target_zgejsv('G', 'U', 'V', 'N', 'N', 'N', m, n, A_g, m, sva_g, &
                           U_g, m, V_g, n, info_g)
        err = max_rel_err_vec(sva_g, sva_r)
        tol = 256.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_r, A_g, U_r, V_r, U_g, V_g, sva_r, sva_g)
    end do
    call report_finalize()
end program test_zgejsv
