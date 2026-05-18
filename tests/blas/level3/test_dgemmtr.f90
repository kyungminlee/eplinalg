program test_dgemmtr
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad
    use target_blas,   only: target_name, target_eps, target_dgemmtr
    use ref_quad_blas, only: dgemmtr
    implicit none

    ! Sweep UPLO × TRANSA × TRANSB = 8 combos.
    integer, parameter :: cases(*)             = [16, 64, 32, 48, 16, 64, 32, 48]
    integer, parameter :: ks(*)                = [10, 32, 24, 20, 10, 32, 24, 20]
    character(len=1), parameter :: uplos(*)    = ['U', 'U', 'U', 'U', 'L', 'L', 'L', 'L']
    character(len=1), parameter :: transas(*)  = ['N', 'N', 'T', 'T', 'N', 'N', 'T', 'T']
    character(len=1), parameter :: transbs(*)  = ['N', 'T', 'N', 'T', 'N', 'T', 'N', 'T']
    integer :: i, n, k, lda, ldb
    real(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=64) :: label

    call report_init('dgemmtr', target_name)
    do i = 1, size(cases)
        n = cases(i); k = ks(i)
        if (transas(i) == 'N') then
            lda = n; call gen_matrix_quad(n, k, A, seed = 2701 + 17 * i)
        else
            lda = k; call gen_matrix_quad(k, n, A, seed = 2701 + 17 * i)
        end if
        if (transbs(i) == 'N') then
            ldb = k; call gen_matrix_quad(k, n, B, seed = 2711 + 17 * i)
        else
            ldb = n; call gen_matrix_quad(n, k, B, seed = 2711 + 17 * i)
        end if
        call gen_matrix_quad(n, n, C0, seed = 2721 + 17 * i)
        alpha = 0.6_ep; beta = 0.4_ep
        allocate(C_ref(n,n), C_got(n,n))
        C_ref = C0; C_got = C0
        call dgemmtr(uplos(i), transas(i), transbs(i), n, k, alpha, A, lda, &
                     B, ldb, beta, C_ref, n)
        call target_dgemmtr(uplos(i), transas(i), transbs(i), n, k, alpha, A, lda, &
                            B, ldb, beta, C_got, n)
        err = max_rel_err_mat(C_got, C_ref)
        tol = 64.0_ep * real(k, ep) * target_eps
        write(label, '(a,a,a,a,a,a,a,i0,a,i0)') 'uplo=', uplos(i), &
            ',ta=', transas(i), ',tb=', transbs(i), ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(A, B, C_ref, C_got)
    end do
    call report_finalize()
end program test_dgemmtr
