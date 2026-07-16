program test_zgemmtr
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat_z
    use test_data,     only: gen_matrix_complex
    use target_blas,   only: target_name, target_eps, target_zgemmtr
    use ref_quad_blas, only: zgemmtr
    implicit none

    ! Sweep UPLO × TRANSA × TRANSB = 2 × 3 × 3 = 18 combos. 'C'
    ! (conjugate-transpose) is its own code path in zgemmtr and the
    ! most bug-prone in complex BLAS.
    integer, parameter :: cases(*)             = [16, 64, 32, 16, 64, 32, 16, 64, 32, 16, 64, 32, 16, 64, 32, 16, 64, 32]
    integer, parameter :: ks(*)                = [10, 32, 24, 10, 32, 24, 10, 32, 24, 10, 32, 24, 10, 32, 24, 10, 32, 24]
    character(len=1), parameter :: uplos(*)    = ['U', 'U', 'U', 'U', 'U', 'U', 'U', 'U', 'U', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L']
    character(len=1), parameter :: transas(*)  = ['N', 'N', 'N', 'T', 'T', 'T', 'C', 'C', 'C', 'N', 'N', 'N', 'T', 'T', 'T', 'C', 'C', 'C']
    character(len=1), parameter :: transbs(*)  = ['N', 'T', 'C', 'N', 'T', 'C', 'N', 'T', 'C', 'N', 'T', 'C', 'N', 'T', 'C', 'N', 'T', 'C']
    integer :: i, n, k, lda, ldb
    complex(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('zgemmtr', target_name)
    do i = 1, size(cases)
        n = cases(i); k = ks(i)
        if (transas(i) == 'N') then
            lda = n; call gen_matrix_complex(n, k, A, seed = 2801 + 17 * i)
        else
            lda = k; call gen_matrix_complex(k, n, A, seed = 2801 + 17 * i)
        end if
        if (transbs(i) == 'N') then
            ldb = k; call gen_matrix_complex(k, n, B, seed = 2811 + 17 * i)
        else
            ldb = n; call gen_matrix_complex(n, k, B, seed = 2811 + 17 * i)
        end if
        call gen_matrix_complex(n, n, C0, seed = 2821 + 17 * i)
        alpha = cmplx(0.6_ep, 0.2_ep, ep)
        beta  = cmplx(0.4_ep, -0.1_ep, ep)
        allocate(C_ref(n,n), C_got(n,n))
        C_ref = C0; C_got = C0
        call zgemmtr(uplos(i), transas(i), transbs(i), n, k, alpha, A, lda, &
                     B, ldb, beta, C_ref, n)
        call target_zgemmtr(uplos(i), transas(i), transbs(i), n, k, alpha, A, lda, &
                            B, ldb, beta, C_got, n)
        err = max_rel_err_mat_z(C_got, C_ref)
        tol = 64.0_ep * real(k, ep) * target_eps
        write(label, '(a,a,a,a,a,a,a,i0,a,i0)') 'uplo=', uplos(i), &
            ',ta=', transas(i), ',tb=', transbs(i), ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(A, B, C_ref, C_got)
    end do
    call report_finalize()
end program test_zgemmtr
