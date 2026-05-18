program test_zher2k
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat_z
    use test_data,     only: gen_matrix_complex
    use target_blas,   only: target_name, target_eps, target_zher2k
    use ref_quad_blas, only: zher2k
    implicit none

    integer, parameter :: cases(*)              = [16, 64, 32, 48]
    integer, parameter :: ks(*)                 = [10, 32, 24, 20]
    ! Sweep UPLO × TRANS = 4 combos.
    character(len=1), parameter :: uplos(*)    = ['U', 'U', 'L', 'L']
    character(len=1), parameter :: transes(*)  = ['N', 'C', 'N', 'C']
    integer :: i, n, k, lda
    complex(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    complex(ep) :: alpha
    real(ep) :: beta, err, tol
    character(len=64) :: label

    call report_init('zher2k', target_name)
    do i = 1, size(cases)
        n = cases(i); k = ks(i)
        if (transes(i) == 'N') then
            lda = n
            call gen_matrix_complex(n, k, A, seed = 3201 + 17 * i)
            call gen_matrix_complex(n, k, B, seed = 3211 + 17 * i)
        else
            lda = k
            call gen_matrix_complex(k, n, A, seed = 3201 + 17 * i)
            call gen_matrix_complex(k, n, B, seed = 3211 + 17 * i)
        end if
        call gen_matrix_complex(n, n, C0, seed = 3221 + 17 * i)
        alpha = cmplx(0.6_ep, 0.2_ep, ep)
        beta  = 0.4_ep
        allocate(C_ref(n,n), C_got(n,n))
        C_ref = C0; C_got = C0
        call zher2k(uplos(i), transes(i), n, k, alpha, A, lda, B, lda, beta, C_ref, n)
        call target_zher2k(uplos(i), transes(i), n, k, alpha, A, lda, B, lda, beta, C_got, n)
        err = max_rel_err_mat_z(C_got, C_ref)
        tol = 64.0_ep * real(k, ep) * target_eps
        write(label, '(a,a,a,a,a,i0,a,i0)') 'uplo=', uplos(i), &
            ',trans=', transes(i), ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(A, B, C_ref, C_got)
    end do
    call report_finalize()
end program test_zher2k
