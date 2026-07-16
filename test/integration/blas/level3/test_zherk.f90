program test_zherk
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat_z
    use test_data,     only: gen_matrix_complex
    use target_blas,   only: target_name, target_eps, target_zherk
    use ref_quad_blas, only: zherk
    implicit none

    integer, parameter :: ns(*) = [4, 32, 80, 48]
    integer, parameter :: ks(*) = [5, 40, 60, 28]
    ! Sweep UPLO × TRANS = 4 combos. ZHERK's 'C' is conjugate-transpose
    ! (the only valid transpose option for Hermitian) — distinct code
    ! from 'N' and the most likely place for a conjugation slip.
    character(len=1), parameter :: uplos(*)   = ['U', 'U', 'L', 'L']
    character(len=1), parameter :: transes(*) = ['N', 'C', 'N', 'C']
    integer :: i, n, k, lda
    complex(ep), allocatable :: A(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call report_init('zherk', target_name)
    do i = 1, size(ns)
        n = ns(i); k = ks(i)
        if (transes(i) == 'N') then
            call gen_matrix_complex(n, k, A, seed = 1061 + 29 * i); lda = n
        else
            call gen_matrix_complex(k, n, A, seed = 1061 + 29 * i); lda = k
        end if
        call gen_matrix_complex(n, n, C0, seed = 1071 + 29 * i)
        alpha = real(0.6_ep, ep)
        beta  = real(0.4_ep, ep)
        allocate(C_ref(n, n), C_got(n, n))
        C_ref = C0
        C_got = C0
        call zherk(uplos(i), transes(i), n, k, alpha, A, lda, beta, C_ref, n)
        call target_zherk(uplos(i), transes(i), n, k, alpha, A, lda, &
                          beta, C_got, n)
        err = max_rel_err_mat_z(C_got, C_ref)
        tol = 32.0_ep * 8.0_ep * real(k, ep) * target_eps
        write(label, '(a,a,a,a,a,i0,a,i0)') 'uplo=', uplos(i), &
            ',trans=', transes(i), ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(A, C0, C_ref, C_got)
    end do
    call report_finalize()
end program test_zherk
