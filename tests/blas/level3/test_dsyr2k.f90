program test_dsyr2k
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad
    use target_blas,   only: target_name, target_eps, target_dsyr2k
    use ref_quad_blas, only: dsyr2k
    implicit none

    integer, parameter :: ns(*) = [4, 32, 100, 64]
    integer, parameter :: ks(*) = [5, 40, 80, 48]
    ! Sweep UPLO × TRANS = 4 combos. TRANS='T' transposes A and B
    ! (reads them as k×n) — independent code from TRANS='N'.
    character(len=1), parameter :: uplos(*)   = ['U', 'U', 'L', 'L']
    character(len=1), parameter :: transes(*) = ['N', 'T', 'N', 'T']
    integer :: i, n, k, lda
    real(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call report_init('dsyr2k', target_name)
    do i = 1, size(ns)
        n = ns(i); k = ks(i)
        if (transes(i) == 'N') then
            call gen_matrix_quad(n, k, A,  seed = 881 + 23 * i)
            call gen_matrix_quad(n, k, B,  seed = 891 + 23 * i)
            lda = n
        else
            call gen_matrix_quad(k, n, A,  seed = 881 + 23 * i)
            call gen_matrix_quad(k, n, B,  seed = 891 + 23 * i)
            lda = k
        end if
        call gen_matrix_quad(n, n, C0, seed = 901 + 23 * i)
        alpha = real(0.6_ep, ep)
        beta  = real(0.4_ep, ep)
        allocate(C_ref(n, n), C_got(n, n))
        C_ref = C0
        C_got = C0
        call dsyr2k(uplos(i), transes(i), n, k, alpha, A, lda, B, lda, beta, C_ref, n)
        call target_dsyr2k(uplos(i), transes(i), n, k, alpha, A, lda, B, lda, beta, C_got, n)
        err = max_rel_err_mat(C_got, C_ref)
        tol = 16.0_ep * 4.0_ep * real(k, ep) * target_eps
        write(label, '(a,a,a,a,a,i0,a,i0)') 'uplo=', uplos(i), &
            ',trans=', transes(i), ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(A, B, C_ref, C_got)
    end do
    call report_finalize()
end program test_dsyr2k
