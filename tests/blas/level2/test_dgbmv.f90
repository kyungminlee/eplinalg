program test_dgbmv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_matrix_quad, gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dgbmv
    use ref_quad_blas, only: dgbmv
    implicit none

    integer, parameter :: cases(*) = [20, 100, 20, 100]
    ! Cycle TRANS so the transpose path is exercised alongside 'N'.
    character(len=1), parameter :: transes(*) = ['N', 'N', 'T', 'T']
    integer :: i, n, kl, ku, lda
    real(ep), allocatable :: A(:,:), x(:), y0(:), y_ref(:), y_got(:)
    real(ep) :: alpha, beta, err, tol
    character(len=32) :: label

    call report_init('dgbmv', target_name)
    do i = 1, size(cases)
        n  = cases(i)
        kl = 2; ku = 3
        lda = kl + ku + 1
        call gen_matrix_quad(lda, n, A, seed = 431 + 17 * i)
        call gen_vector_quad(n, x,  seed = 441 + 17 * i)
        call gen_vector_quad(n, y0, seed = 451 + 17 * i)
        alpha = real(0.5_ep, ep)
        beta  = real(0.7_ep, ep)
        allocate(y_ref(n), y_got(n))
        y_ref = y0
        y_got = y0
        call dgbmv(transes(i), n, n, kl, ku, alpha, A, lda, x, 1, beta, y_ref, 1)
        call target_dgbmv(transes(i), n, n, kl, ku, alpha, A, lda, x, 1, beta, y_got, 1)
        err = max_rel_err_vec(y_got, y_ref)
        tol = 16.0_ep * 2.0_ep * real(lda, ep) * target_eps
        write(label, '(a,a,a,i0)') 'trans=', transes(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(y_ref, y_got)
    end do
    call report_finalize()
end program test_dgbmv
