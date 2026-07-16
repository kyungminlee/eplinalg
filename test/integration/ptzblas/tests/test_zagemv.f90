program test_zagemv
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_vec
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_xagemv
    use test_data,            only: gen_matrix_complex, gen_vector_complex, gen_vector_quad
    use target_ptzblas,       only: target_name, target_eps, target_zagemv
    implicit none

    integer, parameter :: m = 32, n = 24
    character(len=1), parameter :: trans_set(*) = ['N', 'T']
    integer :: i, ny
    character :: trans
    complex(ep), allocatable :: A(:,:), x(:)
    real(ep), allocatable :: y_got(:), y_ref(:)
    real(ep) :: alpha, beta, err, tol
    character(len=32) :: label

    call report_init('zagemv', target_name, 0)
    do i = 1, size(trans_set)
        trans = trans_set(i)
        call gen_matrix_complex(m, n, A, seed = 3500 + 3 * i)
        if (trans == 'N') then
            call gen_vector_complex(n, x,    seed = 3600 + 5 * i)
            call gen_vector_quad   (m, y_got, seed = 3700 + 7 * i)
            ny = m
        else
            call gen_vector_complex(m, x,    seed = 3600 + 5 * i)
            call gen_vector_quad   (n, y_got, seed = 3700 + 7 * i)
            ny = n
        end if
        allocate(y_ref(ny)); y_ref = y_got
        alpha = 1.5_ep; beta = -0.3_ep
        call target_zagemv(trans, m, n, alpha, A, m, x, 1, beta, y_got, 1)
        call ref_xagemv(trans, m, n, alpha, A, x, beta, y_ref)
        err = max_rel_err_vec(y_got, y_ref)
        tol = 64.0_ep * real(max(m, n), ep) * target_eps
        write(label, '(a,a)') 'trans=', trans
        call report_case(trim(label), err, tol)
        deallocate(A, x, y_got, y_ref)
    end do
    call report_finalize()
end program test_zagemv
