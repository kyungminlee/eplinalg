program test_zahemv
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_vec
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_zahemv
    use test_data,            only: gen_matrix_complex, gen_vector_complex, gen_vector_quad
    use target_ptzblas,       only: target_name, target_eps, target_zahemv
    implicit none

    integer, parameter :: n = 32
    character(len=1), parameter :: uplo_set(*) = ['L', 'U']
    real(ep),         parameter :: alpha_set(*) = [0.8_ep, -0.7_ep]
    integer :: i, k
    character :: uplo
    complex(ep), allocatable :: A(:,:), x(:)
    real(ep), allocatable :: y_got(:), y_ref(:)
    real(ep) :: alpha, beta, err, tol
    character(len=64) :: label

    call report_init('zahemv', target_name, 0)
    do i = 1, size(uplo_set)
        uplo = uplo_set(i)
        do k = 1, size(alpha_set)
            call gen_matrix_complex(n, n, A,    seed = 3800 + 3 * i + 11 * k)
            call gen_vector_complex(n, x,        seed = 3900 + 5 * i + 13 * k)
            call gen_vector_quad   (n, y_got,    seed = 4000 + 7 * i + 17 * k)
            allocate(y_ref(n)); y_ref = y_got
            alpha = alpha_set(k); beta = -0.4_ep
            call target_zahemv(uplo, n, alpha, A, n, x, 1, beta, y_got, 1)
            call ref_zahemv(uplo, n, alpha, A, x, beta, y_ref)
            err = max_rel_err_vec(y_got, y_ref)
            tol = 64.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,f5.2)') 'uplo=', uplo, ' alpha=', alpha
            call report_case(trim(label), err, tol)
            deallocate(A, x, y_got, y_ref)
        end do
    end do
    call report_finalize()
end program test_zahemv
