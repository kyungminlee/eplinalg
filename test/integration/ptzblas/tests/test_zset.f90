program test_zset
    use prec_kinds,           only: ep
    use compare,              only: max_rel_err_vec_z
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_zset
    use test_data,            only: gen_vector_complex
    use target_ptzblas,       only: target_name, target_eps, target_zset
    implicit none

    integer, parameter :: cases(*) = [50, 500, 5000]
    integer :: i, n
    complex(ep), allocatable :: x_got(:), x_ref(:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('zset', target_name, 0)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x_got, seed = 2600 + 3 * i)
        x_ref = x_got
        alpha = cmplx(0.6_ep + 0.05_ep * real(i, ep), &
                      -0.2_ep * real(i, ep), ep)
        call target_zset(n, alpha, x_got, 1)
        call ref_zset(n, alpha, x_ref)
        err = max_rel_err_vec_z(x_got, x_ref)
        tol = 4.0_ep * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x_got, x_ref)
    end do
    call report_finalize()
end program test_zset
