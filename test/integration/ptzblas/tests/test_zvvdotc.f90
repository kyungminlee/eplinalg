program test_zvvdotc
    use prec_kinds,            only: ep
    use compare,               only: rel_err_scalar_z
    use ptzblas_prec_report,   only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_blas, only: zdotc
    use test_data,             only: gen_vector_complex
    use target_ptzblas,        only: target_name, target_eps, target_zvvdotc
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer :: i, n
    complex(ep), allocatable :: x(:), y(:)
    complex(ep) :: ref, got
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('zvvdotc', target_name, 0)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x, seed = 2200 + 7 * i)
        call gen_vector_complex(n, y, seed = 2300 + 11 * i)
        got = (0.0_ep, 0.0_ep)
        call target_zvvdotc(n, got, x, 1, y, 1)
        ref = zdotc(n, x, 1, y, 1)
        err = rel_err_scalar_z(got, ref)
        tol = 32.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x, y)
    end do
    call report_finalize()
end program test_zvvdotc
