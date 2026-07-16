! zlanht: norm of a complex Hermitian tridiagonal matrix.
program test_zlanht
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_vector_quad, gen_vector_complex
    use target_lapack,   only: target_name, target_eps, target_zlanht
    use ref_quad_lapack, only: zlanht
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    integer :: i, k, n
    real(ep), allocatable :: d(:)
    complex(ep), allocatable :: e(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=48) :: label

    call report_init('zlanht', target_name)
    do i = 1, size(ns)
        n = ns(i)
        allocate(d(n), e(n-1))
        call gen_vector_quad(n,   d, seed = 19101 + 47 * i)
        call gen_vector_complex(n-1, e, seed = 19111 + 47 * i)
        do k = 1, size(norms)
            ref_val = zlanht(norms(k), n, d, e)
            got_val = target_zlanht(norms(k), n, d, e)
            err = rel_err_scalar(got_val, ref_val)
            tol = 16.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,i0)') 'norm=', norms(k), ',n=', n
            call report_case(trim(label), err, tol)
        end do
        deallocate(d, e)
    end do
    call report_finalize()
end program test_zlanht
