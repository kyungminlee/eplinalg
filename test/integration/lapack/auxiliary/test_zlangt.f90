! zlangt: norm of a complex tridiagonal matrix.
program test_zlangt
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_vector_complex
    use target_lapack,   only: target_name, target_eps, target_zlangt
    use ref_quad_lapack, only: zlangt
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    integer :: i, k, n
    complex(ep), allocatable :: dl(:), d(:), du(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=48) :: label

    call report_init('zlangt', target_name)
    do i = 1, size(ns)
        n = ns(i)
        allocate(dl(n-1), d(n), du(n-1))
        call gen_vector_complex(n-1, dl, seed = 18301 + 47 * i)
        call gen_vector_complex(n,   d,  seed = 18311 + 47 * i)
        call gen_vector_complex(n-1, du, seed = 18321 + 47 * i)
        do k = 1, size(norms)
            ref_val = zlangt(norms(k), n, dl, d, du)
            got_val = target_zlangt(norms(k), n, dl, d, du)
            err = rel_err_scalar(got_val, ref_val)
            tol = 16.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,i0)') 'norm=', norms(k), ',n=', n
            call report_case(trim(label), err, tol)
        end do
        deallocate(dl, d, du)
    end do
    call report_finalize()
end program test_zlangt
