! zrscl: robust reciprocal scaling x <- x / a (complex divisor).
program test_zrscl
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec_z
    use test_data,       only: gen_vector_complex
    use target_lapack,   only: target_name, target_eps, target_zrscl
    use ref_quad_lapack, only: zrscl
    implicit none

    integer, parameter :: ns(*) = [5, 32, 128]
    complex(ep), parameter :: scales(*) = [(1.0_ep, 0.0_ep), (1.0_ep, 1.0_ep), &
                                           (3.0_ep, -2.0_ep), (1.0e8_ep, -1.0e-3_ep)]
    integer :: i, j, n, incx
    complex(ep), allocatable :: x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zrscl', target_name)
    do i = 1, size(ns)
        n = ns(i)
        do incx = 1, 2
            call gen_vector_complex(1 + (n - 1) * incx, x0, seed = 19201 + 79 * i + incx)
            do j = 1, size(scales)
                allocate(x_ref(size(x0)), x_got(size(x0)))
                x_ref = x0; x_got = x0
                call zrscl(n, scales(j), x_ref, incx)
                call target_zrscl(n, scales(j), x_got, incx)
                err = max_rel_err_vec_z(x_got, x_ref)
                tol = 16.0_ep * real(n, ep) * target_eps
                write(label, '(a,i0,a,i0,a,i0)') 'n=', n, ',incx=', incx, ',sc=', j
                call report_case(trim(label), err, tol)
                deallocate(x_ref, x_got)
            end do
            deallocate(x0)
        end do
    end do
    call report_finalize()
end program test_zrscl
