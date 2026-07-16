! drscl: robust reciprocal scaling x <- x / sa.
program test_drscl
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_drscl
    use ref_quad_lapack, only: drscl
    implicit none

    integer, parameter :: ns(*) = [5, 32, 128]
    real(ep), parameter :: scales(*) = [1.0_ep, 0.5_ep, 1.0e10_ep, 1.0e-10_ep]
    integer :: i, j, n, incx
    real(ep), allocatable :: x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('drscl', target_name)
    do i = 1, size(ns)
        n = ns(i)
        do incx = 1, 2
            call gen_vector_quad(1 + (n - 1) * incx, x0, seed = 19001 + 67 * i + incx)
            do j = 1, size(scales)
                allocate(x_ref(size(x0)), x_got(size(x0)))
                x_ref = x0; x_got = x0
                call drscl(n, scales(j), x_ref, incx)
                call target_drscl(n, scales(j), x_got, incx)
                err = max_rel_err_vec(x_got, x_ref)
                tol = 16.0_ep * real(n, ep) * target_eps
                write(label, '(a,i0,a,i0,a,i0)') 'n=', n, ',incx=', incx, ',sc=', j
                call report_case(trim(label), err, tol)
                deallocate(x_ref, x_got)
            end do
            deallocate(x0)
        end do
    end do
    call report_finalize()
end program test_drscl
