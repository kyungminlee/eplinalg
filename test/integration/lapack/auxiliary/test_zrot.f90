! zrot: complex Givens rotation (real c, complex s).
program test_zrot
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec_z
    use test_data,       only: gen_vector_complex
    use target_lapack,   only: target_name, target_eps, target_zrot
    use ref_quad_lapack, only: zrot
    implicit none

    integer, parameter :: ns(*) = [4, 32, 128]
    integer :: i, n, incx, incy
    complex(ep), allocatable :: cx0(:), cy0(:), cx_ref(:), cy_ref(:), cx_got(:), cy_got(:)
    real(ep)    :: c, err1, err2, tol
    complex(ep) :: s
    character(len=48) :: label

    c = 0.6_ep
    s = (0.8_ep, 0.0_ep)
    call report_init('zrot', target_name)
    do i = 1, size(ns)
        n = ns(i)
        do incx = 1, 2
            do incy = 1, 2
                call gen_vector_complex(1 + (n - 1) * incx, cx0, seed = 19301 + 89 * i + incx)
                call gen_vector_complex(1 + (n - 1) * incy, cy0, seed = 19401 + 97 * i + incy)
                allocate(cx_ref(size(cx0)), cy_ref(size(cy0)), cx_got(size(cx0)), cy_got(size(cy0)))
                cx_ref = cx0; cy_ref = cy0; cx_got = cx0; cy_got = cy0
                call zrot(n, cx_ref, incx, cy_ref, incy, c, s)
                call target_zrot(n, cx_got, incx, cy_got, incy, c, s)
                err1 = max_rel_err_vec_z(cx_got, cx_ref)
                err2 = max_rel_err_vec_z(cy_got, cy_ref)
                tol  = 16.0_ep * real(n, ep) * target_eps
                write(label, '(a,i0,a,i0,a,i0)') 'n=', n, ',ix=', incx, ',iy=', incy
                call report_case(trim(label), max(err1, err2), tol)
                deallocate(cx0, cy0, cx_ref, cy_ref, cx_got, cy_got)
            end do
        end do
    end do
    call report_finalize()
end program test_zrot
