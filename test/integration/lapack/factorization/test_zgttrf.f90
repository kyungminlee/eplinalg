program test_zgttrf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec_z
    use test_data,       only: gen_vector_complex
    use target_lapack,   only: target_name, target_eps, target_zgttrf
    use ref_quad_lapack, only: zgttrf
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer :: i, n, info, j
    complex(ep), allocatable :: dl0(:), d0(:), du0(:)
    complex(ep), allocatable :: dl_ref(:), d_ref(:), du_ref(:), du2_ref(:)
    complex(ep), allocatable :: dl_got(:), d_got(:), du_got(:), du2_got(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgttrf', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_complex(n-1, dl0, seed = 98001 + 47 * i)
        call gen_vector_complex(n,   d0,  seed = 98011 + 47 * i)
        call gen_vector_complex(n-1, du0, seed = 98021 + 47 * i)
        do j = 1, n
            d0(j) = d0(j) + cmplx(real(4, ep), 0.0_ep, ep)
        end do
        allocate(dl_ref(n-1), d_ref(n), du_ref(n-1), du2_ref(n-2))
        allocate(dl_got(n-1), d_got(n), du_got(n-1), du2_got(n-2))
        allocate(ipiv_ref(n), ipiv_got(n))
        dl_ref = dl0; d_ref = d0; du_ref = du0
        dl_got = dl0; d_got = d0; du_got = du0
        call zgttrf(n, dl_ref, d_ref, du_ref, du2_ref, ipiv_ref, info)
        call target_zgttrf(n, dl_got, d_got, du_got, du2_got, ipiv_got, info)
        err = max(max_rel_err_vec_z(d_got, d_ref), max_rel_err_vec_z(du_got, du_ref))
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(dl0, d0, du0, dl_ref, d_ref, du_ref, du2_ref)
        deallocate(dl_got, d_got, du_got, du2_got, ipiv_ref, ipiv_got)
    end do
    call report_finalize()
end program test_zgttrf
