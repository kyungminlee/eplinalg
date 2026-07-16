program test_zgtcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_vector_complex
    use target_lapack,   only: target_name, target_eps, target_zgttrf, target_zgtcon
    use ref_quad_lapack, only: zgttrf, zgtcon
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer :: i, n, info, j
    complex(ep), allocatable :: dl0(:), d0(:), du0(:)
    complex(ep), allocatable :: dl_ref(:), d_ref(:), du_ref(:), du2_ref(:)
    complex(ep), allocatable :: dl_got(:), d_got(:), du_got(:), du2_got(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep), allocatable :: work(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('zgtcon', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_complex(n-1, dl0, seed = 117001 + 47 * i)
        call gen_vector_complex(n,   d0,  seed = 117011 + 47 * i)
        call gen_vector_complex(n-1, du0, seed = 117021 + 47 * i)
        do j = 1, n; d0(j) = d0(j) + cmplx(real(4, ep), 0.0_ep, ep); end do
        anorm = 0.0_ep
        do j = 1, n
            anorm = max(anorm, &
                merge(abs(dl0(j-1)), 0.0_ep, j > 1) + abs(d0(j)) + &
                merge(abs(du0(j)),   0.0_ep, j < n))
        end do
        allocate(dl_ref(n-1), d_ref(n), du_ref(n-1), du2_ref(n-2))
        allocate(dl_got(n-1), d_got(n), du_got(n-1), du2_got(n-2))
        allocate(ipiv_ref(n), ipiv_got(n), work(2*n))
        dl_ref = dl0; d_ref = d0; du_ref = du0
        dl_got = dl0; d_got = d0; du_got = du0
        call zgttrf(n, dl_ref, d_ref, du_ref, du2_ref, ipiv_ref, info)
        call target_zgttrf(n, dl_got, d_got, du_got, du2_got, ipiv_got, info)
        call zgtcon('1', n, dl_ref, d_ref, du_ref, du2_ref, ipiv_ref, &
                    anorm, rcond_ref, work, info)
        call target_zgtcon('1', n, dl_got, d_got, du_got, du2_got, ipiv_got, &
                           anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 100.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(dl0, d0, du0, dl_ref, d_ref, du_ref, du2_ref)
        deallocate(dl_got, d_got, du_got, du2_got, ipiv_ref, ipiv_got, work)
    end do
    call report_finalize()
end program test_zgtcon
