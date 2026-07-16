program test_zpotrf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hpd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zpotrf
    use ref_quad_lapack, only: zpotrf
    implicit none

    integer, parameter :: ns(*) = [8, 32, 96]
    integer :: i, j, k, n, info_ref, info_got
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: ju

    call report_init('zpotrf', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hpd_matrix_quad(n, A0, seed = 4101 + 43 * i)
        do ju = 1, size(uplos)
            allocate(A_ref(n, n), A_got(n, n))
            A_ref = A0; A_got = A0
            call zpotrf(uplos(ju), n, A_ref, n, info_ref)
            call target_zpotrf(uplos(ju), n, A_got, n, info_got)
            ! Zero out untouched triangle on both sides — same approach as test_dpotrf.
            if (uplos(ju) == 'U') then
                do j = 1, n
                    do k = j+1, n
                        A_ref(k, j) = (0.0_ep, 0.0_ep)
                        A_got(k, j) = (0.0_ep, 0.0_ep)
                    end do
                end do
            else
                do j = 1, n
                    do k = 1, j-1
                        A_ref(k, j) = (0.0_ep, 0.0_ep)
                        A_got(k, j) = (0.0_ep, 0.0_ep)
                    end do
                end do
            end if
            err = max_rel_err_mat_z(A_got, A_ref)
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got)
        end do
    end do
    call report_finalize()
end program test_zpotrf
