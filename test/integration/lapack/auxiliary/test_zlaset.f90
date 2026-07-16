program test_zlaset
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use target_lapack,   only: target_name, target_eps, target_zlaset
    use ref_quad_lapack, only: zlaset
    implicit none

    integer, parameter :: ms(*) = [8, 32, 96]
    integer, parameter :: ns(*) = [6, 20, 64]
    character(len=1), parameter :: uplos(3) = ['U', 'L', 'A']
    integer :: i, k, m, n
    complex(ep), allocatable :: A_ref(:,:), A_got(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zlaset', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        alpha = cmplx( 0.375_ep,  0.125_ep, ep)
        beta  = cmplx(-1.625_ep, -0.500_ep, ep)
        do k = 1, size(uplos)
            allocate(A_ref(m, n), A_got(m, n))
            A_ref = (1.0_ep, 0.5_ep); A_got = (1.0_ep, 0.5_ep)
            call zlaset(uplos(k), m, n, alpha, beta, A_ref, m)
            call target_zlaset(uplos(k), m, n, alpha, beta, A_got, m)
            err = max_rel_err_mat_z(A_got, A_ref)
            tol = 16.0_ep * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'uplo=', uplos(k), ',m=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got)
        end do
    end do
    call report_finalize()
end program test_zlaset
