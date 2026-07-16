! ztzrzf: RZ factorization of a trapezoidal matrix (complex).
program test_ztzrzf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztzrzf
    use ref_quad_lapack, only: ztzrzf
    implicit none

    integer, parameter :: ms(*) = [6, 16]
    integer, parameter :: ns(*) = [12, 32]
    integer :: i, m, n, info, lwork, j
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: tau_ref(:), tau_got(:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztzrzf', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A0, seed = 20151 + 89 * i)
        do j = 1, m-1
            A0(j+1:m, j) = (0.0_ep, 0.0_ep)
        end do
        allocate(A_ref(m,n), A_got(m,n), tau_ref(m), tau_got(m))
        A_ref = A0; A_got = A0
        call ztzrzf(m, n, A_ref, m, tau_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call ztzrzf(m, n, A_ref, m, tau_ref, work, lwork, info)
        deallocate(work)
        call target_ztzrzf(m, n, A_got, m, tau_got, info)
        do j = 1, m
            A_ref(j+1:m, j) = (0.0_ep, 0.0_ep)
            A_got(j+1:m, j) = (0.0_ep, 0.0_ep)
        end do
        err = max_rel_err_mat_z(A_got(1:m,1:m), A_ref(1:m,1:m))
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_ztzrzf
