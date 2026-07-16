program test_zsyr
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat_z
    use test_data,     only: gen_complex_symmetric_quad, gen_vector_complex
    use target_lapack, only: target_name, target_eps, target_zsyr
    use ref_quad_lapack, only: zsyr
    implicit none
    integer, parameter :: ns(*) = [10, 50, 100]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, ju
    complex(ep), allocatable :: A0(:,:), x(:), A_ref(:,:), A_got(:,:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('zsyr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_complex_symmetric_quad(n, A0, seed = 70101 + 19 * i)
        call gen_vector_complex(n, x, seed = 70111 + 19 * i)
        alpha = cmplx(0.6_ep, 0.2_ep, ep)
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n)); A_ref = A0; A_got = A0
            call zsyr(uplos(ju), n, alpha, x, 1, A_ref, n)
            call target_zsyr(uplos(ju), n, alpha, x, 1, A_got, n)
            err = max_rel_err_mat_z(A_got, A_ref)
            tol = 32.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got)
        end do
    end do
    call report_finalize()
end program test_zsyr
