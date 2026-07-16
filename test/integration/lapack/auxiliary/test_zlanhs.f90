! zlanhs: norm of a complex Hessenberg matrix.
program test_zlanhs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zlanhs
    use ref_quad_lapack, only: zlanhs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    integer :: i, k, n, j
    complex(ep), allocatable :: A(:,:)
    real(ep), allocatable :: work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=48) :: label

    call report_init('zlanhs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 18501 + 47 * i)
        do j = 1, n; A(j+2:n, j) = (0.0_ep, 0.0_ep); end do
        allocate(work(max(1, n)))
        do k = 1, size(norms)
            ref_val = zlanhs(norms(k), n, A, n, work)
            got_val = target_zlanhs(norms(k), n, A, n)
            err = rel_err_scalar(got_val, ref_val)
            tol = 16.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,i0)') 'norm=', norms(k), ',n=', n
            call report_case(trim(label), err, tol)
        end do
        deallocate(A, work)
    end do
    call report_finalize()
end program test_zlanhs
