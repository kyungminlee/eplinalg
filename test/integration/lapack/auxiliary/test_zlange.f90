program test_zlange
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zlange
    use ref_quad_lapack, only: zlange
    implicit none

    integer, parameter :: ms(*) = [8, 32, 96]
    integer, parameter :: ns(*) = [6, 20, 64]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    integer :: i, k, m, n
    complex(ep), allocatable :: A(:,:)
    real(ep),    allocatable :: work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=48) :: label

    call report_init('zlange', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A, seed = 13101 + 83 * i)
        allocate(work(max(1, m)))
        do k = 1, size(norms)
            ref_val = zlange(norms(k), m, n, A, m, work)
            got_val = target_zlange(norms(k), m, n, A, m)
            err = rel_err_scalar(got_val, ref_val)
            tol = 16.0_ep * real(m * n, ep) * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'norm=', norms(k), ',m=', m, ',n=', n
            call report_case(trim(label), err, tol)
        end do
        deallocate(work)
    end do
    call report_finalize()
end program test_zlange
