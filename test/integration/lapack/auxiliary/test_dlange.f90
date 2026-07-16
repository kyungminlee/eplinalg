program test_dlange
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dlange
    use ref_quad_lapack, only: dlange
    implicit none

    integer, parameter :: ms(*) = [8, 32, 96]
    integer, parameter :: ns(*) = [6, 20, 64]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    integer :: i, k, m, n
    real(ep), allocatable :: A(:,:), work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=48) :: label

    call report_init('dlange', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A, seed = 13001 + 83 * i)
        allocate(work(max(1, m)))

        do k = 1, size(norms)
            ref_val = dlange(norms(k), m, n, A, m, work)
            got_val = target_dlange(norms(k), m, n, A, m)
            err = rel_err_scalar(got_val, ref_val)
            tol = 16.0_ep * real(m * n, ep) * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'norm=', norms(k), ',m=', m, ',n=', n
            call report_case(trim(label), err, tol)
        end do

        deallocate(work)
    end do
    call report_finalize()
end program test_dlange
