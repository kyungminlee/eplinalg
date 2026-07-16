! dlantp: norm of a real triangular packed matrix.
program test_dlantp
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_dlantp
    use ref_quad_lapack, only: dlantp
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    character(len=1), parameter :: diags(2) = ['N', 'U']
    integer :: i, k, u, dg, n
    real(ep), allocatable :: AP(:), work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=64) :: label

    call report_init('dlantp', target_name)
    do i = 1, size(ns)
        n = ns(i)
        allocate(AP(n*(n+1)/2))
        call gen_vector_quad(n*(n+1)/2, AP, seed = 19501 + 47 * i)
        allocate(work(max(1, n)))
        do u = 1, size(uplos)
            do dg = 1, size(diags)
                do k = 1, size(norms)
                    ref_val = dlantp(norms(k), uplos(u), diags(dg), n, AP, work)
                    got_val = target_dlantp(norms(k), uplos(u), diags(dg), n, AP)
                    err = rel_err_scalar(got_val, ref_val)
                    tol = 16.0_ep * real(n, ep) * target_eps
                    write(label, '(a,a,a,a,a,a,a,i0)') &
                        'norm=', norms(k), ',uplo=', uplos(u), ',diag=', diags(dg), ',n=', n
                    call report_case(trim(label), err, tol)
                end do
            end do
        end do
        deallocate(AP, work)
    end do
    call report_finalize()
end program test_dlantp
