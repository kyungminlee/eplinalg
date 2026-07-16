! zlantp: norm of a complex triangular packed matrix.
program test_zlantp
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_vector_complex
    use target_lapack,   only: target_name, target_eps, target_zlantp
    use ref_quad_lapack, only: zlantp
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    character(len=1), parameter :: diags(2) = ['N', 'U']
    integer :: i, k, u, dg, n
    complex(ep), allocatable :: AP(:)
    real(ep), allocatable :: work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=64) :: label

    call report_init('zlantp', target_name)
    do i = 1, size(ns)
        n = ns(i)
        allocate(AP(n*(n+1)/2))
        call gen_vector_complex(n*(n+1)/2, AP, seed = 19601 + 47 * i)
        allocate(work(max(1, n)))
        do u = 1, size(uplos)
            do dg = 1, size(diags)
                do k = 1, size(norms)
                    ref_val = zlantp(norms(k), uplos(u), diags(dg), n, AP, work)
                    got_val = target_zlantp(norms(k), uplos(u), diags(dg), n, AP)
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
end program test_zlantp
