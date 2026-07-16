! dlansf: norm of a real symmetric matrix in RFP storage.
program test_dlansf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_symmetric_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dlansf
    use ref_quad_lapack, only: dlansf, dtrttf
    implicit none

    integer, parameter :: ns(*) = [7, 16, 33]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    character(len=1), parameter :: transrs(2) = ['N', 'T']
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, k, t, u, n, info, nt
    real(ep), allocatable :: A(:,:), ARF(:), work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=64) :: label

    call report_init('dlansf', target_name)
    do i = 1, size(ns)
        n = ns(i)
        nt = n * (n + 1) / 2
        call gen_symmetric_matrix_quad(n, A, seed = 17501 + 47 * i)
        allocate(ARF(nt), work(max(1, n)))
        do t = 1, size(transrs)
            do u = 1, size(uplos)
                call dtrttf(transrs(t), uplos(u), n, A, n, ARF, info)
                do k = 1, size(norms)
                    ref_val = dlansf(norms(k), transrs(t), uplos(u), n, ARF, work)
                    got_val = target_dlansf(norms(k), transrs(t), uplos(u), n, ARF)
                    err = rel_err_scalar(got_val, ref_val)
                    tol = 16.0_ep * real(n, ep) * target_eps
                    write(label, '(a,a,a,a,a,a,a,i0)') &
                        'norm=', norms(k), ',transr=', transrs(t), &
                        ',uplo=', uplos(u), ',n=', n
                    call report_case(trim(label), err, tol)
                end do
            end do
        end do
        deallocate(A, ARF, work)
    end do
    call report_finalize()
end program test_dlansf
