! zlansb: norm of a complex symmetric banded matrix.
program test_zlansb
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_complex_symmetric_quad, pack_herm_band_quad
    use target_lapack,   only: target_name, target_eps, target_zlansb
    use ref_quad_lapack, only: zlansb
    implicit none

    integer, parameter :: ns(*) = [16, 32]
    integer, parameter :: kds(*) = [3, 5]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, k, u, n, kd, ldab
    complex(ep), allocatable :: A(:,:), AB(:,:)
    real(ep), allocatable :: work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=64) :: label

    call report_init('zlansb', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_complex_symmetric_quad(n, A, seed = 18751 + 47 * i)
        allocate(AB(ldab, n))
        do u = 1, size(uplos)
            call pack_herm_band_quad(uplos(u), n, kd, A, AB)
            allocate(work(max(1, n)))
            do k = 1, size(norms)
                ref_val = zlansb(norms(k), uplos(u), n, kd, AB, ldab, work)
                got_val = target_zlansb(norms(k), uplos(u), n, kd, AB, ldab)
                err = rel_err_scalar(got_val, ref_val)
                tol = 16.0_ep * real(n, ep) * target_eps
                write(label, '(a,a,a,a,a,i0,a,i0)') &
                    'norm=', norms(k), ',uplo=', uplos(u), ',n=', n, ',k=', kd
                call report_case(trim(label), err, tol)
            end do
            deallocate(work)
        end do
        deallocate(A, AB)
    end do
    call report_finalize()
end program test_zlansb
