! dlantb: norm of a real triangular banded matrix.
program test_dlantb
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dlantb
    use ref_quad_lapack, only: dlantb
    implicit none

    integer, parameter :: ns(*) = [16, 32]
    integer, parameter :: kds(*) = [3, 5]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    character(len=1), parameter :: diags(2) = ['N', 'U']
    integer :: i, k, u, dg, n, kd, ldab
    real(ep), allocatable :: AB(:,:), work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=80) :: label

    call report_init('dlantb', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_matrix_quad(ldab, n, AB, seed = 19301 + 47 * i)
        do u = 1, size(uplos)
            do dg = 1, size(diags)
                allocate(work(max(1, n)))
                do k = 1, size(norms)
                    ref_val = dlantb(norms(k), uplos(u), diags(dg), n, kd, AB, ldab, work)
                    got_val = target_dlantb(norms(k), uplos(u), diags(dg), n, kd, AB, ldab)
                    err = rel_err_scalar(got_val, ref_val)
                    tol = 16.0_ep * real(n, ep) * target_eps
                    write(label, '(a,a,a,a,a,a,a,i0,a,i0)') &
                        'norm=', norms(k), ',uplo=', uplos(u), ',diag=', diags(dg), ',n=', n, ',k=', kd
                    call report_case(trim(label), err, tol)
                end do
                deallocate(work)
            end do
        end do
        deallocate(AB)
    end do
    call report_finalize()
end program test_dlantb
