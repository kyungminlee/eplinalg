! zlangb: norm of a banded complex matrix.
program test_zlangb
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zlangb
    use ref_quad_lapack, only: zlangb
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kls(*) = [2, 4, 8]
    integer, parameter :: kus(*) = [3, 5, 7]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    integer :: i, k, n, kl, ku, ldab
    complex(ep), allocatable :: AB(:,:)
    real(ep), allocatable :: work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=64) :: label

    call report_init('zlangb', target_name)
    do i = 1, size(ns)
        n = ns(i); kl = kls(i); ku = kus(i); ldab = kl + ku + 1
        call gen_matrix_complex(ldab, n, AB, seed = 18101 + 47 * i)
        allocate(work(max(1, n)))
        do k = 1, size(norms)
            ref_val = zlangb(norms(k), n, kl, ku, AB, ldab, work)
            got_val = target_zlangb(norms(k), n, kl, ku, AB, ldab)
            err = rel_err_scalar(got_val, ref_val)
            tol = 16.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,i0,a,i0,a,i0)') 'norm=', norms(k), ',n=', n, ',kl=', kl, ',ku=', ku
            call report_case(trim(label), err, tol)
        end do
        deallocate(AB, work)
    end do
    call report_finalize()
end program test_zlangb
