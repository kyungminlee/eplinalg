! zlanhp: norm of a complex Hermitian packed matrix.
program test_zlanhp
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_hermitian_matrix_quad, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zlanhp
    use ref_quad_lapack, only: zlanhp
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, k, u, n
    complex(ep), allocatable :: A(:,:), AP(:)
    real(ep), allocatable :: work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=48) :: label

    call report_init('zlanhp', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A, seed = 18901 + 47 * i)
        allocate(AP(n*(n+1)/2))
        do u = 1, size(uplos)
            call pack_herm_packed_quad(uplos(u), n, A, AP)
            allocate(work(max(1, n)))
            do k = 1, size(norms)
                ref_val = zlanhp(norms(k), uplos(u), n, AP, work)
                got_val = target_zlanhp(norms(k), uplos(u), n, AP)
                err = rel_err_scalar(got_val, ref_val)
                tol = 16.0_ep * real(n, ep) * target_eps
                write(label, '(a,a,a,a,a,i0)') &
                    'norm=', norms(k), ',uplo=', uplos(u), ',n=', n
                call report_case(trim(label), err, tol)
            end do
            deallocate(work)
        end do
        deallocate(A, AP)
    end do
    call report_finalize()
end program test_zlanhp
