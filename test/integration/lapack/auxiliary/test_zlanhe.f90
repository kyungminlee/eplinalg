program test_zlanhe
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zlanhe
    use ref_quad_lapack, only: zlanhe
    implicit none

    integer, parameter :: ns(*) = [8, 32, 96]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, k, j, n
    complex(ep), allocatable :: A(:,:)
    real(ep),    allocatable :: work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=48) :: label

    call report_init('zlanhe', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A, seed = 13201 + 83 * i)
        allocate(work(max(1, n)))
        do j = 1, size(uplos)
            do k = 1, size(norms)
                ref_val = zlanhe(norms(k), uplos(j), n, A, n, work)
                got_val = target_zlanhe(norms(k), uplos(j), n, A, n)
                err = rel_err_scalar(got_val, ref_val)
                tol = 16.0_ep * real(n * n, ep) * target_eps
                write(label, '(a,a,a,a,a,i0)') 'norm=', norms(k), ',uplo=', &
                    uplos(j), ',n=', n
                call report_case(trim(label), err, tol)
            end do
        end do
        deallocate(work)
    end do
    call report_finalize()
end program test_zlanhe
