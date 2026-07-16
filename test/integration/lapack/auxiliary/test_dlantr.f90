! dlantr: norm of a real trapezoidal matrix.
program test_dlantr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dlantr
    use ref_quad_lapack, only: dlantr
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [12, 24]
    character(len=1), parameter :: norms(4) = ['M', '1', 'I', 'F']
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    character(len=1), parameter :: diags(2) = ['N', 'U']
    integer :: i, k, u, dg, m, n
    real(ep), allocatable :: A(:,:), work(:)
    real(ep) :: ref_val, got_val, err, tol
    character(len=80) :: label

    call report_init('dlantr', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A, seed = 19701 + 47 * i)
        allocate(work(max(1, m)))
        do u = 1, size(uplos)
            do dg = 1, size(diags)
                do k = 1, size(norms)
                    ref_val = dlantr(norms(k), uplos(u), diags(dg), m, n, A, m, work)
                    got_val = target_dlantr(norms(k), uplos(u), diags(dg), m, n, A, m)
                    err = rel_err_scalar(got_val, ref_val)
                    tol = 16.0_ep * real(max(m,n), ep) * target_eps
                    write(label, '(a,a,a,a,a,a,a,i0,a,i0)') &
                        'norm=', norms(k), ',uplo=', uplos(u), ',diag=', diags(dg), ',m=', m, ',n=', n
                    call report_case(trim(label), err, tol)
                end do
            end do
        end do
        deallocate(A, work)
    end do
    call report_finalize()
end program test_dlantr
