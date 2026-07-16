program test_dsyequb
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec, rel_err_scalar
    use test_data,       only: gen_spd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsyequb
    use ref_quad_lapack, only: dsyequb
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info
    real(ep), allocatable :: A(:,:), S_ref(:), S_got(:), work(:)
    real(ep) :: scnd_ref, am_ref, scnd_got, am_got, err, tol
    character(len=48) :: label

    call report_init('dsyequb', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_spd_matrix_quad(n, A, seed = 210031 + 47 * i)
        allocate(S_ref(n), S_got(n), work(2*n))
        call dsyequb('U', n, A, n, S_ref, scnd_ref, am_ref, work, info)
        call target_dsyequb('U', n, A, n, S_got, scnd_got, am_got, info)
        err = max(max_rel_err_vec(S_got, S_ref), &
                  rel_err_scalar(scnd_got, scnd_ref), &
                  rel_err_scalar(am_got, am_ref))
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, S_ref, S_got, work)
    end do
    call report_finalize()
end program test_dsyequb
