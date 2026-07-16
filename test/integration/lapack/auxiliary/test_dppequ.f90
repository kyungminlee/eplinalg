program test_dppequ
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec, rel_err_scalar
    use test_data,       only: gen_spd_matrix_quad, pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, target_dppequ
    use ref_quad_lapack, only: dppequ
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, np, info
    real(ep), allocatable :: A(:,:), AP(:), S_ref(:), S_got(:)
    real(ep) :: scnd_ref, am_ref, scnd_got, am_got, err, tol
    character(len=48) :: label

    call report_init('dppequ', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_spd_matrix_quad(n, A, seed = 200031 + 47 * i)
        allocate(AP(np))
        call pack_sym_packed_quad('U', n, A, AP)
        allocate(S_ref(n), S_got(n))
        call dppequ('U', n, AP, S_ref, scnd_ref, am_ref, info)
        call target_dppequ('U', n, AP, S_got, scnd_got, am_got, info)
        err = max(max_rel_err_vec(S_got, S_ref), &
                  rel_err_scalar(scnd_got, scnd_ref), &
                  rel_err_scalar(am_got, am_ref))
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AP, S_ref, S_got)
    end do
    call report_finalize()
end program test_dppequ
