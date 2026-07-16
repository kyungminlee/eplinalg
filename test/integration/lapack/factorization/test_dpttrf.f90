! Positive-definite tridiagonal: 2 vectors (diag d, subdiag e).
! Make d strongly positive so it's actually SPD.
program test_dpttrf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_dpttrf
    use ref_quad_lapack, only: dpttrf
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer :: i, n, info, j
    real(ep), allocatable :: d0(:), e0(:), d_ref(:), e_ref(:), d_got(:), e_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dpttrf', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n,   d0, seed = 68001 + 47 * i)
        call gen_vector_quad(n-1, e0, seed = 68011 + 47 * i)
        ! Ensure SPD: |d| >= 2*|e| for strict diagonal dominance.
        do j = 1, n
            d0(j) = abs(d0(j)) + real(4, ep)
        end do
        allocate(d_ref(n), e_ref(n-1), d_got(n), e_got(n-1))
        d_ref = d0; e_ref = e0; d_got = d0; e_got = e0
        call dpttrf(n, d_ref, e_ref, info)
        call target_dpttrf(n, d_got, e_got, info)
        err = max(max_rel_err_vec(d_got, d_ref), max_rel_err_vec(e_got, e_ref))
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(d0, e0, d_ref, e_ref, d_got, e_got)
    end do
    call report_finalize()
end program test_dpttrf
