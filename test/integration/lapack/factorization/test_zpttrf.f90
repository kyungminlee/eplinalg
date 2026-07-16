! Hermitian PD tridiagonal: real diag d, complex offdiag e.
program test_zpttrf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec, max_rel_err_vec_z
    use test_data,       only: gen_vector_quad, gen_vector_complex
    use target_lapack,   only: target_name, target_eps, target_zpttrf
    use ref_quad_lapack, only: zpttrf
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer :: i, n, info, j
    real(ep),    allocatable :: d0(:), d_ref(:), d_got(:)
    complex(ep), allocatable :: e0(:), e_ref(:), e_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zpttrf', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n,    d0, seed = 101001 + 47 * i)
        call gen_vector_complex(n-1, e0, seed = 101011 + 47 * i)
        do j = 1, n; d0(j) = abs(d0(j)) + real(4, ep); end do
        allocate(d_ref(n), e_ref(n-1), d_got(n), e_got(n-1))
        d_ref = d0; e_ref = e0; d_got = d0; e_got = e0
        call zpttrf(n, d_ref, e_ref, info)
        call target_zpttrf(n, d_got, e_got, info)
        err = max(max_rel_err_vec(d_got, d_ref), max_rel_err_vec_z(e_got, e_ref))
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(d0, e0, d_ref, e_ref, d_got, e_got)
    end do
    call report_finalize()
end program test_zpttrf
