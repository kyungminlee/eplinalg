program test_zpttrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_vector_quad, gen_vector_complex, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zpttrf, target_zpttrs
    use ref_quad_lapack, only: zpttrf, zpttrs
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer, parameter :: nrhs = 2
    integer :: i, n, info, j
    real(ep),    allocatable :: d0(:), d_ref(:), d_got(:)
    complex(ep), allocatable :: e0(:), e_ref(:), e_got(:)
    complex(ep), allocatable :: B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zpttrs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n,    d0, seed = 102001 + 47 * i)
        call gen_vector_complex(n-1, e0, seed = 102011 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 102021 + 47 * i)
        do j = 1, n; d0(j) = abs(d0(j)) + real(4, ep); end do
        allocate(d_ref(n), e_ref(n-1), d_got(n), e_got(n-1))
        allocate(B_ref(n, nrhs), B_got(n, nrhs))
        d_ref = d0; e_ref = e0; d_got = d0; e_got = e0
        call zpttrf(n, d_ref, e_ref, info)
        call target_zpttrf(n, d_got, e_got, info)
        B_ref = B0; B_got = B0
        call zpttrs('U', n, nrhs, d_ref, e_ref, B_ref, n, info)
        call target_zpttrs('U', n, nrhs, d_got, e_got, B_got, n, info)
        err = max_rel_err_mat_z(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(d0, e0, B0, d_ref, e_ref, d_got, e_got, B_ref, B_got)
    end do
    call report_finalize()
end program test_zpttrs
