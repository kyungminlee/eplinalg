! Packed triangular solve. Build a diagonally-dominant upper-triangular
! matrix in dense form, then convert to packed.
program test_dtptrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad, pack_sym_packed_quad
    use target_lapack,   only: target_name, target_eps, target_dtptrs
    use ref_quad_lapack, only: dtptrs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs  = 3
    integer :: i, n, np, info, j
    real(ep), allocatable :: A(:,:), AP(:), B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtptrs', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_matrix_quad(n, n, A, seed = 106001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + real(2*n, ep); end do
        allocate(AP(np))
        call pack_sym_packed_quad('U', n, A, AP)
        call gen_matrix_quad(n, nrhs, B0, seed = 106011 + 47 * i)
        allocate(B_ref(n, nrhs), B_got(n, nrhs))
        B_ref = B0; B_got = B0
        call dtptrs('U', 'N', 'N', n, nrhs, AP, B_ref, n, info)
        call target_dtptrs('U', 'N', 'N', n, nrhs, AP, B_got, n, info)
        err = max_rel_err_mat(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AP, B0, B_ref, B_got)
    end do
    call report_finalize()
end program test_dtptrs
