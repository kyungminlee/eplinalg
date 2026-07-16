! Banded triangular solve. Build a diagonally-dominant upper-triangular
! banded matrix; pack into AB(kd+1, n) with AB(kd+1+i-j, j) = A(i,j).
program test_dtbtrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtbtrs
    use ref_quad_lapack, only: dtbtrs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kd    = 3
    integer, parameter :: nrhs  = 3
    integer :: i, n, ldab, info, j, k
    real(ep), allocatable :: A(:,:), AB(:,:), B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtbtrs', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = kd + 1
        call gen_matrix_quad(n, n, A, seed = 104001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + real(2*n, ep); end do
        allocate(AB(ldab, n)); AB = 0.0_ep
        do j = 1, n
            do k = max(1, j-kd), j
                AB(kd + 1 + k - j, j) = A(k, j)
            end do
        end do
        call gen_matrix_quad(n, nrhs, B0, seed = 104011 + 47 * i)
        allocate(B_ref(n, nrhs), B_got(n, nrhs))
        B_ref = B0; B_got = B0
        call dtbtrs('U', 'N', 'N', n, kd, nrhs, AB, ldab, B_ref, n, info)
        call target_dtbtrs('U', 'N', 'N', n, kd, nrhs, AB, ldab, B_got, n, info)
        err = max_rel_err_mat(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AB, B0, B_ref, B_got)
    end do
    call report_finalize()
end program test_dtbtrs
