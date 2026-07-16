program test_dtbcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtbcon
    use ref_quad_lapack, only: dtbcon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kd = 3
    integer :: i, n, ldab, info, j, k
    integer, allocatable :: iwork(:)
    real(ep), allocatable :: A(:,:), AB(:,:), work(:)
    real(ep) :: rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('dtbcon', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = kd + 1
        call gen_matrix_quad(n, n, A, seed = 132001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + real(2*n, ep); end do
        allocate(AB(ldab, n), work(3*n), iwork(n))
        AB = 0.0_ep
        do j = 1, n
            do k = max(1, j-kd), j
                AB(kd + 1 + k - j, j) = A(k, j)
            end do
        end do
        call dtbcon('1', 'U', 'N', n, kd, AB, ldab, rcond_ref, work, iwork, info)
        call target_dtbcon('1', 'U', 'N', n, kd, AB, ldab, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AB, work, iwork)
    end do
    call report_finalize()
end program test_dtbcon
