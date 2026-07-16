! dpotrs test: run dpotrf+dpotrs on each side independently and
! compare solution B.
program test_dpotrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_spd_matrix_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, &
                                target_dpotrf, target_dpotrs
    use ref_quad_lapack, only: dpotrf, dpotrs
    implicit none

    integer, parameter :: ns(*) = [8, 32, 64]
    integer, parameter :: nrhs  = 2
    integer :: i, n, info
    real(ep), allocatable :: A0(:,:), B0(:,:)
    real(ep), allocatable :: A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dpotrs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_spd_matrix_quad(n,       A0, seed = 5001 + 47 * i)
        call gen_matrix_quad(n, nrhs, B0, seed = 5011 + 47 * i)

        allocate(A_ref(n, n), B_ref(n, nrhs), A_got(n, n), B_got(n, nrhs))
        A_ref = A0;  B_ref = B0;  A_got = A0;  B_got = B0

        call dpotrf('U', n, A_ref, n, info)
        call dpotrs('U', n, nrhs, A_ref, n, B_ref, n, info)
        call target_dpotrf('U', n, A_got, n, info)
        call target_dpotrs('U', n, nrhs, A_got, n, B_got, n, info)

        err = max_rel_err_mat(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs
        call report_case(trim(label), err, tol)

        deallocate(A_ref, B_ref, A_got, B_got)
    end do
    call report_finalize()
end program test_dpotrs
