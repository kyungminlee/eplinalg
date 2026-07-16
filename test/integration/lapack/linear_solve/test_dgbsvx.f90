program test_dgbsvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, rel_err_scalar
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgbsvx
    use ref_quad_lapack, only: dgbsvx
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer, parameter :: kl    = 2
    integer, parameter :: ku    = 3
    integer, parameter :: nrhs  = 2
    integer :: i, n, ldab, ldafb, info, j, k
    real(ep), allocatable :: A(:,:), B0(:,:)
    real(ep), allocatable :: AB_ref(:,:), AFB_ref(:,:), B_ref(:,:), X_ref(:,:), R_ref(:), C_ref(:)
    real(ep), allocatable :: AB_got(:,:), AFB_got(:,:), B_got(:,:), X_got(:,:), R_got(:), C_got(:)
    real(ep), allocatable :: ferr_ref(:), berr_ref(:), ferr_got(:), berr_got(:), work(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:), iwork(:)
    real(ep) :: rcond_ref, rcond_got, err, tol
    character :: equed_ref, equed_got
    character(len=48) :: label

    call report_init('dgbsvx', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = kl + ku + 1; ldafb = 2*kl + ku + 1
        call gen_matrix_quad(n, n, A, seed = 420001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + real(2*n, ep); end do
        call gen_matrix_quad(n, nrhs, B0, seed = 420011 + 47 * i)
        allocate(AB_ref(ldab, n), AFB_ref(ldafb, n), B_ref(n, nrhs), X_ref(n, nrhs), R_ref(n), C_ref(n))
        allocate(AB_got(ldab, n), AFB_got(ldafb, n), B_got(n, nrhs), X_got(n, nrhs), R_got(n), C_got(n))
        allocate(ferr_ref(nrhs), berr_ref(nrhs), ferr_got(nrhs), berr_got(nrhs))
        allocate(ipiv_ref(n), ipiv_got(n), work(3*n), iwork(n))
        AB_ref = 0.0_ep
        do j = 1, n
            do k = max(1, j - ku), min(n, j + kl)
                AB_ref(ku + 1 + k - j, j) = A(k, j)
            end do
        end do
        AB_got = AB_ref
        B_ref = B0; B_got = B0
        equed_ref = 'N'; equed_got = 'N'
        call dgbsvx('E', 'N', n, kl, ku, nrhs, AB_ref, ldab, AFB_ref, ldafb, ipiv_ref, &
                    equed_ref, R_ref, C_ref, B_ref, n, X_ref, n, rcond_ref, &
                    ferr_ref, berr_ref, work, iwork, info)
        call target_dgbsvx('E', 'N', n, kl, ku, nrhs, AB_got, ldab, AFB_got, ldafb, ipiv_got, &
                           equed_got, R_got, C_got, B_got, n, X_got, n, rcond_got, &
                           ferr_got, berr_got, info)
        err = max(max_rel_err_mat(X_got, X_ref), rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B0, AB_ref, AFB_ref, B_ref, X_ref, R_ref, C_ref, &
                   AB_got, AFB_got, B_got, X_got, R_got, C_got, &
                   ferr_ref, berr_ref, ferr_got, berr_got, ipiv_ref, ipiv_got, work, iwork)
    end do
    call report_finalize()
end program test_dgbsvx
