! dbdsvdx: bidiagonal SVD via Sturm-sequence selection.
program test_dbdsvdx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use target_lapack,   only: target_name, target_eps, target_dbdsvdx
    use ref_quad_lapack, only: dbdsvdx
    implicit none

    integer, parameter :: ns(*) = [10, 32, 64]
    integer :: i, n, ns_ref, ns_got, info, j
    real(ep), allocatable :: D(:), E(:), S_ref(:), S_got(:), Z(:,:), work(:)
    integer, allocatable  :: iwork(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dbdsvdx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        allocate(D(n), E(max(1, n-1)), S_ref(2*n), S_got(2*n), Z(1, 1), &
                 work(14*n), iwork(12*n))
        do j = 1, n
            D(j) = 2.0_ep + 0.01_ep * real(j, ep)
        end do
        do j = 1, n-1
            E(j) = 0.5_ep
        end do
        S_ref = 0.0_ep; S_got = 0.0_ep
        call dbdsvdx('U', 'N', 'A', n, D, E, 0.0_ep, 0.0_ep, 0, 0, &
                     ns_ref, S_ref, Z, 1, work, iwork, info)
        call target_dbdsvdx('U', 'N', 'A', n, D, E, 0.0_ep, 0.0_ep, 0, 0, &
                            ns_got, S_got, Z, 1, info)
        err = max_rel_err_vec(S_got(1:ns_got), S_ref(1:ns_ref))
        if (ns_ref /= ns_got) err = max(err, 1.0_ep)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(D, E, S_ref, S_got, Z, work, iwork)
    end do
    call report_finalize()
end program test_dbdsvdx
