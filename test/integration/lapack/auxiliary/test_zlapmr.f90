! zlapmr: row permutations P*X (complex).
program test_zlapmr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zlapmr
    use ref_quad_lapack, only: zlapmr
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [12, 24]
    logical, parameter :: forwds(2) = [.true., .false.]
    integer :: i, k, m, n, ii
    complex(ep), allocatable :: X0(:,:), X_r(:,:), X_g(:,:)
    integer, allocatable :: K_r(:), K_g(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zlapmr', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, X0, seed = 19101 + 47 * i)
        do k = 1, size(forwds)
            allocate(X_r(m,n), X_g(m,n), K_r(m), K_g(m))
            X_r = X0; X_g = X0
            do ii = 1, m; K_r(ii) = mod(ii * 7 + 3, m) + 1; end do
            K_g = K_r
            call zlapmr(forwds(k), m, n, X_r, m, K_r)
            call target_zlapmr(forwds(k), m, n, X_g, m, K_g)
            err = max_rel_err_mat_z(X_g, X_r)
            tol = 16.0_ep * target_eps
            write(label, '(a,l1,a,i0,a,i0)') 'fw=', forwds(k), ',m=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(X_r, X_g, K_r, K_g)
        end do
        deallocate(X0)
    end do
    call report_finalize()
end program test_zlapmr
