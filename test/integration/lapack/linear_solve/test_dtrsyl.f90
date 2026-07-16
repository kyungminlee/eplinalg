! Sylvester equation AX + XB = sgn*C with A, B upper triangular.
program test_dtrsyl
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, rel_err_scalar
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtrsyl
    use ref_quad_lapack, only: dtrsyl
    implicit none

    integer, parameter :: ms(*) = [12, 24, 32]
    integer, parameter :: ns(*) = [16, 20, 28]
    integer :: i, m, n, info, j, k
    real(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: scale_ref, scale_got, err, tol
    character(len=48) :: label

    call report_init('dtrsyl', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, m, A, seed = 330001 + 47 * i)
        call gen_matrix_quad(n, n, B, seed = 330011 + 47 * i)
        do j = 1, m
            do k = j + 1, m; A(k, j) = 0.0_ep; end do
            A(j, j) = A(j, j) + real(2*m, ep)
        end do
        do j = 1, n
            do k = j + 1, n; B(k, j) = 0.0_ep; end do
            B(j, j) = B(j, j) - real(2*n, ep)
        end do
        call gen_matrix_quad(m, n, C0, seed = 330021 + 47 * i)
        allocate(C_ref(m, n), C_got(m, n))
        C_ref = C0; C_got = C0
        call dtrsyl('N', 'N', 1, m, n, A, m, B, n, C_ref, m, scale_ref, info)
        call target_dtrsyl('N', 'N', 1, m, n, A, m, B, n, C_got, m, scale_got, info)
        err = max(max_rel_err_mat(C_got, C_ref), rel_err_scalar(scale_got, scale_ref))
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, C0, C_ref, C_got)
    end do
    call report_finalize()
end program test_dtrsyl
