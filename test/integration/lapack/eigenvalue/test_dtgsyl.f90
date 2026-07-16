! dtgsyl: solve coupled generalized Sylvester equations. Smoke test
! with upper-triangular (A,D) and (B,E) (valid generalized Schur form).
program test_dtgsyl
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtgsyl
    use ref_quad_lapack, only: dtgsyl
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, m, info, j
    real(ep), allocatable :: A(:,:), B(:,:), D(:,:), E(:,:), C(:,:), F(:,:)
    real(ep), allocatable :: C_r(:,:), F_r(:,:), C_g(:,:), F_g(:,:)
    real(ep) :: scale_r, dif_r, scale_g, dif_g
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtgsyl', target_name)
    do i = 1, size(ns)
        n = ns(i); m = n
        call gen_matrix_quad(m, m, A, seed = 146001 + 47 * i)
        call gen_matrix_quad(n, n, B, seed = 146011 + 47 * i)
        call gen_matrix_quad(m, m, D, seed = 146021 + 47 * i)
        call gen_matrix_quad(n, n, E, seed = 146031 + 47 * i)
        call gen_matrix_quad(m, n, C, seed = 146041 + 47 * i)
        call gen_matrix_quad(m, n, F, seed = 146051 + 47 * i)
        do j = 1, m; A(j+1:m, j) = 0.0_ep; D(j+1:m, j) = 0.0_ep; end do
        do j = 1, n; B(j+1:n, j) = 0.0_ep; E(j+1:n, j) = 0.0_ep; end do
        ! Strengthen diagonal to ensure (A,D)/(B,E) regularity.
        do j = 1, m
            A(j,j) = A(j,j) + 2.0_ep
            D(j,j) = D(j,j) + 2.0_ep
        end do
        do j = 1, n
            B(j,j) = B(j,j) + 2.0_ep
            E(j,j) = E(j,j) + 2.0_ep
        end do
        allocate(C_r(m,n), F_r(m,n), C_g(m,n), F_g(m,n))
        C_r = C; F_r = F; C_g = C; F_g = F
        call dtgsyl_call(C_r, F_r, scale_r, dif_r)
        call target_dtgsyl('N', 0, m, n, A, m, B, n, C_g, m, &
                           D, m, E, n, F_g, m, scale_g, dif_g, info)
        err = max(max_rel_err_mat(C_g, C_r), max_rel_err_mat(F_g, F_r))
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, D, E, C, F, C_r, F_r, C_g, F_g)
    end do
    call report_finalize()
contains
    subroutine dtgsyl_call(Cw, Fw, sc, df)
        real(ep), intent(inout) :: Cw(m,n), Fw(m,n)
        real(ep), intent(out) :: sc, df
        real(ep) :: wopt(1)
        real(ep), allocatable :: work(:)
        integer, allocatable :: iwork(:)
        integer :: lwork, infl
        allocate(iwork(m+n+6))
        call dtgsyl('N', 0, m, n, A, m, B, n, Cw, m, D, m, E, n, Fw, m, &
                    sc, df, wopt, -1, iwork, infl)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dtgsyl('N', 0, m, n, A, m, B, n, Cw, m, D, m, E, n, Fw, m, &
                    sc, df, work, lwork, iwork, infl)
        deallocate(work, iwork)
    end subroutine
end program test_dtgsyl
