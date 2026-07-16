! ztgsyl: complex coupled generalized Sylvester equations. Smoke test
! with upper-triangular (A,D) and (B,E).
program test_ztgsyl
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztgsyl
    use ref_quad_lapack, only: ztgsyl
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, m, info, j
    complex(ep), allocatable :: A(:,:), B(:,:), D(:,:), E(:,:), C(:,:), F(:,:)
    complex(ep), allocatable :: C_r(:,:), F_r(:,:), C_g(:,:), F_g(:,:)
    real(ep) :: scale_r, dif_r, scale_g, dif_g
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztgsyl', target_name)
    do i = 1, size(ns)
        n = ns(i); m = n
        call gen_matrix_complex(m, m, A, seed = 146101 + 47 * i)
        call gen_matrix_complex(n, n, B, seed = 146111 + 47 * i)
        call gen_matrix_complex(m, m, D, seed = 146121 + 47 * i)
        call gen_matrix_complex(n, n, E, seed = 146131 + 47 * i)
        call gen_matrix_complex(m, n, C, seed = 146141 + 47 * i)
        call gen_matrix_complex(m, n, F, seed = 146151 + 47 * i)
        do j = 1, m
            A(j+1:m, j) = (0.0_ep, 0.0_ep)
            D(j+1:m, j) = (0.0_ep, 0.0_ep)
        end do
        do j = 1, n
            B(j+1:n, j) = (0.0_ep, 0.0_ep)
            E(j+1:n, j) = (0.0_ep, 0.0_ep)
        end do
        do j = 1, m
            A(j,j) = A(j,j) + (2.0_ep, 0.0_ep)
            D(j,j) = D(j,j) + (2.0_ep, 0.0_ep)
        end do
        do j = 1, n
            B(j,j) = B(j,j) + (2.0_ep, 0.0_ep)
            E(j,j) = E(j,j) + (2.0_ep, 0.0_ep)
        end do
        allocate(C_r(m,n), F_r(m,n), C_g(m,n), F_g(m,n))
        C_r = C; F_r = F; C_g = C; F_g = F
        call ztgsyl_call(C_r, F_r, scale_r, dif_r)
        call target_ztgsyl('N', 0, m, n, A, m, B, n, C_g, m, &
                           D, m, E, n, F_g, m, scale_g, dif_g, info)
        err = max(max_rel_err_mat_z(C_g, C_r), max_rel_err_mat_z(F_g, F_r))
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, D, E, C, F, C_r, F_r, C_g, F_g)
    end do
    call report_finalize()
contains
    subroutine ztgsyl_call(Cw, Fw, sc, df)
        complex(ep), intent(inout) :: Cw(m,n), Fw(m,n)
        real(ep), intent(out) :: sc, df
        complex(ep) :: wopt(1)
        complex(ep), allocatable :: work(:)
        integer, allocatable :: iwork(:)
        integer :: lwork, infl
        allocate(iwork(m+n+2))
        call ztgsyl('N', 0, m, n, A, m, B, n, Cw, m, D, m, E, n, Fw, m, &
                    sc, df, wopt, -1, iwork, infl)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call ztgsyl('N', 0, m, n, A, m, B, n, Cw, m, D, m, E, n, Fw, m, &
                    sc, df, work, lwork, iwork, infl)
        deallocate(work, iwork)
    end subroutine
end program test_ztgsyl
