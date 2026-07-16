! zhsein: complex inverse iteration. SIDE='R', EIGSRC='N', INITV='N'.
program test_zhsein
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zhsein
    use ref_quad_lapack, only: zgees, zhsein
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, sdim, m_r, m_g
    complex(ep), allocatable :: A(:,:), H(:,:), W(:), W_r(:), VS(:,:), work(:)
    complex(ep), allocatable :: VL(:,:), VR_r(:,:), VR_g(:,:)
    real(ep),    allocatable :: rwork(:)
    integer,     allocatable :: ifaill(:), ifailr(:)
    logical,     allocatable :: bwork(:), select(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zhsein', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 113001 + 47 * i)
        allocate(H(n,n), W(n), W_r(n), VS(n,n), rwork(n), bwork(n))
        H = A
        call zgees('N', 'N', sel_all_c, n, H, n, sdim, W, VS, n, &
                   wopt, -1, rwork, bwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgees('N', 'N', sel_all_c, n, H, n, sdim, W, VS, n, &
                   work, lwork, rwork, bwork, info)
        deallocate(work)
        ! Restore H from A; zgees overwrites with Schur form.
        H = A
        W_r = W
        allocate(VL(1, n), VR_r(n, n), VR_g(n, n))
        allocate(ifaill(n), ifailr(n), select(n))
        select = .true.
        VR_r = (0.0_ep, 0.0_ep); VR_g = (0.0_ep, 0.0_ep)
        allocate(work(n * n))
        if (size(rwork) < n) then
            deallocate(rwork); allocate(rwork(n))
        end if
        call zhsein('R', 'N', 'N', select, n, H, n, W_r, &
                    VL, 1, VR_r, n, n, m_r, work, rwork, ifaill, ifailr, info)
        deallocate(work)
        W_r = W
        select = .true.
        call target_zhsein('R', 'N', 'N', select, n, H, n, W_r, &
                           VL, 1, VR_g, n, n, m_g, ifaill, ifailr, info)
        err = max_rel_err_mat_z(VR_g, VR_r)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, H, W, W_r, VS, rwork, bwork, VL, VR_r, VR_g, &
                   ifaill, ifailr, select)
    end do
    call report_finalize()
contains
    logical function sel_all_c(z)
        complex(ep), intent(in) :: z
        sel_all_c = .true.
        if (z == z) continue
    end function
end program test_zhsein
