! ztrevc3: blocked TREVC for complex Schur eigenvectors. SIDE='R',
! HOWMNY='A'.
program test_ztrevc3
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztrevc3
    use ref_quad_lapack, only: zgees, ztrevc3
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, lrwork, sdim, m_r, m_g
    complex(ep), allocatable :: A(:,:), T(:,:), W(:), VS(:,:), work(:)
    complex(ep), allocatable :: VL(:,:), VR_r(:,:), VR_g(:,:)
    real(ep),    allocatable :: rwork(:)
    logical,     allocatable :: bwork(:), select(:)
    complex(ep) :: wopt(1)
    real(ep) :: rwopt(1), err, tol
    character(len=48) :: label

    call report_init('ztrevc3', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 111001 + 47 * i)
        allocate(T(n,n), W(n), VS(n,n), rwork(n), bwork(n))
        T = A
        call zgees('N', 'N', sel_all_c, n, T, n, sdim, W, VS, n, &
                   wopt, -1, rwork, bwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgees('N', 'N', sel_all_c, n, T, n, sdim, W, VS, n, &
                   work, lwork, rwork, bwork, info)
        deallocate(work)
        allocate(VL(1, n), VR_r(n, n), VR_g(n, n), select(n))
        select = .false.
        call ztrevc3('R', 'A', select, n, T, n, VL, 1, VR_r, n, n, m_r, &
                     wopt, -1, rwopt, -1, info)
        lwork  = max(1, int(real(wopt(1), ep)))
        lrwork = max(1, int(rwopt(1)))
        allocate(work(lwork))
        if (size(rwork) < lrwork) then
            deallocate(rwork); allocate(rwork(lrwork))
        end if
        call ztrevc3('R', 'A', select, n, T, n, VL, 1, VR_r, n, n, m_r, &
                     work, lwork, rwork, lrwork, info)
        deallocate(work)
        call target_ztrevc3('R', 'A', select, n, T, n, VL, 1, VR_g, n, n, &
                            m_g, info)
        err = max_rel_err_mat_z(VR_g, VR_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, T, W, VS, rwork, bwork, VL, VR_r, VR_g, select)
    end do
    call report_finalize()
contains
    logical function sel_all_c(z)
        complex(ep), intent(in) :: z
        sel_all_c = .true.
        if (z == z) continue
    end function
end program test_ztrevc3
