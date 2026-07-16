program test_ztrevc
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztrevc
    use ref_quad_lapack, only: zgehrd, zhseqr, ztrevc
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, m_ref, m_got
    complex(ep), allocatable :: A(:,:), tau(:), work(:)
    complex(ep), allocatable :: T_ref(:,:), T_got(:,:), W(:), Z(:,:)
    complex(ep), allocatable :: VL(:,:), VR_ref(:,:), VR_got(:,:)
    logical,  allocatable :: sel(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztrevc', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 350011 + 47 * i)
        allocate(tau(n-1))
        call zgehrd(n, 1, n, A, n, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgehrd(n, 1, n, A, n, tau, work, lwork, info)
        deallocate(work)
        allocate(T_ref(n,n), W(n), Z(1,1))
        T_ref = A
        call zhseqr('E', 'N', n, 1, n, T_ref, n, W, Z, 1, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zhseqr('E', 'N', n, 1, n, T_ref, n, W, Z, 1, work, lwork, info)
        deallocate(work)
        allocate(T_got(n,n)); T_got = T_ref
        allocate(VL(1, n), VR_ref(n, n), VR_got(n, n))
        allocate(sel(n)); sel = .false.
        allocate(work(2*n))
        block
            real(ep), allocatable :: rwork(:)
            allocate(rwork(n))
            call ztrevc('R', 'A', sel, n, T_ref, n, VL, 1, VR_ref, n, n, m_ref, &
                        work, rwork, info)
            deallocate(rwork)
        end block
        call target_ztrevc('R', 'A', sel, n, T_got, n, VL, 1, VR_got, n, n, m_got, info)
        err = max_rel_err_mat_z(VR_got, VR_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, tau, T_ref, T_got, W, Z, VL, VR_ref, VR_got, sel, work)
    end do
    call report_finalize()
end program test_ztrevc
