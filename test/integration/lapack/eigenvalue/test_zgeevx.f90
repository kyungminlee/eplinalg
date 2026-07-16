! zgeevx: complex non-Hermitian eig. BALANC='N', JOBV*='N', SENSE='N'.
program test_zgeevx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgeevx
    use ref_quad_lapack, only: zgeevx
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, ilo_r, ihi_r, ilo_g, ihi_g, j
    complex(ep), allocatable :: A(:,:), Aref(:,:), Agot(:,:)
    complex(ep), allocatable :: W_r(:), W_g(:), VL(:,:), VR(:,:), work(:)
    real(ep),    allocatable :: scale_r(:), scale_g(:), rce_r(:), rcv_r(:)
    real(ep),    allocatable :: rce_g(:), rcv_g(:), rwork(:), mods_r(:), mods_g(:)
    complex(ep) :: wopt(1)
    real(ep) :: abnrm_r, abnrm_g, err, tol
    character(len=48) :: label

    call report_init('zgeevx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 105001 + 47 * i)
        allocate(Aref(n,n), Agot(n,n), W_r(n), W_g(n))
        allocate(VL(1,1), VR(1,1), scale_r(n), scale_g(n))
        allocate(rce_r(n), rcv_r(n), rce_g(n), rcv_g(n))
        allocate(rwork(2*n), mods_r(n), mods_g(n))
        Aref = A; Agot = A
        call zgeevx('N', 'N', 'N', 'N', n, Aref, n, W_r, &
                    VL, 1, VR, 1, ilo_r, ihi_r, scale_r, abnrm_r, &
                    rce_r, rcv_r, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgeevx('N', 'N', 'N', 'N', n, Aref, n, W_r, &
                    VL, 1, VR, 1, ilo_r, ihi_r, scale_r, abnrm_r, &
                    rce_r, rcv_r, work, lwork, rwork, info)
        deallocate(work)
        call target_zgeevx('N', 'N', 'N', 'N', n, Agot, n, W_g, &
                           VL, 1, VR, 1, ilo_g, ihi_g, scale_g, abnrm_g, &
                           rce_g, rcv_g, info)
        do j = 1, n
            mods_r(j) = abs(W_r(j))
            mods_g(j) = abs(W_g(j))
        end do
        call sort_desc(mods_r, n); call sort_desc(mods_g, n)
        err = max_rel_err_vec(mods_g, mods_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(Aref, Agot, W_r, W_g, VL, VR, &
                   scale_r, scale_g, rce_r, rcv_r, rce_g, rcv_g, &
                   rwork, mods_r, mods_g)
    end do
    call report_finalize()
contains
    subroutine sort_desc(x, m)
        real(ep), intent(inout) :: x(:)
        integer,  intent(in)    :: m
        integer :: ii, jj
        real(ep) :: tt
        do ii = 1, m - 1
            do jj = ii + 1, m
                if (x(ii) < x(jj)) then
                    tt = x(ii); x(ii) = x(jj); x(jj) = tt
                end if
            end do
        end do
    end subroutine
end program test_zgeevx
