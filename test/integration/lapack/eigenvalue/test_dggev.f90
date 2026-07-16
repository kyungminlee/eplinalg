! dggev: generalized non-symmetric eigenvalue problem A*x = lambda*B*x.
! JOBVL='N', JOBVR='N' — eigenvalues only via (alphar, alphai, beta).
program test_dggev
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dggev
    use ref_quad_lapack, only: dggev
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork
    real(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:)
    real(ep), allocatable :: ar_ref(:), ai_ref(:), be_ref(:)
    real(ep), allocatable :: ar_got(:), ai_got(:), be_got(:)
    real(ep), allocatable :: VL(:,:), VR(:,:), work(:)
    real(ep), allocatable :: ev_ref(:), ev_got(:)
    real(ep) :: wopt(1), err, tol
    integer  :: j
    character(len=48) :: label

    call report_init('dggev', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A0, seed = 360001 + 47 * i)
        call gen_matrix_quad(n, n, B0, seed = 360011 + 47 * i)
        allocate(A_ref(n,n), B_ref(n,n), A_got(n,n), B_got(n,n))
        allocate(ar_ref(n), ai_ref(n), be_ref(n), ar_got(n), ai_got(n), be_got(n))
        allocate(VL(1,1), VR(1,1))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        call dggev('N', 'N', n, A_ref, n, B_ref, n, ar_ref, ai_ref, be_ref, &
                   VL, 1, VR, 1, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dggev('N', 'N', n, A_ref, n, B_ref, n, ar_ref, ai_ref, be_ref, &
                   VL, 1, VR, 1, work, lwork, info)
        call target_dggev('N', 'N', n, A_got, n, B_got, n, ar_got, ai_got, be_got, &
                          VL, 1, VR, 1, info)
        ! Compare sorted |alpha/beta| (eigenvalue moduli) for order-independence.
        allocate(ev_ref(n), ev_got(n))
        do j = 1, n
            ev_ref(j) = sqrt(ar_ref(j)**2 + ai_ref(j)**2) / max(abs(be_ref(j)), tiny(1.0_ep))
            ev_got(j) = sqrt(ar_got(j)**2 + ai_got(j)**2) / max(abs(be_got(j)), tiny(1.0_ep))
        end do
        call sort_asc(ev_ref, n); call sort_asc(ev_got, n)
        err = max_rel_err_vec(ev_got, ev_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, A_got, B_got, &
                   ar_ref, ai_ref, be_ref, ar_got, ai_got, be_got, &
                   VL, VR, work, ev_ref, ev_got)
    end do
    call report_finalize()
contains
    subroutine sort_asc(x, n)
        real(ep), intent(inout) :: x(:)
        integer,  intent(in)    :: n
        integer :: ii, jj
        real(ep) :: tt
        do ii = 1, n - 1
            do jj = ii + 1, n
                if (x(ii) > x(jj)) then
                    tt = x(ii); x(ii) = x(jj); x(jj) = tt
                end if
            end do
        end do
    end subroutine
end program test_dggev
