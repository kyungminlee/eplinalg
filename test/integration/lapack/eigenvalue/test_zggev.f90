program test_zggev
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zggev
    use ref_quad_lapack, only: zggev
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, j
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:)
    complex(ep), allocatable :: al_ref(:), be_ref(:), al_got(:), be_got(:)
    complex(ep), allocatable :: VL(:,:), VR(:,:), work(:)
    real(ep), allocatable :: rwork(:), ev_ref(:), ev_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zggev', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A0, seed = 360021 + 47 * i)
        call gen_matrix_complex(n, n, B0, seed = 360031 + 47 * i)
        allocate(A_ref(n,n), B_ref(n,n), A_got(n,n), B_got(n,n))
        allocate(al_ref(n), be_ref(n), al_got(n), be_got(n))
        allocate(VL(1,1), VR(1,1), rwork(8*n))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        call zggev('N', 'N', n, A_ref, n, B_ref, n, al_ref, be_ref, &
                   VL, 1, VR, 1, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zggev('N', 'N', n, A_ref, n, B_ref, n, al_ref, be_ref, &
                   VL, 1, VR, 1, work, lwork, rwork, info)
        call target_zggev('N', 'N', n, A_got, n, B_got, n, al_got, be_got, &
                          VL, 1, VR, 1, info)
        allocate(ev_ref(n), ev_got(n))
        do j = 1, n
            ev_ref(j) = abs(al_ref(j)) / max(abs(be_ref(j)), tiny(1.0_ep))
            ev_got(j) = abs(al_got(j)) / max(abs(be_got(j)), tiny(1.0_ep))
        end do
        call sort_asc(ev_ref, n); call sort_asc(ev_got, n)
        err = max_rel_err_vec(ev_got, ev_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, A_got, B_got, al_ref, be_ref, al_got, be_got, &
                   VL, VR, work, rwork, ev_ref, ev_got)
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
end program test_zggev
