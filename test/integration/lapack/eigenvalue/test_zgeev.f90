! zgeev: complex non-Hermitian eigenvalue. Compare via sorted moduli
! to dodge per-eigenvector sign/phase ambiguity.
program test_zgeev
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgeev
    use ref_quad_lapack, only: zgeev
    implicit none

    integer, parameter :: ns(*) = [8, 24, 48]
    integer :: i, n, info, lwork, j, k
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), w_ref(:), w_got(:), vl(:,:), vr(:,:), work(:)
    real(ep),    allocatable :: rwork(:), mods_ref(:), mods_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: tmp, err, tol
    character(len=48) :: label

    call report_init('zgeev', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A0, seed = 29001 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), w_ref(n), w_got(n))
        allocate(vl(1,1), vr(1,1), rwork(2*n))
        A_ref = A0; A_got = A0
        call zgeev('N', 'N', n, A_ref, n, w_ref, vl, 1, vr, 1, &
                   wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgeev('N', 'N', n, A_ref, n, w_ref, vl, 1, vr, 1, &
                   work, lwork, rwork, info)
        deallocate(work)
        call target_zgeev('N', 'N', n, A_got, n, w_got, vl, 1, vr, 1, info)
        allocate(mods_ref(n), mods_got(n))
        do j = 1, n
            mods_ref(j) = abs(w_ref(j))
            mods_got(j) = abs(w_got(j))
        end do
        do j = 2, n
            tmp = mods_ref(j); k = j - 1
            do while (k >= 1)
                if (mods_ref(k) >= tmp) exit
                mods_ref(k+1) = mods_ref(k); k = k - 1
            end do
            mods_ref(k+1) = tmp
        end do
        do j = 2, n
            tmp = mods_got(j); k = j - 1
            do while (k >= 1)
                if (mods_got(k) >= tmp) exit
                mods_got(k+1) = mods_got(k); k = k - 1
            end do
            mods_got(k+1) = tmp
        end do
        err = max_rel_err_vec(mods_got, mods_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, w_ref, w_got, vl, vr, rwork, mods_ref, mods_got)
    end do
    call report_finalize()
end program test_zgeev
