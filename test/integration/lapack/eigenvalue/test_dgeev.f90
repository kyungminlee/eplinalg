! dgeev: general non-symmetric eigenvalue. Real-arithmetic eigenvalues
! come back as (wr, wi) pairs. Compare via the sign-and-permutation-
! immune metric: sort moduli of (wr + i*wi) and compare descending.
program test_dgeev
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgeev
    use ref_quad_lapack, only: dgeev
    implicit none

    integer, parameter :: ns(*) = [8, 24, 48]
    integer :: i, n, info, lwork, j, k
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: wr_ref(:), wi_ref(:), wr_got(:), wi_got(:)
    real(ep), allocatable :: vl(:,:), vr(:,:), work(:)
    real(ep), allocatable :: mods_ref(:), mods_got(:)
    real(ep) :: wopt(1), tmp, err, tol
    character(len=48) :: label

    call report_init('dgeev', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A0, seed = 20001 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n))
        allocate(wr_ref(n), wi_ref(n), wr_got(n), wi_got(n))
        allocate(vl(1,1), vr(1,1))
        A_ref = A0; A_got = A0
        call dgeev('N', 'N', n, A_ref, n, wr_ref, wi_ref, vl, 1, vr, 1, &
                   wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgeev('N', 'N', n, A_ref, n, wr_ref, wi_ref, vl, 1, vr, 1, &
                   work, lwork, info)
        deallocate(work)
        call target_dgeev('N', 'N', n, A_got, n, wr_got, wi_got, vl, 1, vr, 1, info)
        allocate(mods_ref(n), mods_got(n))
        do j = 1, n
            mods_ref(j) = sqrt(wr_ref(j)**2 + wi_ref(j)**2)
            mods_got(j) = sqrt(wr_got(j)**2 + wi_got(j)**2)
        end do
        ! Insertion sort, descending.
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
        deallocate(A_ref, A_got, wr_ref, wi_ref, wr_got, wi_got, vl, vr, &
                   mods_ref, mods_got)
    end do
    call report_finalize()
end program test_dgeev
