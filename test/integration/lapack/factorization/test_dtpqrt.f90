! dtpqrt: pentagonal QR factorization of [A; B].
program test_dtpqrt
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtpqrt
    use ref_quad_lapack, only: dtpqrt
    implicit none

    integer, parameter :: ms(*) = [12, 24]
    integer, parameter :: ns(*) = [6,  12]
    integer, parameter :: ls(*) = [0,  3]
    integer, parameter :: nbs(*) = [3, 4]
    integer :: i, m, n, l, nb, j, info
    real(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: B_ref(:,:), B_got(:,:), T_ref(:,:), T_got(:,:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('dtpqrt', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); l = ls(i); nb = nbs(i)
        call gen_matrix_quad(n, n, A0, seed = 22001 + 71 * i)
        ! upper triangular A
        do j = 1, n-1
            A0(j+1:n, j) = 0.0_ep
        end do
        call gen_matrix_quad(m, n, B0, seed = 22011 + 71 * i)
        allocate(A_ref(n, n), A_got(n, n), B_ref(m, n), B_got(m, n), &
                 T_ref(nb, n), T_got(nb, n))
        A_ref = A0; A_got = A0; B_ref = B0; B_got = B0
        T_ref = 0.0_ep; T_got = 0.0_ep
        block
            real(ep), allocatable :: work(:)
            allocate(work(nb*n))
            call dtpqrt(m, n, l, nb, A_ref, n, B_ref, m, T_ref, nb, work, info)
            deallocate(work)
        end block
        call target_dtpqrt(m, n, l, nb, A_got, n, B_got, m, T_got, nb, info)
        err = max(max_rel_err_mat(A_got, A_ref), max_rel_err_mat(B_got, B_ref), &
                  max_rel_err_mat(T_got, T_ref))
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',l=', l, ',nb=', nb
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, A_got, B_ref, B_got, T_ref, T_got)
    end do
    call report_finalize()
end program test_dtpqrt
