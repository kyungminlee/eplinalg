! ztpqrt2: unblocked pentagonal QR factorization of [A; B] (complex).
program test_ztpqrt2
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztpqrt2
    use ref_quad_lapack, only: ztpqrt2
    implicit none

    integer, parameter :: ms(*) = [12, 24]
    integer, parameter :: ns(*) = [6,  12]
    integer, parameter :: ls(*) = [0,  3]
    integer :: i, m, n, l, j, info
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: B_ref(:,:), B_got(:,:), T_ref(:,:), T_got(:,:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztpqrt2', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); l = ls(i)
        call gen_matrix_complex(n, n, A0, seed = 22351 + 79 * i)
        do j = 1, n-1
            A0(j+1:n, j) = (0.0_ep, 0.0_ep)
        end do
        call gen_matrix_complex(m, n, B0, seed = 22361 + 79 * i)
        allocate(A_ref(n, n), A_got(n, n), B_ref(m, n), B_got(m, n), &
                 T_ref(n, n), T_got(n, n))
        A_ref = A0; A_got = A0; B_ref = B0; B_got = B0
        T_ref = (0.0_ep, 0.0_ep); T_got = (0.0_ep, 0.0_ep)
        call ztpqrt2(m, n, l, A_ref, n, B_ref, m, T_ref, n, info)
        call target_ztpqrt2(m, n, l, A_got, n, B_got, m, T_got, n, info)
        err = max(max_rel_err_mat_z(A_got, A_ref), max_rel_err_mat_z(B_got, B_ref), &
                  max_rel_err_mat_z(T_got, T_ref))
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',l=', l
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, A_got, B_ref, B_got, T_ref, T_got)
    end do
    call report_finalize()
end program test_ztpqrt2
