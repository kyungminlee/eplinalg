program test_zgebal
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec, max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgebal
    use ref_quad_lapack, only: zgebal
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, ilo_ref, ihi_ref, ilo_got, ihi_got
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: s_ref(:), s_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgebal', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A0, seed = 300021 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), s_ref(n), s_got(n))
        A_ref = A0; A_got = A0
        call zgebal('B', n, A_ref, n, ilo_ref, ihi_ref, s_ref, info)
        call target_zgebal('B', n, A_got, n, ilo_got, ihi_got, s_got, info)
        err = max(max_rel_err_mat_z(A_got, A_ref), max_rel_err_vec(s_got, s_ref))
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, s_ref, s_got)
    end do
    call report_finalize()
end program test_zgebal
