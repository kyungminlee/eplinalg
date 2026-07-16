program test_zggbal
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec, max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zggbal
    use ref_quad_lapack, only: zggbal
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, ilo_ref, ihi_ref, ilo_got, ihi_got
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:)
    real(ep), allocatable :: ls_ref(:), rs_ref(:), ls_got(:), rs_got(:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zggbal', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A0, seed = 310021 + 47 * i)
        call gen_matrix_complex(n, n, B0, seed = 310031 + 47 * i)
        allocate(A_ref(n,n), B_ref(n,n), A_got(n,n), B_got(n,n))
        allocate(ls_ref(n), rs_ref(n), ls_got(n), rs_got(n), work(6*n))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        call zggbal('B', n, A_ref, n, B_ref, n, ilo_ref, ihi_ref, ls_ref, rs_ref, work, info)
        call target_zggbal('B', n, A_got, n, B_got, n, ilo_got, ihi_got, ls_got, rs_got, info)
        err = max(max_rel_err_mat_z(A_got, A_ref), max_rel_err_mat_z(B_got, B_ref), &
                  max_rel_err_vec(ls_got, ls_ref), max_rel_err_vec(rs_got, rs_ref))
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, A_got, B_got, ls_ref, rs_ref, ls_got, rs_got, work)
    end do
    call report_finalize()
end program test_zggbal
