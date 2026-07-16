program test_zggbak
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zggbak
    use ref_quad_lapack, only: zggbal, zggbak
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer, parameter :: ncols = 4
    integer :: i, n, info, ilo, ihi
    complex(ep), allocatable :: A(:,:), B(:,:), V0(:,:), V_ref(:,:), V_got(:,:)
    real(ep), allocatable :: ls(:), rs(:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zggbak', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 310071 + 47 * i)
        call gen_matrix_complex(n, n, B, seed = 310081 + 47 * i)
        call gen_matrix_complex(n, ncols, V0, seed = 310091 + 47 * i)
        allocate(ls(n), rs(n), work(6*n))
        call zggbal('B', n, A, n, B, n, ilo, ihi, ls, rs, work, info)
        allocate(V_ref(n, ncols), V_got(n, ncols))
        V_ref = V0; V_got = V0
        call zggbak('B', 'R', n, ilo, ihi, ls, rs, ncols, V_ref, n, info)
        call target_zggbak('B', 'R', n, ilo, ihi, ls, rs, ncols, V_got, n, info)
        err = max_rel_err_mat_z(V_got, V_ref)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, V0, ls, rs, work, V_ref, V_got)
    end do
    call report_finalize()
end program test_zggbak
