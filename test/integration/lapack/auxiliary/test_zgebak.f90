program test_zgebak
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgebak
    use ref_quad_lapack, only: zgebal, zgebak
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer, parameter :: ncols = 4
    integer :: i, n, info, ilo, ihi
    complex(ep), allocatable :: A(:,:), V0(:,:), V_ref(:,:), V_got(:,:)
    real(ep), allocatable :: scale(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgebak', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 300051 + 47 * i)
        call gen_matrix_complex(n, ncols, V0, seed = 300061 + 47 * i)
        allocate(scale(n))
        call zgebal('B', n, A, n, ilo, ihi, scale, info)
        allocate(V_ref(n, ncols), V_got(n, ncols))
        V_ref = V0; V_got = V0
        call zgebak('B', 'R', n, ilo, ihi, scale, ncols, V_ref, n, info)
        call target_zgebak('B', 'R', n, ilo, ihi, scale, ncols, V_got, n, info)
        err = max_rel_err_mat_z(V_got, V_ref)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, V0, scale, V_ref, V_got)
    end do
    call report_finalize()
end program test_zgebak
