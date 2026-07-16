program test_zgghrd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgghrd
    use ref_quad_lapack, only: zgghrd
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:)
    complex(ep), allocatable :: Q(:,:), Z(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgghrd', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A0, seed = 310121 + 47 * i)
        call gen_matrix_complex(n, n, B0, seed = 310131 + 47 * i)
        allocate(A_ref(n,n), B_ref(n,n), A_got(n,n), B_got(n,n), Q(1,1), Z(1,1))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        call zgghrd('N', 'N', n, 1, n, A_ref, n, B_ref, n, Q, 1, Z, 1, info)
        call target_zgghrd('N', 'N', n, 1, n, A_got, n, B_got, n, Q, 1, Z, 1, info)
        err = max(max_rel_err_mat_z(A_got, A_ref), max_rel_err_mat_z(B_got, B_ref))
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, A_got, B_got, Q, Z)
    end do
    call report_finalize()
end program test_zgghrd
