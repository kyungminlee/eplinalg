! dgebak applies the balance from dgebal back to a matrix of
! eigenvectors. Feed both implementations the reference balance scale
! and ilo/ihi to isolate dgebak from dgebal's variance.
program test_dgebak
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgebak
    use ref_quad_lapack, only: dgebal, dgebak
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer, parameter :: ncols = 4
    integer :: i, n, info, ilo, ihi
    real(ep), allocatable :: A(:,:), V0(:,:), V_ref(:,:), V_got(:,:), scale(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgebak', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 300031 + 47 * i)
        call gen_matrix_quad(n, ncols, V0, seed = 300041 + 47 * i)
        allocate(scale(n))
        call dgebal('B', n, A, n, ilo, ihi, scale, info)
        allocate(V_ref(n, ncols), V_got(n, ncols))
        V_ref = V0; V_got = V0
        call dgebak('B', 'R', n, ilo, ihi, scale, ncols, V_ref, n, info)
        call target_dgebak('B', 'R', n, ilo, ihi, scale, ncols, V_got, n, info)
        err = max_rel_err_mat(V_got, V_ref)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, V0, scale, V_ref, V_got)
    end do
    call report_finalize()
end program test_dgebak
