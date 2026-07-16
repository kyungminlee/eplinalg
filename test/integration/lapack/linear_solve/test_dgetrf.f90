program test_dgetrf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgetrf
    use ref_quad_lapack, only: dgetrf
    implicit none

    integer, parameter :: ms(*) = [8, 32, 64]
    integer, parameter :: ns(*) = [8, 24, 80]
    integer :: i, m, n, info_ref, info_got
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgetrf', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 2001 + 37 * i)

        allocate(A_ref(m, n), A_got(m, n))
        allocate(ipiv_ref(min(m, n)), ipiv_got(min(m, n)))
        A_ref = A0;  A_got = A0

        call dgetrf(m, n, A_ref, m, ipiv_ref, info_ref)
        call target_dgetrf(m, n, A_got, m, ipiv_got, info_got)

        err = max_rel_err_mat(A_got, A_ref)
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)

        deallocate(A_ref, A_got, ipiv_ref, ipiv_got)
    end do
    call report_finalize()
end program test_dgetrf
