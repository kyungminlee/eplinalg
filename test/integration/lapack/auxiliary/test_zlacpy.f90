program test_zlacpy
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zlacpy
    use ref_quad_lapack, only: zlacpy
    implicit none

    integer, parameter :: ms(*) = [8, 32, 96]
    integer, parameter :: ns(*) = [6, 20, 64]
    character(len=1), parameter :: uplos(3) = ['U', 'L', 'A']
    integer :: i, k, m, n
    complex(ep), allocatable :: A(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zlacpy', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A, seed = 14101 + 89 * i)
        do k = 1, size(uplos)
            allocate(B_ref(m, n), B_got(m, n))
            B_ref = (0.0_ep, 0.0_ep); B_got = (0.0_ep, 0.0_ep)
            call zlacpy(uplos(k), m, n, A, m, B_ref, m)
            call target_zlacpy(uplos(k), m, n, A, m, B_got, m)
            err = max_rel_err_mat_z(B_got, B_ref)
            tol = 16.0_ep * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'uplo=', uplos(k), ',m=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(B_ref, B_got)
        end do
    end do
    call report_finalize()
end program test_zlacpy
