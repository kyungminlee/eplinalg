program test_dggrqf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dggrqf
    use ref_quad_lapack, only: dggrqf
    implicit none

    integer, parameter :: ms(*) = [12, 20, 28]
    integer, parameter :: ps(*) = [8,  16, 24]
    integer, parameter :: ns(*) = [16, 24, 32]
    integer :: i, m, p, n, info, lwork
    real(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:)
    real(ep), allocatable :: taua_ref(:), taub_ref(:), taua_got(:), taub_got(:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dggrqf', target_name)
    do i = 1, size(ms)
        m = ms(i); p = ps(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 320041 + 47 * i)
        call gen_matrix_quad(p, n, B0, seed = 320051 + 47 * i)
        allocate(A_ref(m,n), B_ref(p,n), A_got(m,n), B_got(p,n))
        allocate(taua_ref(min(m,n)), taub_ref(min(p,n)), &
                 taua_got(min(m,n)), taub_got(min(p,n)))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        call dggrqf(m, p, n, A_ref, m, taua_ref, B_ref, p, taub_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dggrqf(m, p, n, A_ref, m, taua_ref, B_ref, p, taub_ref, work, lwork, info)
        call target_dggrqf(m, p, n, A_got, m, taua_got, B_got, p, taub_got, info)
        err = max(max_rel_err_mat(A_got, A_ref), max_rel_err_mat(B_got, B_ref), &
                  max_rel_err_vec(taua_got, taua_ref), max_rel_err_vec(taub_got, taub_ref))
        tol = 16.0_ep * real(max(m,p,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',p=', p, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, A_got, B_got, taua_ref, taub_ref, &
                   taua_got, taub_got, work)
    end do
    call report_finalize()
end program test_dggrqf
