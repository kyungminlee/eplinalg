program test_dggqrf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dggqrf
    use ref_quad_lapack, only: dggqrf
    implicit none

    integer, parameter :: ns(*) = [16, 24, 32]
    integer, parameter :: ms(*) = [12, 20, 28]
    integer, parameter :: ps(*) = [8,  16, 24]
    integer :: i, n, m, p, info, lwork
    real(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:)
    real(ep), allocatable :: taua_ref(:), taub_ref(:), taua_got(:), taub_got(:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dggqrf', target_name)
    do i = 1, size(ns)
        n = ns(i); m = ms(i); p = ps(i)
        call gen_matrix_quad(n, m, A0, seed = 320001 + 47 * i)
        call gen_matrix_quad(n, p, B0, seed = 320011 + 47 * i)
        allocate(A_ref(n,m), B_ref(n,p), A_got(n,m), B_got(n,p))
        allocate(taua_ref(min(n,m)), taub_ref(min(n,p)), &
                 taua_got(min(n,m)), taub_got(min(n,p)))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        call dggqrf(n, m, p, A_ref, n, taua_ref, B_ref, n, taub_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dggqrf(n, m, p, A_ref, n, taua_ref, B_ref, n, taub_ref, work, lwork, info)
        call target_dggqrf(n, m, p, A_got, n, taua_got, B_got, n, taub_got, info)
        err = max(max_rel_err_mat(A_got, A_ref), max_rel_err_mat(B_got, B_ref), &
                  max_rel_err_vec(taua_got, taua_ref), max_rel_err_vec(taub_got, taub_ref))
        tol = 16.0_ep * real(max(n,m,p), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'n=', n, ',m=', m, ',p=', p
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, A_got, B_got, taua_ref, taub_ref, &
                   taua_got, taub_got, work)
    end do
    call report_finalize()
end program test_dggqrf
