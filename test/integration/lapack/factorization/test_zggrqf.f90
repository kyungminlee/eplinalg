program test_zggrqf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, max_rel_err_vec_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zggrqf
    use ref_quad_lapack, only: zggrqf
    implicit none

    integer, parameter :: ms(*) = [12, 20, 28]
    integer, parameter :: ps(*) = [8,  16, 24]
    integer, parameter :: ns(*) = [16, 24, 32]
    integer :: i, m, p, n, info, lwork
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:)
    complex(ep), allocatable :: taua_ref(:), taub_ref(:), taua_got(:), taub_got(:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zggrqf', target_name)
    do i = 1, size(ms)
        m = ms(i); p = ps(i); n = ns(i)
        call gen_matrix_complex(m, n, A0, seed = 320061 + 47 * i)
        call gen_matrix_complex(p, n, B0, seed = 320071 + 47 * i)
        allocate(A_ref(m,n), B_ref(p,n), A_got(m,n), B_got(p,n))
        allocate(taua_ref(min(m,n)), taub_ref(min(p,n)), &
                 taua_got(min(m,n)), taub_got(min(p,n)))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        call zggrqf(m, p, n, A_ref, m, taua_ref, B_ref, p, taub_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zggrqf(m, p, n, A_ref, m, taua_ref, B_ref, p, taub_ref, work, lwork, info)
        call target_zggrqf(m, p, n, A_got, m, taua_got, B_got, p, taub_got, info)
        err = max(max_rel_err_mat_z(A_got, A_ref), max_rel_err_mat_z(B_got, B_ref), &
                  max_rel_err_vec_z(taua_got, taua_ref), max_rel_err_vec_z(taub_got, taub_ref))
        tol = 16.0_ep * real(max(m,p,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',p=', p, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, A_got, B_got, taua_ref, taub_ref, &
                   taua_got, taub_got, work)
    end do
    call report_finalize()
end program test_zggrqf
