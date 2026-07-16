program test_zunmrq
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgerqf, target_zunmrq
    use ref_quad_lapack, only: zgerqf, zunmrq
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ks(*) = [8, 16]
    integer, parameter :: nrhs  = 4
    character(len=1), parameter :: transes(2) = ['N', 'C']
    integer :: i, m, k, info, lwork, jt
    complex(ep), allocatable :: A0(:,:), C0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: tau_ref(:), tau_got(:), C_ref(:,:), C_got(:,:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zunmrq', target_name)
    do i = 1, size(ms)
        m = ms(i); k = ks(i)
        call gen_matrix_complex(k, m,    A0, seed = 71001 + 47 * i)
        call gen_matrix_complex(m, nrhs, C0, seed = 71011 + 47 * i)
        allocate(A_ref(k,m), A_got(k,m), tau_ref(k), tau_got(k))
        A_ref = A0; A_got = A0
        call zgerqf(k, m, A_ref, k, tau_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgerqf(k, m, A_ref, k, tau_ref, work, lwork, info)
        deallocate(work)
        call target_zgerqf(k, m, A_got, k, tau_got, info)
        do jt = 1, size(transes)
            allocate(C_ref(m,nrhs), C_got(m,nrhs))
            C_ref = C0; C_got = C0
            call zunmrq('L', transes(jt), m, nrhs, k, A_ref, k, tau_ref, &
                        C_ref, m, wopt, -1, info)
            lwork = max(1, int(real(wopt(1), ep)))
            allocate(work(lwork))
            call zunmrq('L', transes(jt), m, nrhs, k, A_ref, k, tau_ref, &
                        C_ref, m, work, lwork, info)
            deallocate(work)
            call target_zunmrq('L', transes(jt), m, nrhs, k, A_got, k, tau_got, &
                               C_got, m, info)
            err = max_rel_err_mat_z(C_got, C_ref)
            tol = 16.0_ep * real(m, ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'trans=', transes(jt), ',m=', m, ',k=', k
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end do
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_zunmrq
