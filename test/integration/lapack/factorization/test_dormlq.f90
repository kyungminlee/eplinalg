! dormlq: apply Q from LQ.
program test_dormlq
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgelqf, target_dormlq
    use ref_quad_lapack, only: dgelqf, dormlq
    implicit none

    integer, parameter :: ms(*) = [8, 16]
    integer, parameter :: ns(*) = [16, 32]
    integer, parameter :: nrhs  = 4
    character(len=1), parameter :: transes(2) = ['N', 'T']
    integer :: i, m, n, info, lwork, jt
    real(ep), allocatable :: A0(:,:), C0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: tau_ref(:), tau_got(:), C_ref(:,:), C_got(:,:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dormlq', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n,    A0, seed = 19901 + 67 * i)
        call gen_matrix_quad(n, nrhs, C0, seed = 19911 + 67 * i)
        allocate(A_ref(m,n), A_got(m,n), tau_ref(m), tau_got(m))
        A_ref = A0; A_got = A0
        call dgelqf(m, n, A_ref, m, tau_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgelqf(m, n, A_ref, m, tau_ref, work, lwork, info)
        deallocate(work)
        call target_dgelqf(m, n, A_got, m, tau_got, info)
        do jt = 1, size(transes)
            allocate(C_ref(n,nrhs), C_got(n,nrhs))
            C_ref = C0; C_got = C0
            call dormlq('L', transes(jt), n, nrhs, m, A_ref, m, tau_ref, C_ref, n, wopt, -1, info)
            lwork = max(1, int(wopt(1))); allocate(work(lwork))
            call dormlq('L', transes(jt), n, nrhs, m, A_ref, m, tau_ref, C_ref, n, work, lwork, info)
            deallocate(work)
            call target_dormlq('L', transes(jt), n, nrhs, m, A_got, m, tau_got, C_got, n, info)
            err = max_rel_err_mat(C_got, C_ref)
            tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'trans=', transes(jt), ',m=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end do
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_dormlq
