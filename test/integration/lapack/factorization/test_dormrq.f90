! dormrq: apply Q from an RQ factorization to a matrix C.
! Reflectors of order m are produced by dgerqf on a k-by-m matrix
! (k <= m), applied to an m-by-nrhs C from the left.
program test_dormrq
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgerqf, target_dormrq
    use ref_quad_lapack, only: dgerqf, dormrq
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ks(*) = [8, 16]
    integer, parameter :: nrhs  = 4
    character(len=1), parameter :: transes(2) = ['N', 'T']
    integer :: i, m, k, info, lwork, jt
    real(ep), allocatable :: A0(:,:), C0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: tau_ref(:), tau_got(:), C_ref(:,:), C_got(:,:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dormrq', target_name)
    do i = 1, size(ms)
        m = ms(i); k = ks(i)
        call gen_matrix_quad(k, m,    A0, seed = 68001 + 47 * i)
        call gen_matrix_quad(m, nrhs, C0, seed = 68011 + 47 * i)
        allocate(A_ref(k,m), A_got(k,m), tau_ref(k), tau_got(k))
        A_ref = A0; A_got = A0
        call dgerqf(k, m, A_ref, k, tau_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgerqf(k, m, A_ref, k, tau_ref, work, lwork, info)
        deallocate(work)
        call target_dgerqf(k, m, A_got, k, tau_got, info)
        do jt = 1, size(transes)
            allocate(C_ref(m,nrhs), C_got(m,nrhs))
            C_ref = C0; C_got = C0
            call dormrq('L', transes(jt), m, nrhs, k, A_ref, k, tau_ref, &
                        C_ref, m, wopt, -1, info)
            lwork = max(1, int(wopt(1)))
            allocate(work(lwork))
            call dormrq('L', transes(jt), m, nrhs, k, A_ref, k, tau_ref, &
                        C_ref, m, work, lwork, info)
            deallocate(work)
            call target_dormrq('L', transes(jt), m, nrhs, k, A_got, k, tau_got, &
                               C_got, m, info)
            err = max_rel_err_mat(C_got, C_ref)
            tol = 16.0_ep * real(m, ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'trans=', transes(jt), ',m=', m, ',k=', k
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end do
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_dormrq
