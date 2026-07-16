! dormrz: apply Q from RZ factorization (dtzrzf). Use the ref factorization
! for both ref- and target-apply so the apply alone is exercised.
program test_dormrz
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dormrz
    use ref_quad_lapack, only: dtzrzf, dormrz
    implicit none

    integer, parameter :: ms(*) = [6, 16]
    integer, parameter :: ns(*) = [12, 32]
    integer, parameter :: nrhs  = 4
    character(len=1), parameter :: transes(2) = ['N', 'T']
    integer :: i, m, n, info, lwork, jt, k, l, j
    real(ep), allocatable :: A0(:,:), C0(:,:), A_fact(:,:)
    real(ep), allocatable :: tau(:), C_ref(:,:), C_got(:,:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dormrz', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        k = m            ! number of reflectors
        l = n - m        ! length of trapezoidal "extra" region
        call gen_matrix_quad(m, n,    A0, seed = 20201 + 97 * i)
        call gen_matrix_quad(n, nrhs, C0, seed = 20211 + 97 * i)
        do j = 1, m-1
            A0(j+1:m, j) = 0.0_ep
        end do
        allocate(A_fact(m,n), tau(m))
        A_fact = A0
        call dtzrzf(m, n, A_fact, m, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dtzrzf(m, n, A_fact, m, tau, work, lwork, info)
        deallocate(work)
        do jt = 1, size(transes)
            allocate(C_ref(n,nrhs), C_got(n,nrhs))
            C_ref = C0; C_got = C0
            call dormrz('L', transes(jt), n, nrhs, k, l, A_fact, m, tau, &
                        C_ref, n, wopt, -1, info)
            lwork = max(1, int(wopt(1))); allocate(work(lwork))
            call dormrz('L', transes(jt), n, nrhs, k, l, A_fact, m, tau, &
                        C_ref, n, work, lwork, info)
            deallocate(work)
            call target_dormrz('L', transes(jt), n, nrhs, k, l, A_fact, m, tau, &
                               C_got, n, info)
            err = max_rel_err_mat(C_got, C_ref)
            tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'trans=', transes(jt), ',m=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end do
        deallocate(A_fact, tau)
    end do
    call report_finalize()
end program test_dormrz
