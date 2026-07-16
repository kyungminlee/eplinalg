! zunmrz: apply Q from RZ factorization (ztzrzf, complex).
program test_zunmrz
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zunmrz
    use ref_quad_lapack, only: ztzrzf, zunmrz
    implicit none

    integer, parameter :: ms(*) = [6, 16]
    integer, parameter :: ns(*) = [12, 32]
    integer, parameter :: nrhs  = 4
    character(len=1), parameter :: transes(2) = ['N', 'C']
    integer :: i, m, n, info, lwork, jt, k, l, j
    complex(ep), allocatable :: A0(:,:), C0(:,:), A_fact(:,:)
    complex(ep), allocatable :: tau(:), C_ref(:,:), C_got(:,:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zunmrz', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        k = m
        l = n - m
        call gen_matrix_complex(m, n,    A0, seed = 20251 + 101 * i)
        call gen_matrix_complex(n, nrhs, C0, seed = 20261 + 101 * i)
        do j = 1, m-1
            A0(j+1:m, j) = (0.0_ep, 0.0_ep)
        end do
        allocate(A_fact(m,n), tau(m))
        A_fact = A0
        call ztzrzf(m, n, A_fact, m, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call ztzrzf(m, n, A_fact, m, tau, work, lwork, info)
        deallocate(work)
        do jt = 1, size(transes)
            allocate(C_ref(n,nrhs), C_got(n,nrhs))
            C_ref = C0; C_got = C0
            call zunmrz('L', transes(jt), n, nrhs, k, l, A_fact, m, tau, &
                        C_ref, n, wopt, -1, info)
            lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
            call zunmrz('L', transes(jt), n, nrhs, k, l, A_fact, m, tau, &
                        C_ref, n, work, lwork, info)
            deallocate(work)
            call target_zunmrz('L', transes(jt), n, nrhs, k, l, A_fact, m, tau, &
                               C_got, n, info)
            err = max_rel_err_mat_z(C_got, C_ref)
            tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'trans=', transes(jt), ',m=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end do
        deallocate(A_fact, tau)
    end do
    call report_finalize()
end program test_zunmrz
