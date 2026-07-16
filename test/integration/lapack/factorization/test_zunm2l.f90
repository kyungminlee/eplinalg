! zunm2l: apply Q from QL (unblocked, complex).
program test_zunm2l
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgeqlf, target_zunm2l
    use ref_quad_lapack, only: zgeqlf, zunm2l
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [8, 16]
    integer, parameter :: nrhs  = 4
    character(len=1), parameter :: transes(2) = ['N', 'C']
    integer :: i, m, n, info, lwork, jt
    complex(ep), allocatable :: A0(:,:), C0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: tau_ref(:), tau_got(:), C_ref(:,:), C_got(:,:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zunm2l', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n,    A0, seed = 19851 + 61 * i)
        call gen_matrix_complex(m, nrhs, C0, seed = 19861 + 61 * i)
        allocate(A_ref(m,n), A_got(m,n), tau_ref(n), tau_got(n))
        A_ref = A0; A_got = A0
        call zgeqlf(m, n, A_ref, m, tau_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgeqlf(m, n, A_ref, m, tau_ref, work, lwork, info)
        deallocate(work)
        call target_zgeqlf(m, n, A_got, m, tau_got, info)
        do jt = 1, size(transes)
            allocate(C_ref(m,nrhs), C_got(m,nrhs), work(nrhs))
            C_ref = C0; C_got = C0
            call zunm2l('L', transes(jt), m, nrhs, n, A_ref, m, tau_ref, C_ref, m, work, info)
            call target_zunm2l('L', transes(jt), m, nrhs, n, A_got, m, tau_got, C_got, m, info)
            err = max_rel_err_mat_z(C_got, C_ref)
            tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'trans=', transes(jt), ',m=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got, work)
        end do
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_zunm2l
