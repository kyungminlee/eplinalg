! zunmhr: apply Q from Hessenberg reduction (complex).
program test_zunmhr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgehrd, target_zunmhr
    use ref_quad_lapack, only: zgehrd, zunmhr
    implicit none

    integer, parameter :: ns(*) = [12, 32]
    integer, parameter :: nrhs  = 4
    character(len=1), parameter :: transes(2) = ['N', 'C']
    integer :: i, n, info, lwork, jt, ilo, ihi
    complex(ep), allocatable :: A0(:,:), C0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: tau_ref(:), tau_got(:), C_ref(:,:), C_got(:,:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zunmhr', target_name)
    do i = 1, size(ns)
        n = ns(i); ilo = 1; ihi = n
        call gen_matrix_complex(n, n,    A0, seed = 20051 + 79 * i)
        call gen_matrix_complex(n, nrhs, C0, seed = 20061 + 79 * i)
        allocate(A_ref(n,n), A_got(n,n), tau_ref(n-1), tau_got(n-1))
        A_ref = A0; A_got = A0
        call zgehrd(n, ilo, ihi, A_ref, n, tau_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgehrd(n, ilo, ihi, A_ref, n, tau_ref, work, lwork, info)
        deallocate(work)
        call target_zgehrd(n, ilo, ihi, A_got, n, tau_got, info)
        do jt = 1, size(transes)
            allocate(C_ref(n,nrhs), C_got(n,nrhs))
            C_ref = C0; C_got = C0
            call zunmhr('L', transes(jt), n, nrhs, ilo, ihi, A_ref, n, tau_ref, &
                        C_ref, n, wopt, -1, info)
            lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
            call zunmhr('L', transes(jt), n, nrhs, ilo, ihi, A_ref, n, tau_ref, &
                        C_ref, n, work, lwork, info)
            deallocate(work)
            call target_zunmhr('L', transes(jt), n, nrhs, ilo, ihi, A_got, n, tau_got, &
                               C_got, n, info)
            err = max_rel_err_mat_z(C_got, C_ref)
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'trans=', transes(jt), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end do
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_zunmhr
