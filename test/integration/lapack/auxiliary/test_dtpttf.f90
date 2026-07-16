! dtpttf: triangular packed -> rectangular full packed (RFP).
program test_dtpttf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtpttf
    use ref_quad_lapack, only: dtrttp, dtpttf
    implicit none

    integer, parameter :: ns(*) = [7, 16, 33]
    character(len=1), parameter :: transrs(2) = ['N', 'T']
    character(len=1), parameter :: uplos(2)   = ['U', 'L']
    integer :: i, t, u, n, info, info_got, nt
    real(ep), allocatable :: A0(:,:), AP(:), ARF_ref(:), ARF_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('dtpttf', target_name)
    do i = 1, size(ns)
        n = ns(i)
        nt = n*(n+1)/2
        call gen_matrix_quad(n, n, A0, seed = 17401 + 107 * i)
        do u = 1, size(uplos)
            allocate(AP(nt))
            call dtrttp(uplos(u), n, A0, n, AP, info)
            do t = 1, size(transrs)
                allocate(ARF_ref(nt), ARF_got(nt))
                ARF_ref = 0.0_ep; ARF_got = 0.0_ep
                call dtpttf(transrs(t), uplos(u), n, AP, ARF_ref, info)
                call target_dtpttf(transrs(t), uplos(u), n, AP, ARF_got, info_got)
                err = max_rel_err_vec(ARF_got, ARF_ref)
                tol = 16.0_ep * real(n, ep) * target_eps
                write(label, '(a,a,a,a,a,i0)') 'transr=', transrs(t), &
                    ',uplo=', uplos(u), ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(ARF_ref, ARF_got)
            end do
            deallocate(AP)
        end do
        deallocate(A0)
    end do
    call report_finalize()
end program test_dtpttf
