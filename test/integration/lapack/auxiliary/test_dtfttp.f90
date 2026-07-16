! dtfttp: rectangular full packed (RFP) -> triangular packed.
program test_dtfttp
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtfttp
    use ref_quad_lapack, only: dtrttf, dtfttp
    implicit none

    integer, parameter :: ns(*) = [7, 16, 33]
    character(len=1), parameter :: transrs(2) = ['N', 'T']
    character(len=1), parameter :: uplos(2)   = ['U', 'L']
    integer :: i, t, u, n, info, info_got, nt
    real(ep), allocatable :: A0(:,:), ARF(:), AP_ref(:), AP_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('dtfttp', target_name)
    do i = 1, size(ns)
        n = ns(i)
        nt = n*(n+1)/2
        call gen_matrix_quad(n, n, A0, seed = 17501 + 113 * i)
        do u = 1, size(uplos)
            do t = 1, size(transrs)
                allocate(ARF(nt), AP_ref(nt), AP_got(nt))
                call dtrttf(transrs(t), uplos(u), n, A0, n, ARF, info)
                AP_ref = 0.0_ep; AP_got = 0.0_ep
                call dtfttp(transrs(t), uplos(u), n, ARF, AP_ref, info)
                call target_dtfttp(transrs(t), uplos(u), n, ARF, AP_got, info_got)
                err = max_rel_err_vec(AP_got, AP_ref)
                tol = 16.0_ep * real(n, ep) * target_eps
                write(label, '(a,a,a,a,a,i0)') 'transr=', transrs(t), &
                    ',uplo=', uplos(u), ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(ARF, AP_ref, AP_got)
            end do
        end do
        deallocate(A0)
    end do
    call report_finalize()
end program test_dtfttp
