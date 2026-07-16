! dtfttr: rectangular full packed (RFP) -> triangular full.
program test_dtfttr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtfttr
    use ref_quad_lapack, only: dtrttf, dtfttr
    implicit none

    integer, parameter :: ns(*) = [7, 16, 33]
    character(len=1), parameter :: transrs(2) = ['N', 'T']
    character(len=1), parameter :: uplos(2)   = ['U', 'L']
    integer :: i, t, u, n, info, info_got, nt
    real(ep), allocatable :: A0(:,:), ARF(:), A_ref(:,:), A_got(:,:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('dtfttr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        nt = n*(n+1)/2
        call gen_matrix_quad(n, n, A0, seed = 17301 + 101 * i)
        do t = 1, size(transrs)
            do u = 1, size(uplos)
                allocate(ARF(nt), A_ref(n, n), A_got(n, n))
                call dtrttf(transrs(t), uplos(u), n, A0, n, ARF, info)
                A_ref = 0.0_ep; A_got = 0.0_ep
                call dtfttr(transrs(t), uplos(u), n, ARF, A_ref, n, info)
                call target_dtfttr(transrs(t), uplos(u), n, ARF, A_got, n, info_got)
                err = max_rel_err_mat(A_got, A_ref)
                tol = 16.0_ep * real(n, ep) * target_eps
                write(label, '(a,a,a,a,a,i0)') 'transr=', transrs(t), &
                    ',uplo=', uplos(u), ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(ARF, A_ref, A_got)
            end do
        end do
        deallocate(A0)
    end do
    call report_finalize()
end program test_dtfttr
