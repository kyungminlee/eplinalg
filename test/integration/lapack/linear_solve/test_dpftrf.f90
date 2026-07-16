! dpftrf: Cholesky factor of an SPD matrix in RFP storage.
program test_dpftrf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_spd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dpftrf
    use ref_quad_lapack, only: dtrttf, dpftrf
    implicit none

    integer, parameter :: ns(*) = [7, 16, 33]
    character(len=1), parameter :: transrs(2) = ['N', 'T']
    character(len=1), parameter :: uplos(2)   = ['U', 'L']
    integer :: i, t, u, n, info, nt
    real(ep), allocatable :: A(:,:), ARF0(:), AR_ref(:), AR_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('dpftrf', target_name)
    do i = 1, size(ns)
        n = ns(i); nt = n*(n+1)/2
        call gen_spd_matrix_quad(n, A, seed = 21001 + 89 * i)
        do t = 1, size(transrs)
            do u = 1, size(uplos)
                allocate(ARF0(nt), AR_ref(nt), AR_got(nt))
                call dtrttf(transrs(t), uplos(u), n, A, n, ARF0, info)
                AR_ref = ARF0; AR_got = ARF0
                call dpftrf(transrs(t), uplos(u), n, AR_ref, info)
                call target_dpftrf(transrs(t), uplos(u), n, AR_got, info)
                err = max_rel_err_vec(AR_got, AR_ref)
                tol = 16.0_ep * real(n, ep)**2 * target_eps
                write(label, '(a,a,a,a,a,i0)') 'transr=', transrs(t), &
                    ',uplo=', uplos(u), ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(ARF0, AR_ref, AR_got)
            end do
        end do
        deallocate(A)
    end do
    call report_finalize()
end program test_dpftrf
