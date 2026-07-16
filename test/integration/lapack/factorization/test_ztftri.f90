! ztftri: triangular inverse where A is in RFP form (complex).
program test_ztftri
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec_z
    use test_data,       only: gen_hpd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_ztftri
    use ref_quad_lapack, only: ztrttf, ztftri
    implicit none

    integer, parameter :: ns(*) = [7, 16, 33]
    character(len=1), parameter :: transrs(2) = ['N', 'C']
    character(len=1), parameter :: uplos(2)   = ['U', 'L']
    character(len=1), parameter :: diags(2)   = ['N', 'U']
    integer :: i, t, u, d, n, info, nt
    complex(ep), allocatable :: A(:,:), ARF0(:), AR_ref(:), AR_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztftri', target_name)
    do i = 1, size(ns)
        n = ns(i); nt = n*(n+1)/2
        call gen_hpd_matrix_quad(n, A, seed = 21751 + 83 * i)
        do t = 1, size(transrs)
            do u = 1, size(uplos)
                allocate(ARF0(nt))
                call ztrttf(transrs(t), uplos(u), n, A, n, ARF0, info)
                do d = 1, size(diags)
                    allocate(AR_ref(nt), AR_got(nt))
                    AR_ref = ARF0; AR_got = ARF0
                    call ztftri(transrs(t), uplos(u), diags(d), n, AR_ref, info)
                    call target_ztftri(transrs(t), uplos(u), diags(d), n, AR_got, info)
                    err = max_rel_err_vec_z(AR_got, AR_ref)
                    tol = 16.0_ep * real(n, ep)**3 * target_eps
                    write(label, '(a,a,a,a,a,a,a,i0)') 'tR=', transrs(t), &
                        ',u=', uplos(u), ',d=', diags(d), ',n=', n
                    call report_case(trim(label), err, tol)
                    deallocate(AR_ref, AR_got)
                end do
                deallocate(ARF0)
            end do
        end do
        deallocate(A)
    end do
    call report_finalize()
end program test_ztftri
