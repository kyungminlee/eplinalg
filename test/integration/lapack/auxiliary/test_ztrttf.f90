! ztrttf: triangular full -> rectangular full packed (RFP, complex).
program test_ztrttf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztrttf
    use ref_quad_lapack, only: ztrttf
    implicit none

    integer, parameter :: ns(*) = [7, 16, 33]
    character(len=1), parameter :: transrs(2) = ['N', 'C']
    character(len=1), parameter :: uplos(2)   = ['U', 'L']
    integer :: i, t, u, n, info_ref, info_got, nt
    complex(ep), allocatable :: A(:,:), ARF_ref(:), ARF_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztrttf', target_name)
    do i = 1, size(ns)
        n = ns(i)
        nt = n*(n+1)/2
        call gen_matrix_complex(n, n, A, seed = 17251 + 97 * i)
        do t = 1, size(transrs)
            do u = 1, size(uplos)
                allocate(ARF_ref(nt), ARF_got(nt))
                ARF_ref = (0.0_ep, 0.0_ep); ARF_got = (0.0_ep, 0.0_ep)
                call ztrttf(transrs(t), uplos(u), n, A, n, ARF_ref, info_ref)
                call target_ztrttf(transrs(t), uplos(u), n, A, n, ARF_got, info_got)
                err = max_rel_err_vec_z(ARF_got, ARF_ref)
                tol = 16.0_ep * real(n, ep) * target_eps
                write(label, '(a,a,a,a,a,i0)') 'transr=', transrs(t), &
                    ',uplo=', uplos(u), ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(ARF_ref, ARF_got)
            end do
        end do
        deallocate(A)
    end do
    call report_finalize()
end program test_ztrttf
