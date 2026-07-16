! zpftrs: solve A*X = B from Hermitian RFP Cholesky factor.
program test_zpftrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hpd_matrix_quad, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zpftrs
    use ref_quad_lapack, only: ztrttf, zpftrf, zpftrs
    implicit none

    integer, parameter :: ns(*) = [7, 16, 33]
    integer, parameter :: nrhs = 4
    character(len=1), parameter :: transrs(2) = ['N', 'C']
    character(len=1), parameter :: uplos(2)   = ['U', 'L']
    integer :: i, t, u, n, info, nt
    complex(ep), allocatable :: A(:,:), B0(:,:), ARF(:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('zpftrs', target_name)
    do i = 1, size(ns)
        n = ns(i); nt = n*(n+1)/2
        call gen_hpd_matrix_quad(n, A, seed = 21251 + 73 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 21261 + 73 * i)
        do t = 1, size(transrs)
            do u = 1, size(uplos)
                allocate(ARF(nt), B_ref(n, nrhs), B_got(n, nrhs))
                call ztrttf(transrs(t), uplos(u), n, A, n, ARF, info)
                call zpftrf(transrs(t), uplos(u), n, ARF, info)
                B_ref = B0; B_got = B0
                call zpftrs(transrs(t), uplos(u), n, nrhs, ARF, B_ref, n, info)
                call target_zpftrs(transrs(t), uplos(u), n, nrhs, ARF, B_got, n, info)
                err = max_rel_err_mat_z(B_got, B_ref)
                tol = 16.0_ep * real(n, ep)**2 * target_eps
                write(label, '(a,a,a,a,a,i0)') 'transr=', transrs(t), &
                    ',uplo=', uplos(u), ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(ARF, B_ref, B_got)
            end do
        end do
        deallocate(A, B0)
    end do
    call report_finalize()
end program test_zpftrs
