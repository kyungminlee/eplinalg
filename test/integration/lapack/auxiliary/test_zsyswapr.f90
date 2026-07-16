program test_zsyswapr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_complex_symmetric_quad
    use target_lapack,   only: target_name, target_eps, target_zsyswapr
    use ref_quad_lapack, only: zsyswapr
    implicit none
    integer, parameter :: ns(*) = [10, 32, 64]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, ju, i1, i2
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('zsyswapr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_complex_symmetric_quad(n, A0, seed = 71201 + 19 * i)
        i1 = 2; i2 = n - 1
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n)); A_ref = A0; A_got = A0
            call zsyswapr(uplos(ju), n, A_ref, n, i1, i2)
            call target_zsyswapr(uplos(ju), n, A_got, n, i1, i2)
            err = max_rel_err_mat_z(A_got, A_ref)
            tol = 16.0_ep * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got)
        end do
    end do
    call report_finalize()
end program test_zsyswapr
