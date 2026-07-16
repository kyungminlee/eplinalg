program test_zhpevx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zhpevx
    use ref_quad_lapack, only: zhpevx
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    character(len=1), parameter :: jobzs(2) = ['N', 'V']
    integer :: i, n, info, jz, np, mref, mgot
    integer, allocatable :: iwork(:), ifail(:)
    complex(ep), allocatable :: A0(:,:), AP_ref(:), AP_got(:)
    complex(ep), allocatable :: z_ref(:,:), z_got(:,:), work(:)
    real(ep),    allocatable :: w_ref(:), w_got(:), rwork(:)
    real(ep) :: abstol, vl, vu, err, tol
    character(len=48) :: label

    call report_init('zhpevx', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_hermitian_matrix_quad(n, A0, seed = 58001 + 47 * i)
        abstol = 0.0_ep; vl = 0.0_ep; vu = 0.0_ep
        do jz = 1, size(jobzs)
            allocate(AP_ref(np), AP_got(np))
            allocate(w_ref(n), w_got(n), z_ref(n, n), z_got(n, n))
            allocate(work(max(1, 2*n)), rwork(7*n), iwork(5*n), ifail(n))
            call pack_herm_packed_quad('U', n, A0, AP_ref)
            AP_got = AP_ref
            call zhpevx(jobzs(jz), 'A', 'U', n, AP_ref, vl, vu, 1, n, &
                        abstol, mref, w_ref, z_ref, n, &
                        work, rwork, iwork, ifail, info)
            call target_zhpevx(jobzs(jz), 'A', 'U', n, AP_got, vl, vu, 1, n, &
                               abstol, mgot, w_got, z_got, n, ifail, info)
            err = max_rel_err_vec(w_got(1:mgot), w_ref(1:mref))
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'jobz=', jobzs(jz), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(AP_ref, AP_got, w_ref, w_got, z_ref, z_got, &
                       work, rwork, iwork, ifail)
        end do
    end do
    call report_finalize()
end program test_zhpevx
