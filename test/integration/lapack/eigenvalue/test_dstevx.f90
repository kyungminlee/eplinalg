program test_dstevx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_dstevx
    use ref_quad_lapack, only: dstevx
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    character(len=1), parameter :: jobzs(2) = ['N', 'V']
    integer :: i, n, info, jz, mref, mgot
    integer, allocatable :: iwork(:), ifail(:)
    real(ep), allocatable :: d0(:), e0(:), d_ref(:), e_ref(:), d_got(:), e_got(:)
    real(ep), allocatable :: w_ref(:), w_got(:), z_ref(:,:), z_got(:,:), work(:)
    real(ep) :: abstol, vl, vu, err, tol
    character(len=48) :: label

    call report_init('dstevx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n,   d0, seed = 45001 + 47 * i)
        call gen_vector_quad(n-1, e0, seed = 45011 + 47 * i)
        abstol = 0.0_ep; vl = 0.0_ep; vu = 0.0_ep
        do jz = 1, size(jobzs)
            allocate(d_ref(n), e_ref(n-1), d_got(n), e_got(n-1))
            allocate(w_ref(n), w_got(n), z_ref(n,n), z_got(n,n))
            allocate(work(5*n), iwork(5*n), ifail(n))
            d_ref = d0; e_ref = e0; d_got = d0; e_got = e0
            call dstevx(jobzs(jz), 'A', n, d_ref, e_ref, vl, vu, 1, n, &
                        abstol, mref, w_ref, z_ref, n, work, iwork, ifail, info)
            call target_dstevx(jobzs(jz), 'A', n, d_got, e_got, vl, vu, 1, n, &
                               abstol, mgot, w_got, z_got, n, ifail, info)
            err = max_rel_err_vec(w_got(1:mgot), w_ref(1:mref))
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'jobz=', jobzs(jz), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(d_ref, e_ref, d_got, e_got, w_ref, w_got, &
                       z_ref, z_got, work, iwork, ifail)
        end do
        deallocate(d0, e0)
    end do
    call report_finalize()
end program test_dstevx
