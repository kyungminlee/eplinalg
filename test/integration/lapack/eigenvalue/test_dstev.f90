program test_dstev
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_dstev
    use ref_quad_lapack, only: dstev
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    character(len=1), parameter :: jobzs(2) = ['N', 'V']
    integer :: i, n, info, jz
    real(ep), allocatable :: d0(:), e0(:), d_ref(:), e_ref(:), d_got(:), e_got(:)
    real(ep), allocatable :: z_ref(:,:), z_got(:,:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dstev', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n,   d0, seed = 43001 + 47 * i)
        call gen_vector_quad(n-1, e0, seed = 43011 + 47 * i)
        do jz = 1, size(jobzs)
            allocate(d_ref(n), e_ref(n-1), d_got(n), e_got(n-1))
            allocate(z_ref(n,n), z_got(n,n), work(max(1, 2*n - 2)))
            d_ref = d0; e_ref = e0; d_got = d0; e_got = e0
            call dstev(jobzs(jz), n, d_ref, e_ref, z_ref, n, work, info)
            call target_dstev(jobzs(jz), n, d_got, e_got, z_got, n, info)
            err = max_rel_err_vec(d_got, d_ref)
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'jobz=', jobzs(jz), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(d_ref, e_ref, d_got, e_got, z_ref, z_got, work)
        end do
        deallocate(d0, e0)
    end do
    call report_finalize()
end program test_dstev
