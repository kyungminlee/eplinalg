program test_zsytrf_aa_2stage
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, max_rel_err_vec_z
    use test_data,       only: gen_complex_symmetric_quad
    use target_lapack,   only: target_name, target_eps, target_zsytrf_aa_2stage
    use ref_quad_lapack, only: zsytrf_aa_2stage
    implicit none
    integer, parameter :: ns(*) = [16, 32, 64]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, info, ju, lwork, ltb, j
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), TB_ref(:), TB_got(:), work(:)
    integer,     allocatable :: ipiv_ref(:), ipiv_got(:), ipiv2_ref(:), ipiv2_got(:)
    complex(ep) :: wopt(1), tbopt(1)
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('zsytrf_aa_2stage', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_complex_symmetric_quad(n, A0, seed = 54001 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + cmplx(real(n, ep), 0.0_ep, ep); end do
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n), ipiv2_ref(n), ipiv2_got(n))
            A_ref = A0; A_got = A0
            call zsytrf_aa_2stage(uplos(ju), n, A_ref, n, tbopt, -1, ipiv_ref, ipiv2_ref, wopt, -1, info)
            ltb = max(1, int(real(tbopt(1), ep))); lwork = max(1, int(real(wopt(1), ep)))
            allocate(TB_ref(ltb), TB_got(ltb), work(lwork))
            TB_ref = (0.0_ep, 0.0_ep); TB_got = (0.0_ep, 0.0_ep)
            call zsytrf_aa_2stage(uplos(ju), n, A_ref, n, TB_ref, ltb, ipiv_ref, ipiv2_ref, work, lwork, info)
            deallocate(work)
            call target_zsytrf_aa_2stage(uplos(ju), n, A_got, n, TB_got, ltb, ipiv_got, ipiv2_got, info)
            err = real(abs(info), ep)
            tol = 0.5_ep
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, ipiv_ref, ipiv_got, ipiv2_ref, ipiv2_got, TB_ref, TB_got)
        end do
    end do
    call report_finalize()
end program test_zsytrf_aa_2stage
