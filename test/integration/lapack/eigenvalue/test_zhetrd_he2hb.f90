program test_zhetrd_he2hb
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zhetrd_he2hb
    use ref_quad_lapack, only: zhetrd_he2hb
    implicit none
    integer, parameter :: ns(*) = [32, 64]
    integer, parameter :: kd    = 4
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, info, ju, lwork, ldab, j
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), AB_r(:,:), AB_g(:,:), tau_r(:), tau_g(:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('zhetrd_he2hb', target_name)
    do i = 1, size(ns)
        n = ns(i)
        ldab = kd + 1
        call gen_hermitian_matrix_quad(n, A0, seed = 111101 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + cmplx(real(n, ep), 0.0_ep, ep); end do
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n), AB_r(ldab, n), AB_g(ldab, n), tau_r(max(1,n-kd)), tau_g(max(1,n-kd)))
            A_ref = A0; A_got = A0
            call zhetrd_he2hb(uplos(ju), n, kd, A_ref, n, AB_r, ldab, tau_r, wopt, -1, info)
            lwork = max(1, int(real(wopt(1), ep)))
            allocate(work(lwork))
            call zhetrd_he2hb(uplos(ju), n, kd, A_ref, n, AB_r, ldab, tau_r, work, lwork, info)
            deallocate(work)
            call target_zhetrd_he2hb(uplos(ju), n, kd, A_got, n, AB_g, ldab, tau_g, info)
            err = real(abs(info), ep)
            tol = 0.5_ep
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, AB_r, AB_g, tau_r, tau_g)
        end do
    end do
    call report_finalize()
end program test_zhetrd_he2hb
