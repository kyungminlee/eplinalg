program test_zpbtrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hpd_matrix_quad, gen_matrix_complex, pack_herm_band_quad
    use target_lapack,   only: target_name, target_eps, target_zpbtrf, target_zpbtrs
    use ref_quad_lapack, only: zpbtrf, zpbtrs
    implicit none

    integer, parameter :: ns(*)  = [16, 32, 64]
    integer, parameter :: kds(*) = [2, 4, 8]
    integer, parameter :: nrhs = 3
    integer :: i, n, kd, ldab, info
    complex(ep), allocatable :: A0(:,:), AB_ref(:,:), AB_got(:,:)
    complex(ep), allocatable :: B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zpbtrs', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_hpd_matrix_quad(n, A0, seed = 96001 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 96011 + 47 * i)
        allocate(AB_ref(ldab, n), AB_got(ldab, n), B_ref(n, nrhs), B_got(n, nrhs))
        call pack_herm_band_quad('U', n, kd, A0, AB_ref); AB_got = AB_ref
        call zpbtrf('U', n, kd, AB_ref, ldab, info)
        call target_zpbtrf('U', n, kd, AB_got, ldab, info)
        B_ref = B0; B_got = B0
        call zpbtrs('U', n, kd, nrhs, AB_ref, ldab, B_ref, n, info)
        call target_zpbtrs('U', n, kd, nrhs, AB_got, ldab, B_got, n, info)
        err = max_rel_err_mat_z(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'n=', n, ',kd=', kd
        call report_case(trim(label), err, tol)
        deallocate(AB_ref, AB_got, B_ref, B_got)
    end do
    call report_finalize()
end program test_zpbtrs
