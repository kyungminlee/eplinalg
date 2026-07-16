! zpbstf: split-Cholesky factor of a Hermitian PD banded matrix.
program test_zpbstf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hpd_matrix_quad, pack_herm_band_quad
    use target_lapack,   only: target_name, target_eps, target_zpbstf
    use ref_quad_lapack, only: zpbstf
    implicit none

    integer, parameter :: ns(*)  = [16, 32]
    integer, parameter :: kds(*) = [3,  5]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, u, n, kd, info, ldab
    complex(ep), allocatable :: A(:,:), AB0(:,:), AB_ref(:,:), AB_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zpbstf', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_hpd_matrix_quad(n, A, seed = 21451 + 79 * i)
        do u = 1, size(uplos)
            allocate(AB0(ldab, n), AB_ref(ldab, n), AB_got(ldab, n))
            call pack_herm_band_quad(uplos(u), n, kd, A, AB0)
            AB_ref = AB0; AB_got = AB0
            call zpbstf(uplos(u), n, kd, AB_ref, ldab, info)
            call target_zpbstf(uplos(u), n, kd, AB_got, ldab, info)
            err = max_rel_err_mat_z(AB_got, AB_ref)
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'uplo=', uplos(u), ',n=', n, ',kd=', kd
            call report_case(trim(label), err, tol)
            deallocate(AB0, AB_ref, AB_got)
        end do
        deallocate(A)
    end do
    call report_finalize()
end program test_zpbstf
