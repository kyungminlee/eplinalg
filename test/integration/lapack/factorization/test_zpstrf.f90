! zpstrf: pivoted Cholesky of a Hermitian PD matrix.
program test_zpstrf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hpd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zpstrf
    use ref_quad_lapack, only: zpstrf
    implicit none

    integer, parameter :: ns(*) = [8, 24, 48]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, u, n, info_ref, info_got, rank_ref, rank_got
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep),    allocatable :: work(:)
    integer, allocatable :: piv_ref(:), piv_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('zpstrf', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hpd_matrix_quad(n, A0, seed = 21351 + 71 * i)
        do u = 1, size(uplos)
            allocate(A_ref(n, n), A_got(n, n), piv_ref(n), piv_got(n), work(2*n))
            A_ref = A0; A_got = A0
            call zpstrf(uplos(u), n, A_ref, n, piv_ref, rank_ref, -1.0_ep, work, info_ref)
            call target_zpstrf(uplos(u), n, A_got, n, piv_got, rank_got, -1.0_ep, info_got)
            err = max_rel_err_mat_z(A_got, A_ref)
            if (any(piv_ref(1:n) /= piv_got(1:n)) .or. rank_ref /= rank_got) then
                err = max(err, 1.0_ep)
            end if
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(u), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, piv_ref, piv_got, work)
        end do
        deallocate(A0)
    end do
    call report_finalize()
end program test_zpstrf
