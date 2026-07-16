program test_dposv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_spd_matrix_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dposv
    use ref_quad_lapack, only: dposv
    implicit none

    integer, parameter :: ns(*) = [8, 32, 64]
    integer, parameter :: nrhs  = 2
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, info, ju
    real(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), A_got(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dposv', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_spd_matrix_quad(n, A0, seed = 14001 + 47 * i)
        call gen_matrix_quad(n, nrhs, B0, seed = 14011 + 47 * i)
        do ju = 1, size(uplos)
            allocate(A_ref(n,n), A_got(n,n), B_ref(n,nrhs), B_got(n,nrhs))
            A_ref = A0; A_got = A0; B_ref = B0; B_got = B0
            call dposv(uplos(ju), n, nrhs, A_ref, n, B_ref, n, info)
            call target_dposv(uplos(ju), n, nrhs, A_got, n, B_got, n, info)
            err = max_rel_err_mat(B_got, B_ref)
            tol = 16.0_ep * real(n, ep)**3 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, B_ref, B_got)
        end do
    end do
    call report_finalize()
end program test_dposv
