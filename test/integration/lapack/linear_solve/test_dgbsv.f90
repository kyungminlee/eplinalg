program test_dgbsv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgbsv
    use ref_quad_lapack, only: dgbsv
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs  = 2
    integer, parameter :: kl    = 2
    integer, parameter :: ku    = 3
    integer :: i, n, ldab, info, j, k
    real(ep), allocatable :: Adense(:,:), AB0(:,:), B0(:,:)
    real(ep), allocatable :: AB_ref(:,:), AB_got(:,:), B_ref(:,:), B_got(:,:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgbsv', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = 2*kl + ku + 1
        call gen_matrix_quad(n, n, Adense, seed = 62001 + 47 * i)
        do j = 1, n
            Adense(j, j) = Adense(j, j) + real(2*n, ep)
        end do
        call gen_matrix_quad(n, nrhs, B0, seed = 62011 + 47 * i)
        allocate(AB0(ldab, n))
        AB0 = 0.0_ep
        do j = 1, n
            do k = max(1, j-ku), min(n, j+kl)
                AB0(kl + ku + 1 + k - j, j) = Adense(k, j)
            end do
        end do
        allocate(AB_ref(ldab,n), AB_got(ldab,n), B_ref(n,nrhs), B_got(n,nrhs))
        allocate(ipiv_ref(n), ipiv_got(n))
        AB_ref = AB0; AB_got = AB0; B_ref = B0; B_got = B0
        call dgbsv(n, kl, ku, nrhs, AB_ref, ldab, ipiv_ref, B_ref, n, info)
        call target_dgbsv(n, kl, ku, nrhs, AB_got, ldab, ipiv_got, B_got, n, info)
        err = max_rel_err_mat(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(Adense, AB0, B0, AB_ref, AB_got, B_ref, B_got, ipiv_ref, ipiv_got)
    end do
    call report_finalize()
end program test_dgbsv
