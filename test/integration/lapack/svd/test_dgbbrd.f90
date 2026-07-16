! dgbbrd: banded bidiagonal reduction.
program test_dgbbrd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgbbrd
    use ref_quad_lapack, only: dgbbrd
    implicit none

    integer, parameter :: ms(*) = [16, 24]
    integer, parameter :: ns(*) = [16, 24]
    integer, parameter :: kls(*) = [2, 3]
    integer, parameter :: kus(*) = [2, 3]
    integer :: i, m, n, kl, ku, ldab, info, j, k
    real(ep), allocatable :: A(:,:), AB_r(:,:), AB_g(:,:)
    real(ep), allocatable :: D_r(:), D_g(:), E_r(:), E_g(:)
    real(ep), allocatable :: Q(:,:), Pt(:,:), C(:,:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgbbrd', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); kl = kls(i); ku = kus(i)
        ldab = kl + ku + 1
        call gen_matrix_quad(m, n, A, seed = 25201 + 67 * i)
        ! Pack as banded
        allocate(AB_r(ldab, n), AB_g(ldab, n))
        AB_r = 0.0_ep
        do j = 1, n
            do k = max(1, j - ku), min(m, j + kl)
                AB_r(ku + 1 + k - j, j) = A(k, j)
            end do
        end do
        AB_g = AB_r
        allocate(D_r(min(m, n)), D_g(min(m, n)), E_r(max(1, min(m, n)-1)), E_g(max(1, min(m, n)-1)))
        allocate(Q(m, m), Pt(n, n), C(1, 1), work(2*max(m, n)))
        call dgbbrd('N', m, n, 0, kl, ku, AB_r, ldab, D_r, E_r, Q, m, Pt, n, C, 1, work, info)
        call target_dgbbrd('N', m, n, 0, kl, ku, AB_g, ldab, D_g, E_g, Q, m, Pt, n, C, 1, info)
        err = max(max_rel_err_vec(D_g, D_r), max_rel_err_vec(E_g, E_r))
        tol = 64.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',kl=', kl
        call report_case(trim(label), err, tol)
        deallocate(A, AB_r, AB_g, D_r, D_g, E_r, E_g, Q, Pt, C, work)
    end do
    call report_finalize()
end program test_dgbbrd
