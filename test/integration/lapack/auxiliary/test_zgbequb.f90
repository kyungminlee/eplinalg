program test_zgbequb
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec, rel_err_scalar
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgbequb
    use ref_quad_lapack, only: zgbequb
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kl    = 2
    integer, parameter :: ku    = 3
    integer :: i, n, ldab, info, j, k
    complex(ep), allocatable :: A(:,:), AB(:,:)
    real(ep), allocatable :: R_ref(:), C_ref(:), R_got(:), C_got(:)
    real(ep) :: rcnd_ref, ccnd_ref, am_ref, rcnd_got, ccnd_got, am_got
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgbequb', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = kl + ku + 1
        call gen_matrix_complex(n, n, A, seed = 210111 + 47 * i)
        allocate(AB(ldab, n))
        AB = (0.0_ep, 0.0_ep)
        do j = 1, n
            do k = max(1, j-ku), min(n, j+kl)
                AB(ku + 1 + k - j, j) = A(k, j)
            end do
        end do
        allocate(R_ref(n), C_ref(n), R_got(n), C_got(n))
        call zgbequb(n, n, kl, ku, AB, ldab, R_ref, C_ref, rcnd_ref, ccnd_ref, am_ref, info)
        call target_zgbequb(n, n, kl, ku, AB, ldab, R_got, C_got, rcnd_got, ccnd_got, am_got, info)
        err = max(max_rel_err_vec(R_got, R_ref), max_rel_err_vec(C_got, C_ref), &
                  rel_err_scalar(rcnd_got, rcnd_ref), rel_err_scalar(ccnd_got, ccnd_ref), &
                  rel_err_scalar(am_got, am_ref))
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AB, R_ref, C_ref, R_got, C_got)
    end do
    call report_finalize()
end program test_zgbequb
