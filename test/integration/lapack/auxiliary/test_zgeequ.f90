program test_zgeequ
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec, rel_err_scalar
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgeequ
    use ref_quad_lapack, only: zgeequ
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info
    complex(ep), allocatable :: A(:,:)
    real(ep), allocatable :: R_ref(:), C_ref(:), R_got(:), C_got(:)
    real(ep) :: rcnd_ref, ccnd_ref, am_ref, rcnd_got, ccnd_got, am_got
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgeequ', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 200101 + 47 * i)
        allocate(R_ref(n), C_ref(n), R_got(n), C_got(n))
        call zgeequ(n, n, A, n, R_ref, C_ref, rcnd_ref, ccnd_ref, am_ref, info)
        call target_zgeequ(n, n, A, n, R_got, C_got, rcnd_got, ccnd_got, am_got, info)
        err = max(max_rel_err_vec(R_got, R_ref), max_rel_err_vec(C_got, C_ref), &
                  rel_err_scalar(rcnd_got, rcnd_ref), rel_err_scalar(ccnd_got, ccnd_ref), &
                  rel_err_scalar(am_got, am_ref))
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, R_ref, C_ref, R_got, C_got)
    end do
    call report_finalize()
end program test_zgeequ
