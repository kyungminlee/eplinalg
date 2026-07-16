program test_dtrtrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtrtrs
    use ref_quad_lapack, only: dtrtrs
    implicit none

    integer, parameter :: ns(*) = [8, 32, 64]
    integer, parameter :: nrhs  = 2
    character(len=1), parameter :: uplos(2)  = ['U', 'L']
    character(len=1), parameter :: transes(2) = ['N', 'T']
    character(len=1), parameter :: diags(2)  = ['N', 'U']
    integer :: i, n, info, ju, jt, jd, j
    real(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('dtrtrs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A0, seed = 25001 + 47 * i)
        do j = 1, n
            A0(j, j) = A0(j, j) + real(2*n, ep)
        end do
        call gen_matrix_quad(n, nrhs, B0, seed = 25011 + 47 * i)
        do ju = 1, size(uplos)
            do jt = 1, size(transes)
                do jd = 1, size(diags)
                    allocate(A_ref(n,n), B_ref(n,nrhs), B_got(n,nrhs))
                    A_ref = A0; B_ref = B0; B_got = B0
                    call dtrtrs(uplos(ju), transes(jt), diags(jd), n, nrhs, &
                                A_ref, n, B_ref, n, info)
                    call target_dtrtrs(uplos(ju), transes(jt), diags(jd), n, nrhs, &
                                A_ref, n, B_got, n, info)
                    err = max_rel_err_mat(B_got, B_ref)
                    tol = 16.0_ep * real(n, ep)**2 * target_eps
                    write(label, '(a,a,a,a,a,a,a,i0)') 'uplo=', uplos(ju), &
                        ',trans=', transes(jt), ',diag=', diags(jd), ',n=', n
                    call report_case(trim(label), err, tol)
                    deallocate(A_ref, B_ref, B_got)
                end do
            end do
        end do
    end do
    call report_finalize()
end program test_dtrtrs
