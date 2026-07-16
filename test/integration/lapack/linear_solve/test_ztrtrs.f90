program test_ztrtrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztrtrs
    use ref_quad_lapack, only: ztrtrs
    implicit none

    integer, parameter :: ns(*) = [8, 32, 64]
    integer, parameter :: nrhs  = 2
    character(len=1), parameter :: uplos(2)  = ['U', 'L']
    character(len=1), parameter :: transes(3) = ['N', 'T', 'C']
    character(len=1), parameter :: diags(2)  = ['N', 'U']
    integer :: i, n, info, ju, jt, jd, j
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztrtrs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A0, seed = 25101 + 47 * i)
        do j = 1, n
            A0(j, j) = A0(j, j) + cmplx(real(2*n, ep), 0.0_ep, ep)
        end do
        call gen_matrix_complex(n, nrhs, B0, seed = 25111 + 47 * i)
        do ju = 1, size(uplos)
            do jt = 1, size(transes)
                do jd = 1, size(diags)
                    allocate(A_ref(n,n), B_ref(n,nrhs), B_got(n,nrhs))
                    A_ref = A0; B_ref = B0; B_got = B0
                    call ztrtrs(uplos(ju), transes(jt), diags(jd), n, nrhs, &
                                A_ref, n, B_ref, n, info)
                    call target_ztrtrs(uplos(ju), transes(jt), diags(jd), n, nrhs, &
                                A_ref, n, B_got, n, info)
                    err = max_rel_err_mat_z(B_got, B_ref)
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
end program test_ztrtrs
