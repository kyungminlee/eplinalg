program test_zgetrs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgetrf, target_zgetrs
    use ref_quad_lapack, only: zgetrf, zgetrs
    implicit none

    integer, parameter :: ns(*)   = [16, 32, 48]
    integer, parameter :: nrhs    = 2
    character(len=1), parameter :: transes(*) = ['N', 'C', 'T']
    integer :: i, n, info, ti
    complex(ep), allocatable :: A(:,:), B0(:,:), B_ref(:,:), B_got(:,:), Aref(:,:), Agot(:,:)
    integer,     allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgetrs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n,    A,  seed = 5201 + 23 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 5211 + 23 * i)
        do ti = 1, size(transes)
            allocate(Aref(n, n), Agot(n, n), B_ref(n, nrhs), B_got(n, nrhs))
            allocate(ipiv_ref(n), ipiv_got(n))
            Aref = A; Agot = A
            call zgetrf(n, n, Aref, n, ipiv_ref, info)
            call zgetrf(n, n, Agot, n, ipiv_got, info)
            B_ref = B0; B_got = B0
            call zgetrs(transes(ti), n, nrhs, Aref, n, ipiv_ref, B_ref, n, info)
            call target_zgetrs(transes(ti), n, nrhs, Agot, n, ipiv_got, B_got, n, info)
            err = max_rel_err_mat_z(B_got, B_ref)
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'trans=', transes(ti), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(Aref, Agot, B_ref, B_got, ipiv_ref, ipiv_got)
        end do
    end do
    call report_finalize()
end program test_zgetrs
