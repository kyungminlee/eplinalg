program test_zgtsv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_vector_complex, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgtsv
    use ref_quad_lapack, only: zgtsv
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer, parameter :: nrhs = 2
    integer :: i, n, info, j
    complex(ep), allocatable :: dl0(:), d0(:), du0(:), B0(:,:)
    complex(ep), allocatable :: dl_ref(:), d_ref(:), du_ref(:)
    complex(ep), allocatable :: dl_got(:), d_got(:), du_got(:)
    complex(ep), allocatable :: B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgtsv', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_complex(n-1, dl0, seed = 100001 + 47 * i)
        call gen_vector_complex(n,   d0,  seed = 100011 + 47 * i)
        call gen_vector_complex(n-1, du0, seed = 100021 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 100031 + 47 * i)
        do j = 1, n; d0(j) = d0(j) + cmplx(real(4, ep), 0.0_ep, ep); end do
        allocate(dl_ref(n-1), d_ref(n), du_ref(n-1))
        allocate(dl_got(n-1), d_got(n), du_got(n-1))
        allocate(B_ref(n, nrhs), B_got(n, nrhs))
        dl_ref = dl0; d_ref = d0; du_ref = du0
        dl_got = dl0; d_got = d0; du_got = du0
        B_ref = B0; B_got = B0
        call zgtsv(n, nrhs, dl_ref, d_ref, du_ref, B_ref, n, info)
        call target_zgtsv(n, nrhs, dl_got, d_got, du_got, B_got, n, info)
        err = max_rel_err_mat_z(B_got, B_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(dl0, d0, du0, B0, dl_ref, d_ref, du_ref, dl_got, d_got, du_got, B_ref, B_got)
    end do
    call report_finalize()
end program test_zgtsv
