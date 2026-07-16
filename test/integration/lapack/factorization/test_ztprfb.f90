! ztprfb: block reflector applied via pentagonal V (complex).
program test_ztprfb
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztprfb
    use ref_quad_lapack, only: ztpqrt, ztprfb
    implicit none

    integer, parameter :: nrhs = 4
    integer, parameter :: mps(*) = [4, 8]
    integer, parameter :: nps(*) = [3, 6]
    character(len=1), parameter :: transes(2) = ['N', 'C']
    integer :: i, mp, np, info, jt, nb
    complex(ep), allocatable :: A0(:,:), B0(:,:), V(:,:), T_arr(:,:)
    complex(ep), allocatable :: Aap0(:,:), Bap0(:,:), Aap_ref(:,:), Bap_ref(:,:)
    complex(ep), allocatable :: Aap_got(:,:), Bap_got(:,:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztprfb', target_name)
    do i = 1, size(mps)
        mp = mps(i); np = nps(i); nb = np
        call gen_matrix_complex(np, np, A0, seed = 22451 + 71 * i)
        block
            integer :: j
            do j = 1, np-1
                A0(j+1:np, j) = (0.0_ep, 0.0_ep)
            end do
        end block
        call gen_matrix_complex(mp, np, B0, seed = 22461 + 71 * i)
        allocate(V(mp, np), T_arr(nb, np))
        V = B0; T_arr = (0.0_ep, 0.0_ep)
        block
            complex(ep), allocatable :: work(:), A_t(:,:)
            allocate(work(nb*np), A_t(np, np))
            A_t = A0
            call ztpqrt(mp, np, 0, nb, A_t, np, V, mp, T_arr, nb, work, info)
            deallocate(work, A_t)
        end block
        call gen_matrix_complex(np, nrhs, Aap0, seed = 22471 + 71 * i)
        call gen_matrix_complex(mp, nrhs, Bap0, seed = 22481 + 71 * i)
        do jt = 1, size(transes)
            allocate(Aap_ref(np, nrhs), Bap_ref(mp, nrhs), &
                     Aap_got(np, nrhs), Bap_got(mp, nrhs))
            Aap_ref = Aap0; Aap_got = Aap0
            Bap_ref = Bap0; Bap_got = Bap0
            block
                complex(ep), allocatable :: work2(:,:)
                allocate(work2(np, nrhs))
                call ztprfb('L', transes(jt), 'F', 'C', mp, nrhs, np, 0, &
                            V, mp, T_arr, nb, Aap_ref, np, Bap_ref, mp, work2, np)
                deallocate(work2)
            end block
            call target_ztprfb('L', transes(jt), 'F', 'C', mp, nrhs, np, 0, &
                               V, mp, T_arr, nb, Aap_got, np, Bap_got, mp, np)
            err = max(max_rel_err_mat_z(Aap_got, Aap_ref), &
                      max_rel_err_mat_z(Bap_got, Bap_ref))
            tol = 16.0_ep * real(max(mp, np, nrhs), ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'trans=', transes(jt), ',mp=', mp, ',np=', np
            call report_case(trim(label), err, tol)
            deallocate(Aap_ref, Bap_ref, Aap_got, Bap_got)
        end do
        deallocate(A0, B0, V, T_arr, Aap0, Bap0)
    end do
    call report_finalize()
end program test_ztprfb
