! dtpmqrt: apply Q from pentagonal QR.
program test_dtpmqrt
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtpmqrt
    use ref_quad_lapack, only: dtpqrt, dtpmqrt
    implicit none

    integer, parameter :: nrhs = 4
    integer, parameter :: mps(*) = [4, 8]    ! pentagonal B rows
    integer, parameter :: nps(*) = [3, 6]    ! cols
    integer, parameter :: nbs(*) = [2, 3]
    character(len=1), parameter :: transes(2) = ['N', 'T']
    integer :: i, mp, np, nb, info, jt
    real(ep), allocatable :: A0(:,:), B0(:,:), V(:,:), T_arr(:,:)
    real(ep), allocatable :: Aap0(:,:), Bap0(:,:), Aap_ref(:,:), Bap_ref(:,:)
    real(ep), allocatable :: Aap_got(:,:), Bap_got(:,:)
    real(ep) :: err, tol
    character(len=64) :: label
    character(len=1)  :: side

    call report_init('dtpmqrt', target_name)
    side = 'L'
    do i = 1, size(mps)
        mp = mps(i); np = nps(i); nb = nbs(i)
        ! Build factorization with reference: tpqrt(mp, np, 0, nb, A, np, B, mp, T, nb)
        call gen_matrix_quad(np, np, A0, seed = 22201 + 67 * i)
        block
            integer :: j
            do j = 1, np-1
                A0(j+1:np, j) = 0.0_ep
            end do
        end block
        call gen_matrix_quad(mp, np, B0, seed = 22211 + 67 * i)
        allocate(V(mp, np), T_arr(nb, np))
        V = B0; T_arr = 0.0_ep
        block
            real(ep), allocatable :: work(:), A_t(:,:)
            allocate(work(nb*np), A_t(np, np))
            A_t = A0
            call dtpqrt(mp, np, 0, nb, A_t, np, V, mp, T_arr, nb, work, info)
            deallocate(work, A_t)
        end block
        ! Now apply Q to (Aap; Bap) where Aap is np x nrhs, Bap is mp x nrhs
        call gen_matrix_quad(np, nrhs, Aap0, seed = 22221 + 67 * i)
        call gen_matrix_quad(mp, nrhs, Bap0, seed = 22231 + 67 * i)
        do jt = 1, size(transes)
            allocate(Aap_ref(np, nrhs), Bap_ref(mp, nrhs), &
                     Aap_got(np, nrhs), Bap_got(mp, nrhs))
            Aap_ref = Aap0; Aap_got = Aap0
            Bap_ref = Bap0; Bap_got = Bap0
            block
                real(ep), allocatable :: work(:)
                allocate(work(nb*max(mp, nrhs)))
                call dtpmqrt(side, transes(jt), mp, nrhs, np, 0, nb, V, mp, &
                             T_arr, nb, Aap_ref, np, Bap_ref, mp, work, info)
                deallocate(work)
            end block
            call target_dtpmqrt(side, transes(jt), mp, nrhs, np, 0, nb, V, mp, &
                                T_arr, nb, Aap_got, np, Bap_got, mp, info)
            err = max(max_rel_err_mat(Aap_got, Aap_ref), max_rel_err_mat(Bap_got, Bap_ref))
            tol = 16.0_ep * real(max(mp, np, nrhs), ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'trans=', transes(jt), ',mp=', mp, ',np=', np
            call report_case(trim(label), err, tol)
            deallocate(Aap_ref, Bap_ref, Aap_got, Bap_got)
        end do
        deallocate(A0, B0, V, T_arr, Aap0, Bap0)
    end do
    call report_finalize()
end program test_dtpmqrt
