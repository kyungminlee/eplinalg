! dtpmlqt: apply Q from pentagonal LQ.
program test_dtpmlqt
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtpmlqt
    use ref_quad_lapack, only: dtplqt, dtpmlqt
    implicit none

    integer, parameter :: nrhs = 4
    integer, parameter :: mps(*) = [3, 6]
    integer, parameter :: nps(*) = [4, 8]
    integer, parameter :: mbs(*) = [2, 3]
    character(len=1), parameter :: transes(2) = ['N', 'T']
    integer :: i, mp, np, mb, info, jt
    real(ep), allocatable :: A0(:,:), B0(:,:), V(:,:), T_arr(:,:)
    real(ep), allocatable :: Aap0(:,:), Bap0(:,:), Aap_ref(:,:), Bap_ref(:,:)
    real(ep), allocatable :: Aap_got(:,:), Bap_got(:,:)
    real(ep) :: err, tol
    character(len=64) :: label
    character(len=1)  :: side

    call report_init('dtpmlqt', target_name)
    side = 'L'
    do i = 1, size(mps)
        mp = mps(i); np = nps(i); mb = mbs(i)
        call gen_matrix_quad(mp, mp, A0, seed = 22301 + 73 * i)
        block
            integer :: j
            do j = 2, mp
                A0(1:j-1, j) = 0.0_ep
            end do
        end block
        call gen_matrix_quad(mp, np, B0, seed = 22311 + 73 * i)
        allocate(V(mp, np), T_arr(mb, mp))
        V = B0; T_arr = 0.0_ep
        block
            real(ep), allocatable :: work(:), A_t(:,:)
            allocate(work(mb*mp), A_t(mp, mp))
            A_t = A0
            call dtplqt(mp, np, 0, mb, A_t, mp, V, mp, T_arr, mb, work, info)
            deallocate(work, A_t)
        end block
        ! Apply: K = mp reflectors, M_arg = mp (B's rows), N_arg = nrhs
        call gen_matrix_quad(mp, nrhs, Aap0, seed = 22321 + 73 * i)
        call gen_matrix_quad(np, nrhs, Bap0, seed = 22331 + 73 * i)
        do jt = 1, size(transes)
            allocate(Aap_ref(mp, nrhs), Bap_ref(np, nrhs), &
                     Aap_got(mp, nrhs), Bap_got(np, nrhs))
            Aap_ref = Aap0; Aap_got = Aap0
            Bap_ref = Bap0; Bap_got = Bap0
            block
                real(ep), allocatable :: work(:)
                allocate(work(mb*max(mp, nrhs)))
                call dtpmlqt(side, transes(jt), np, nrhs, mp, 0, mb, V, mp, &
                             T_arr, mb, Aap_ref, mp, Bap_ref, np, work, info)
                deallocate(work)
            end block
            call target_dtpmlqt(side, transes(jt), np, nrhs, mp, 0, mb, V, mp, &
                                T_arr, mb, Aap_got, mp, Bap_got, np, info)
            err = max(max_rel_err_mat(Aap_got, Aap_ref), max_rel_err_mat(Bap_got, Bap_ref))
            tol = 16.0_ep * real(max(mp, np, nrhs), ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'trans=', transes(jt), ',mp=', mp, ',np=', np
            call report_case(trim(label), err, tol)
            deallocate(Aap_ref, Bap_ref, Aap_got, Bap_got)
        end do
        deallocate(A0, B0, V, T_arr, Aap0, Bap0)
    end do
    call report_finalize()
end program test_dtpmlqt
