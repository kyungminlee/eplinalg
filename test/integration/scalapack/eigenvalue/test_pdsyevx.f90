program test_pdsyevx
    ! PDSYEVX: bisection + inverse-iteration symmetric eigensolver.
    ! RANGE='A', JOBZ='N' (eigenvalues only). Compare to dsyev.
    !
    ! Previously blocked by a heap-corruption bug rooted in upstream
    ! PJLAENV's strict S/D/C/Z precision-letter gate; see the
    ! recipes/scalapack/extras/pjlaenv_ep.f override.
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: dsyev
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix
    use target_scalapack,  only: target_name, target_eps, target_pdsyevx
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, liwork, m_got, nz_got
    integer :: locm, locn, lld
    integer :: desca(9), descZ(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_sym(:,:), A_ref(:,:)
    real(ep), allocatable :: Z_loc(:,:)
    real(ep), allocatable :: w(:), w_ref(:), gap(:), work(:), work_ref(:)
    integer,  allocatable :: iwork(:), ifail(:), iclustr(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdsyevx', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 35101 + 31*i)

        allocate(A_sym(n, n))
        A_sym = 0.5_ep * (A_glob + transpose(A_glob))
        if (size(A_loc, 1) > 0 .and. size(A_loc, 2) > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) A_loc(il, jl) = A_sym(ig, jg)
                end do
            end do
        end if

        locm = numroc_local(n, mb, my_row, 0, my_nprow)
        locn = numroc_local(n, nb, my_col, 0, my_npcol); lld = max(1, locm)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld, info)
        call descinit_local(descZ, n, n, mb, nb, 0, 0, my_context, lld, info)

        allocate(Z_loc(lld, max(1, locn))); Z_loc = 0.0_ep
        allocate(w(n), gap(my_nprow * my_npcol), &
                 ifail(n), iclustr(2 * my_nprow * my_npcol))

        allocate(work(1), iwork(1))
        call target_pdsyevx('N', 'A', 'U', n, A_loc, 1, 1, desca, &
                            0.0_ep, 0.0_ep, 0, 0, -1.0_ep, m_got, nz_got, &
                            w, 1.0e-3_ep, Z_loc, 1, 1, descZ, &
                            work, -1, iwork, -1, ifail, iclustr, gap, info)
        lwork  = max(1, int(work(1)))
        liwork = max(1, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call target_pdsyevx('N', 'A', 'U', n, A_loc, 1, 1, desca, &
                            0.0_ep, 0.0_ep, 0, 0, -1.0_ep, m_got, nz_got, &
                            w, 1.0e-3_ep, Z_loc, 1, 1, descZ, &
                            work, lwork, iwork, liwork, ifail, iclustr, gap, info)

        if (my_rank == 0) then
            allocate(A_ref(n, n), w_ref(n), work_ref(max(1, 64 * n)))
            A_ref = A_sym
            call dsyev('N', 'U', n, A_ref, n, w_ref, &
                       work_ref, size(work_ref), info_ref)
            ! Both produce eigenvalues in ascending order.
            err = max_rel_err_vec(w(1:m_got), w_ref(1:m_got))
            tol = 64.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',m=', m_got
            call report_case(trim(label), err, tol)
            deallocate(A_ref, w_ref, work_ref)
        end if
        deallocate(A_loc, A_glob, A_sym, Z_loc, w, gap, &
                   ifail, iclustr, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdsyevx
