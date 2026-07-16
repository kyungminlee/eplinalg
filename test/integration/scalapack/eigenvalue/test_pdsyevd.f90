program test_pdsyevd
    ! Symmetric eigensolver via divide-and-conquer; computes eigenvalues
    ! and (always) eigenvectors. Reference path is dsyevd at quad
    ! precision.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dsyevd
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdsyevd
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, liwork, j
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9), descz(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    real(ep), allocatable :: A_loc(:,:), Z_loc(:,:), Z_glob(:,:)
    real(ep), allocatable :: A_glob(:,:), A_sym(:,:), A_ref(:,:)
    real(ep), allocatable :: w_got(:), w_ref(:), work(:), work_ref(:)
    integer,  allocatable :: iwork(:), iwork_ref(:)
    real(ep), allocatable :: AZ(:,:), ZD(:,:)
    real(ep) :: err, tol, anorm, rnorm
    character(len=1) :: jobz

    call grid_init()
    call report_init('pdsyevd', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i); jobz = 'V'
        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 13001 + 31*i)

        allocate(A_sym(n, n))
        A_sym = 0.5_ep * (A_glob + transpose(A_glob))

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        if (locm_a > 0 .and. locn_a > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) A_loc(il, jl) = A_sym(ig, jg)
                end do
            end do
        end if
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descz, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(Z_loc(max(1, locm_a), max(1, locn_a)))
        Z_loc = 0.0_ep
        allocate(w_got(n), work(1), iwork(1))
        call target_pdsyevd(jobz, uplos(i), n, A_loc, 1, 1, desca, w_got, &
                            Z_loc, 1, 1, descz, work, -1, iwork, -1, info)
        lwork  = max(1, int(work(1)))
        liwork = max(1, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call target_pdsyevd(jobz, uplos(i), n, A_loc, 1, 1, desca, w_got, &
                            Z_loc, 1, 1, descz, work, lwork, &
                            iwork, liwork, info)
        call gather_matrix(n, n, mb, nb, Z_loc, Z_glob)

        if (my_rank == 0) then
            ! dsyevd JOBZ='V' wants LWORK >= 1 + 6*N + 2*N**2; LIWORK >=
            ! 3 + 5*N. The previous "n*64" sizing was tight enough that
            ! n=32 came up short and dsyevd returned info=-8, so the
            ! report_case calls below never fired and the JSON came back
            ! header-only.
            allocate(A_ref(n, n), w_ref(n), &
                     work_ref(max(1, 1 + 6 * n + 2 * n * n)), &
                     iwork_ref(max(1, 3 + 5 * n)))
            A_ref = A_sym
            call dsyevd(jobz, uplos(i), n, A_ref, n, w_ref, &
                        work_ref, size(work_ref), iwork_ref, size(iwork_ref), info_ref)
            err = max_rel_err_vec(w_got, w_ref)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            block
                character(len=48) :: label
                write(label, '(a,a,a,i0,a)') 'uplo=', uplos(i), ',n=', n, ',out=W'
                call report_case(trim(label), err, tol)
            end block
            ! Eigenvector residual.
            allocate(AZ(n, n), ZD(n, n))
            AZ = matmul(A_sym, Z_glob)
            do j = 1, n
                ZD(:, j) = Z_glob(:, j) * w_got(j)
            end do
            anorm = maxval(abs(A_sym))
            rnorm = maxval(abs(AZ - ZD))
            err = rnorm / max(anorm, tiny(1.0_ep))
            block
                character(len=48) :: label
                write(label, '(a,a,a,i0,a)') 'uplo=', uplos(i), ',n=', n, ',out=residual'
                call report_case(trim(label), err, tol)
            end block
            deallocate(AZ, ZD)
            deallocate(A_ref, w_ref, work_ref, iwork_ref)
        end if
        if (allocated(Z_glob)) deallocate(Z_glob)
        deallocate(A_loc, Z_loc, A_glob, A_sym, w_got, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdsyevd
